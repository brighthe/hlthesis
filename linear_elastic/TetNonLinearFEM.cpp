
#include "TetNonLinearFEM.h"
#include <cmath>
#include <iostream>

using namespace Eigen;
namespace TOPT {
namespace FEM {

void TetNonLinearFEM::init_data()
{
    // 初始化点元数组和体元数组的形状
    //  * 点元数组（NN, 3）
    //  * 体元数组（NC, 4）
    // 注意，这里假设 m_nodes 和 m_cells 中已经填入网格数据
    m_nodes.init_shape(); 
    m_cells.init_shape();

    m_gflag = false; // 没有重力
    m_d = {0.0, 0.0, 0.0};
    m_g = 0.0;

    auto NN = number_of_nodes();
    auto NC = number_of_cells();

    // 初始化数据内存，
    m_u.init(NN, 3); // 点元上的位移
    m_isddof.init(NN, 3);  // Dirichlet 位移边界的标记

    // 体元应变数据，这里做了一维展开 
    // [epsilon_XX, epsilon_YY, epsilon_ZZ, epsilon_YZ, epsilon_XZ, epsilon_XY]
    m_strain.init(NC, 6); 

    // 体元应力数据，这里做了一维展开 
    // [sigma_XX, sigma_YY, sigma_ZZ, sigma_YZ, sigma_XZ, sigma_XY]
    m_stress.init(NC, 6);

    // 体元主轴应力数据 
    // [sigma_0, sigma_1, sigma_2]
    // sigma_0 < sigma_1 < sigma_2
    m_pstress.init(NC, 3);

    m_gphi.init(NC, 4, 3);  //体元线性基函数导数
    m_vols.resize(NC, 0.0); //所有网格体元的体积

    // 初始化体元的体积
    init_cell_measure();

    // 初始化体元上的线性基函数的导数
    init_grad_basis();

    // 初始化网格拓扑
    init_mesh_top();
}

TetNonLinearFEM::Face TetNonLinearFEM::cell_local_sorted_face(size_type c, unsigned i)
{
  Face f = 
  {
      m_cells(c, m_localface[i][0]), 
      m_cells(c, m_localface[i][1]), 
      m_cells(c, m_localface[i][2])
  };
  // 从小到大排序， 这样共享同一个三角形面两个四面体单元得到是相同的数组
  std::sort(f.begin(), f.end()); 
  return  f;
}

void TetNonLinearFEM::init_mesh_top()
{
    auto NN = m_nodes.size();
    auto NC = m_cells.size();
    m_face2cell.reserve(2*NC); // 2*NF = 4*NC
    std::map<Face, unsigned> fidxmap; // sorted face 到整体编号的映射

    unsigned NF = 0;
    // 遍历所有体元
    for(unsigned c = 0; c < NC; c++)
    {
      // 遍历 4 个三角形面元
      for(unsigned i = 0; i < 4; i++)
      {
         auto f = cell_local_sorted_face(c, i); // 第 c 个体元的第 i 个 sorted face
         auto it = fidxmap.find(f);
         if(it == fidxmap.end())
         {
            fidxmap.insert(std::pair<Face, unsigned>(f, NF));
            m_face2cell.push_back({c, c, i, i});
            NF++;
         }
         else
         {
            m_face2cell[it->second][1] = c;
            m_face2cell[it->second][3] = i;
         }
     }
    }
    fidxmap.clear();

    m_faces.reserve(NF);
    for(unsigned f = 0; f < NF; f++)
    {
        auto i = m_face2cell[f][2];
        m_faces.push_back({
            m_cells(m_face2cell[f][0], m_localface[i][0]), 
            m_cells(m_face2cell[f][0], m_localface[i][1]),
            m_cells(m_face2cell[f][0], m_localface[i][2])});
    }
}

template<typename Writer>
void TetNonLinearFEM::output(Writer & writer, std::string &fname)
{
  auto NN = m_nodes.size();
  writer.set_mesh(m_nodes, m_cells, 10); // 10 表示四面体元
  writer.set_point_data(m_u.data(), 3, "displacement");
  writer.set_point_data(m_isddof.data(), 3, "isDDof");

  auto s11 = cell_stress(0); // XX
  auto s22 = cell_stress(1); // YY
  auto s33 = cell_stress(2); // ZZ
  auto s23 = cell_stress(3); // YZ
  auto s13 = cell_stress(4); // XZ
  auto s12 = cell_stress(5); // XY 

  writer.set_cell_data(s11, 1, "s11");
  writer.set_cell_data(s22, 1, "s22");
  writer.set_cell_data(s33, 1, "s33");
  writer.set_cell_data(s23, 1, "s23");
  writer.set_cell_data(s13, 1, "s13");
  writer.set_cell_data(s12, 1, "s12");

  auto S11 = cell_data_to_node_data(s11, 1);
  auto S22 = cell_data_to_node_data(s22, 1);
  auto S33 = cell_data_to_node_data(s33, 1);
  auto S23 = cell_data_to_node_data(s23, 1);
  auto S13 = cell_data_to_node_data(s13, 1);
  auto S12 = cell_data_to_node_data(s12, 1);

  writer.set_point_data(S11, 1, "S11");
  writer.set_point_data(S22, 1, "S22");
  writer.set_point_data(S33, 1, "S33");
  writer.set_point_data(S23, 1, "S23");
  writer.set_point_data(S13, 1, "S13");
  writer.set_point_data(S12, 1, "S12");

  auto ms = cell_mises_stress();
  auto Mises = cell_data_to_node_data(ms, 1);
  writer.set_point_data(Mises, 1, "Mises");

  auto ts = cell_tresca_stress();
  auto Tresca = cell_data_to_node_data(ts, 1);
  writer.set_point_data(Tresca, 1, "Tresca");

  writer.write(fname);
}

void TetNonLinearFEM::set_disp_condtion(std::vector<unsigned> nidxs, double val)
{
  for(auto i : nidxs)
  {
    m_u(i, 0) = val;
    m_u(i, 1) = val;
    m_u(i, 2) = val;

    m_isddof(i, 0) = true;
    m_isddof(i, 1) = true;
    m_isddof(i, 2) = true;
  }
}


std::vector<double> TetNonLinearFEM::cell_mises_stress()
{
  auto NC = number_of_cells();
  std::vector<double> s(NC, 0.0);
  for(auto i = 0; i < NC; i++)
  {
    auto val = m_stress[i][0] - m_stress[i][1];
    s[i] += val*val;

    val = m_stress[i][0] - m_stress[i][2];
    s[i] += val*val;

    val = m_stress[i][1] - m_stress[i][2];
    s[i] += val*val;

    s[i] += 6*(m_stress[i][3]*m_stress[i][3] + m_stress[i][4]*m_stress[i][4] + m_stress[i][5]*m_stress[i][5]);   

    s[i] = std::sqrt(s[i]/2.0);
  }
  return s;
}

std::vector<double> TetNonLinearFEM::cell_tresca_stress()
{
  auto NC = number_of_cells();
  std::vector<double> s(NC, 0.0);
  for(auto i = 0; i < NC; i++)
  {
    s[i] = m_pstress[i][2] - m_pstress[i][0];
  }
  return s;
}

std::vector<double> TetNonLinearFEM::cell_max_abs_principal_stress()
{
  auto NC = number_of_cells();
  std::vector<double> s(NC, 0.0);
  for(auto i = 0; i < NC; i++)
  {
    s[i] = std::abs(m_pstress[i][0]);
    auto val = std::abs(m_pstress[i][1]);
    if(s[i] < val)
      s[i] = val;
    val = std::abs(m_pstress[i][2]);
    if(s[i] < val)
      s[i] = val;
  }
  return s;
}

std::vector<double> TetNonLinearFEM::cell_data_to_node_data(std::vector<double> & data, unsigned nc)
{
  auto NN = number_of_nodes();
  auto NC = number_of_cells();
  std::vector<double> nd(nc*NN, 0.0);
  std::vector<unsigned> nn(NN, 0);

  for(auto i = 0; i < NC; i++)
  {
    for(auto j = 0; j < 4; j++)
    {
      auto m = m_cells[i][j];
      for(auto k = 0; k < nc; k++)
      {
        nd[nc*m + k] += data[nc*i + k];
      }
      nn[m] += 1;
    }
  }

  for(auto i = 0; i < NN; i++)
  {
    if(nn[i] > 0)
    {
      for(auto j = 0; j < nc; j++)
      {
        nd[nc*i + j] /= nn[i];
      }
    }
  }
  return nd;
}

void TetNonLinearFEM::compute_strain()
{
  auto NC = number_of_cells();
  for(auto i = 0; i < NC; i++)
  {
    m_strain[i][0] = 0.0;
    m_strain[i][1] = 0.0;
    m_strain[i][2] = 0.0;
    m_strain[i][3] = 0.0;
    m_strain[i][4] = 0.0;
    m_strain[i][5] = 0.0;

    double lambda;
    double mu;
    compute_cell_elastic_modulus(i, lambda, mu);

    // 计算应变
    for(auto j = 0; j < 4; j++)
    {
      m_strain[i][0] += m_u(m_cells[i][j], 0)*m_gphi(i, j, 0);
      m_strain[i][1] += m_u(m_cells[i][j], 1)*m_gphi(i, j, 1);
      m_strain[i][2] += m_u(m_cells[i][j], 2)*m_gphi(i, j, 2);

      m_strain[i][3] += m_u(m_cells[i][j], 2)*m_gphi(i, j, 1); // (w_y + v_z)/2
      m_strain[i][3] += m_u(m_cells[i][j], 1)*m_gphi(i, j, 2);

      m_strain[i][4] += m_u(m_cells[i][j], 2)*m_gphi(i, j, 0);
      m_strain[i][4] += m_u(m_cells[i][j], 0)*m_gphi(i, j, 2);

      m_strain[i][5] += m_u(m_cells[i][j], 1)*m_gphi(i, j, 0);
      m_strain[i][5] += m_u(m_cells[i][j], 0)*m_gphi(i, j, 1);
    }

    m_strain[i][3] /= 2.0;
    m_strain[i][4] /= 2.0;
    m_strain[i][5] /= 2.0;
  }
}
  
  
double TetNonLinearFEM::compute_young_modulus(double strain)
{
  if(strain < 1e-10)
    return m_ustress[1]/m_ustrain[1];
  
  int k = 0;
  for( ; k < m_ustrain.size()-1; k++)
  {
    if(strain < m_ustrain[k+1])
    {
      break;
    }
  }

  if(k == m_ustrain.size()-1)
  {
    return m_ustress[k]/strain;
  }
  else
  {
    double l = (strain-m_ustrain[k])/(m_ustrain[k+1]-m_ustrain[k]);
    double stress = (1-l)*m_ustrain[k] + l*m_ustrain[k+1];      
    return stress/strain;
  }
}


void TetNonLinearFEM::compute_cell_elastic_modulus(unsigned i, double & lambda, double & mu)
{
  double s11 = m_strain[i][0];
  double s22 = m_strain[i][1];
  double s33 = m_strain[i][2];
  double s23 = m_strain[i][3];
  double s13 = m_strain[i][4];
  double s12 = m_strain[i][5];
  double s1 = s11 + s22 + s33;
  double s = s1/3;
  s11 -= s;
  s22 -= s;
  s33 -= s;
  s = s11*s11 + s22*s22 + s33*s33 + 2*(s23*s23 + s13*s13 + s12*s12);
  double gamma = std::sqrt(2*s);
  if(m_type)
  {
    double epsilon = std::sqrt(3)*gamma/(1+m_nu)/2.0;
    double E = compute_young_modulus(epsilon);
    lambda = m_nu*E/(1+m_nu)/(1-2*m_nu);
    mu = E/(1+m_nu)/2.0;
  }
  else
  {
    double epsilon = gamma/std::sqrt(3)+std::abs(s1)/3.0;
    double E = compute_young_modulus(epsilon);
    mu = 3*m_K*E/(9*m_K-E);
    lambda = m_K-mu*2/3.0;
  }
}

void TetNonLinearFEM::compute_stress()
{
  auto NC = number_of_cells();
  for(auto i = 0; i < NC; i++)
  {
    //TODO: 计算当前单元上的 mu 和 lambda
    // 计算应力
    double lambda;
    double mu;
    compute_cell_elastic_modulus(i, lambda, mu);

    auto c = 2*mu + lambda;
    m_stress[i][0] = c*m_strain[i][0] + lambda*(m_strain[i][1] + m_strain[i][2]); 
    m_stress[i][1] = c*m_strain[i][1] + lambda*(m_strain[i][2] + m_strain[i][0]);
    m_stress[i][2] = c*m_strain[i][2] + lambda*(m_strain[i][0] + m_strain[i][1]);
    m_stress[i][3] = 2*mu*m_strain[i][3];
    m_stress[i][4] = 2*mu*m_strain[i][4];
    m_stress[i][5] = 2*mu*m_strain[i][5];

    Matrix3d A;
    A(0, 0) = m_stress[i][0];
    A(0, 1) = m_stress[i][5];
    A(0, 2) = m_stress[i][4];

    A(1, 0) = m_stress[i][5];
    A(1, 1) = m_stress[i][1];
    A(1, 2) = m_stress[i][3];

    A(2, 0) = m_stress[i][4]; 
    A(2, 1) = m_stress[i][3];
    A(2, 2) = m_stress[i][2];

    EigenSolver es(A);
    auto v = es.eigenvalues(); 
    m_pstress[i][0] = v[0];
    m_pstress[i][1] = v[1];
    m_pstress[i][2] = v[2];
  }
}

void TetNonLinearFEM::compute_strain_and_stress()
{
  compute_strain();
  compute_stress();
}

void TetNonLinearFEM::apply_load_boundary_condition(VectorXd & b)
{
  auto NN = m_nodes.size();
  
  // 施加以点元全局编号设定的载荷
  for(auto & load: m_nloads)
  {
    for(auto n: load.nidxs)
    {
      b[3*n    ] += load.value[0];
      b[3*n + 1] += load.value[1];
      b[3*n + 2] += load.value[2];
    }
  }


  for(auto & load: m_celoads)
  {
    auto N = load.ceidxs.size()/2;
    for(auto i = 0; i < N; i++)
    {
      auto c = load.ceidxs[2*i];
      auto e = load.ceidxs[2*i+1];

      auto v0 = m_cells[c][m_localedge[e][0]];
      auto v1 = m_cells[c][m_localedge[e][1]];

      auto x0 = m_nodes[v1][0] - m_nodes[v0][0];
      auto y0 = m_nodes[v1][1] - m_nodes[v0][1];
      auto z0 = m_nodes[v1][2] - m_nodes[v0][2];

      auto l = std::sqrt(x0*x0 + y0*y0 + z0*z0);
      Vector force = {
        l*load.value[0]/2.0, 
        l*load.value[1]/2.0,
        l*load.value[2]/2.0
      };

      b[3*v0    ] += force[0];
      b[3*v0 + 1] += force[1];
      b[3*v0 + 2] += force[2];

      b[3*v1    ] += force[0];
      b[3*v1 + 1] += force[1];
      b[3*v1 + 2] += force[2];
    }
  }

  // 施加以体元全局编号及其局部面编号的载荷条件
  for(auto & load: m_cfloads)
  {
    auto N = load.cfidxs.size()/2;
    for(auto i = 0; i < N; i++)
    {
      auto c = load.cfidxs[2*i];
      auto f = load.cfidxs[2*i+1];

      auto v0 = m_cells[c][m_localface[f][0]];
      auto v1 = m_cells[c][m_localface[f][1]];
      auto v2 = m_cells[c][m_localface[f][2]];

      auto x1 = m_nodes[v1][0] - m_nodes[v0][0];
      auto y1 = m_nodes[v1][1] - m_nodes[v0][1];
      auto z1 = m_nodes[v1][2] - m_nodes[v0][2];

      auto x2 = m_nodes[v2][0] - m_nodes[v0][0];
      auto y2 = m_nodes[v2][1] - m_nodes[v0][1];
      auto z2 = m_nodes[v2][2] - m_nodes[v0][2];

      Vector force = {
        (z1*y2 - y1*z2)*load.value[0]/6.0, 
        (x1*z2 - z1*x2)*load.value[0]/6.0,
        (y1*x2 - x1*y2)*load.value[0]/6.0
      };

      b[3*v0    ] += force[0];
      b[3*v0 + 1] += force[1];
      b[3*v0 + 2] += force[2];

      b[3*v1    ] += force[0];
      b[3*v1 + 1] += force[1];
      b[3*v1 + 2] += force[2];

      b[3*v2    ] += force[0];
      b[3*v2 + 1] += force[1];
      b[3*v2 + 2] += force[2];
    }
  }

  // 施加以面元全局编号设定的载荷
  for(auto & load:m_floads)
  {
    for(auto f:load.fidxs)
    {
      auto x1 = m_nodes[m_faces[f][1]][0] - m_nodes[m_faces[f][0]][0];
      auto y1 = m_nodes[m_faces[f][1]][1] - m_nodes[m_faces[f][0]][1];
      auto z1 = m_nodes[m_faces[f][1]][2] - m_nodes[m_faces[f][0]][2];

      auto x2 = m_nodes[m_faces[f][2]][0] - m_nodes[m_faces[f][0]][0];
      auto y2 = m_nodes[m_faces[f][2]][1] - m_nodes[m_faces[f][0]][1];
      auto z2 = m_nodes[m_faces[f][2]][2] - m_nodes[m_faces[f][0]][2];

      Vector force = {
        (z1*y2 - y1*z2)*load.value[0]/6.0, 
        (x1*z2 - z1*x2)*load.value[0]/6.0,
        (y1*x2 - x1*y2)*load.value[0]/6.0
      };

      std::cout << force[0] << ", " << force[1] << ", " << force[2] << std::endl;

      for(unsigned i = 0; i < 3; i++)
      {
        b[3*m_faces[f][i]    ] += force[0];
        b[3*m_faces[f][i] + 1] += force[1];
        b[3*m_faces[f][i] + 2] += force[2];
      }
    }
  }
}


void TetNonLinearFEM::apply_disp_boundary_condition(VectorXd & x, CSRMatrix & K,  VectorXd & b)
{

  auto gdof = 3*m_nodes.size();

  auto data = K.valuePtr(); // 非零元数组
  auto indices = K.innerIndexPtr(); // 非零元对应的列指标数组
  auto indptr = K.outerIndexPtr();  // 非零元的起始位置数组

  b -= K*x;
  for(unsigned i = 0; i < gdof; i++)
  {
    if(m_isddof.data()[i])
    {
      b[i] = m_u.data()[i];
      for(auto k = indptr[i]; k < indptr[i+1]; k++)
      {
        auto j = indices[k];
        if( i == j)
        {
          data[k] = 1.0;
        }
        else
        {
          data[k] = 0.0;
        }
      }
    }
    else
    {
      for(auto k = indptr[i]; k < indptr[i+1]; k++)
      {
        auto j = indices[k];
        if(m_isddof.data()[j])
          data[k] = 0.0; 
      }
    }
  }
}

double TetNonLinearFEM::cell_measure(unsigned c)
{
  auto x1 = m_nodes[m_cells[c][1]][0] - m_nodes[m_cells[c][0]][0];
  auto y1 = m_nodes[m_cells[c][1]][1] - m_nodes[m_cells[c][0]][1];
  auto z1 = m_nodes[m_cells[c][1]][2] - m_nodes[m_cells[c][0]][2];

  auto x2 = m_nodes[m_cells[c][2]][0] - m_nodes[m_cells[c][0]][0];
  auto y2 = m_nodes[m_cells[c][2]][1] - m_nodes[m_cells[c][0]][1];
  auto z2 = m_nodes[m_cells[c][2]][2] - m_nodes[m_cells[c][0]][2];

  auto x3 = m_nodes[m_cells[c][3]][0] - m_nodes[m_cells[c][0]][0];
  auto y3 = m_nodes[m_cells[c][3]][1] - m_nodes[m_cells[c][0]][1];
  auto z3 = m_nodes[m_cells[c][3]][2] - m_nodes[m_cells[c][0]][2];

  double vol = x3*(y1*z2 - z1*y2);
  vol += y3*(z1*x2 - x1*z2);
  vol += z3*(x1*y2 - y1*x2);
  return vol/6.0;

}

void TetNonLinearFEM::init_cell_measure()
{
  for(unsigned c = 0; c < m_cells.size(); c++)
  {
    m_vols[c] = cell_measure(c);
  }
}

void TetNonLinearFEM::init_grad_basis()
{
  auto NC = m_cells.size();
  for(unsigned c = 0; c < NC; c++)
  {
    auto x12 = m_nodes[m_cells[c][2]][0] - m_nodes[m_cells[c][1]][0];
    auto y12 = m_nodes[m_cells[c][2]][1] - m_nodes[m_cells[c][1]][1];
    auto z12 = m_nodes[m_cells[c][2]][2] - m_nodes[m_cells[c][1]][2];

    auto x13 = m_nodes[m_cells[c][3]][0] - m_nodes[m_cells[c][1]][0];
    auto y13 = m_nodes[m_cells[c][3]][1] - m_nodes[m_cells[c][1]][1]; 
    auto z13 = m_nodes[m_cells[c][3]][2] - m_nodes[m_cells[c][1]][2];

    m_gphi(c, 0, 0) = (z12*y13 - y12*z13)/m_vols[c]/6.0;
    m_gphi(c, 0, 1) = (x12*z13 - z12*x13)/m_vols[c]/6.0;
    m_gphi(c, 0, 2) = (y12*x13 - x12*y13)/m_vols[c]/6.0;

    auto x01 = m_nodes[m_cells[c][1]][0] - m_nodes[m_cells[c][0]][0];
    auto y01 = m_nodes[m_cells[c][1]][1] - m_nodes[m_cells[c][0]][1];
    auto z01 = m_nodes[m_cells[c][1]][2] - m_nodes[m_cells[c][0]][2];

    auto x02 = m_nodes[m_cells[c][2]][0] - m_nodes[m_cells[c][0]][0];
    auto y02 = m_nodes[m_cells[c][2]][1] - m_nodes[m_cells[c][0]][1]; 
    auto z02 = m_nodes[m_cells[c][2]][2] - m_nodes[m_cells[c][0]][2];

    auto x03 = m_nodes[m_cells[c][3]][0] - m_nodes[m_cells[c][0]][0];
    auto y03 = m_nodes[m_cells[c][3]][1] - m_nodes[m_cells[c][0]][1];
    auto z03 = m_nodes[m_cells[c][3]][2] - m_nodes[m_cells[c][0]][2];

    m_gphi(c, 1, 0) = (z03*y02 - y03*z02)/m_vols[c]/6.0;
    m_gphi(c, 1, 1) = (x03*z02 - z03*x02)/m_vols[c]/6.0;
    m_gphi(c, 1, 2) = (y03*x02 - x03*y02)/m_vols[c]/6.0;

    m_gphi(c, 2, 0) = (z01*y03 - y01*z03)/m_vols[c]/6.0;
    m_gphi(c, 2, 1) = (x01*z03 - z01*x03)/m_vols[c]/6.0;
    m_gphi(c, 2, 2) = (y01*x03 - x01*y03)/m_vols[c]/6.0;

    m_gphi(c, 3, 0) = (z02*y01 - y02*z01)/m_vols[c]/6.0;
    m_gphi(c, 3, 1) = (x02*z01 - z02*x01)/m_vols[c]/6.0;
    m_gphi(c, 3, 2) = (y02*x01 - x02*y01)/m_vols[c]/6.0;
  }
}


void TetNonLinearFEM::construct_cell_stiff_matrix(Data3d & CM)
{
  auto NN = m_nodes.size();
  auto NC = m_cells.size();

  CM.init(NC, 12, 12);

  for(unsigned c = 0; c < NC; c++)
  {
    auto vol = m_vols[c];
    double lambda;
    double mu;
    compute_cell_elastic_modulus(c, lambda, mu);
    for(auto m = 0; m < 4; m++)
    {
      for(auto n = 0; n < 4; n++)
      {
        double xx = m_gphi(c, m, 0)*m_gphi(c, n, 0);
        double yy = m_gphi(c, m, 1)*m_gphi(c, n, 1);
        double zz = m_gphi(c, m, 2)*m_gphi(c, n, 2);

        // K00
        CM(c, 3*m, 3*n) = ((2*mu+lambda)*xx + mu*yy + mu*zz)*vol;

        // K11
        CM(c, 3*m+1, 3*n+1) = (mu*xx + (2*mu+lambda)*yy + mu*zz)*vol;

        // K22
        CM(c, 3*m+2, 3*n+2) = (mu*xx + mu*yy + (2*mu+lambda)*zz)*vol;

        // K01 and K10
        CM(c, 3*m, 3*n+1) = (lambda*m_gphi(c, m, 0)*m_gphi(c, n, 1) + mu*m_gphi(c, m, 1)*m_gphi(c, n, 0))*vol;
        CM(c, 3*n+1, 3*m) = CM(c, 3*m, 3*n+1);

        // K12  and K21
        CM(c, 3*m+1, 3*n+2) = (lambda*m_gphi(c, m, 1)*m_gphi(c, n, 2) + mu*m_gphi(c, m, 2)*m_gphi(c, n, 1))*vol;
        CM(c, 3*n+2, 3*m+1) = CM(c, 3*m+1, 3*n+2);

        // K02 and K20
        CM(c, 3*m, 3*n+2) = (lambda*m_gphi(c, m, 0)*m_gphi(c, n, 2) + mu*m_gphi(c, m, 2)*m_gphi(c, n, 0))*vol;
        CM(c, 3*n+2, 3*m) = CM(c, 3*m, 3*n+2);
      }
    }
  }
}


void TetNonLinearFEM::to_coo(Data3d & M, std::vector<Triplet> & tlist)
{
  auto NN = m_nodes.size();
  auto NC = m_cells.size();
  auto & shape = M.shape();
  auto nnz = shape[0]*shape[1]*shape[2];
  tlist.reserve(nnz);
  for(unsigned c = 0; c < NC; c++)
  {
    for(unsigned m = 0; m < 4; m++)
    {
      for(unsigned n = 0; n < 4; n++)
      {
        unsigned i = 3*m_cells[c][m];
        unsigned j = 3*m_cells[c][n];

        tlist.push_back(Triplet(i, j, M(c, 3*m, 3*n)));

        i = 3*m_cells[c][m] + 1;
        j = 3*m_cells[c][n] + 1;
        tlist.push_back(Triplet(i, j, M(c, 3*m+1, 3*n+1)));

        i = 3*m_cells[c][m] + 2;
        j = 3*m_cells[c][n] + 2;
        tlist.push_back(Triplet(i, j, M(c, 3*m+2, 3*n+2)));

        i = 3*m_cells[c][m];
        j = 3*m_cells[c][n] + 1;
        tlist.push_back(Triplet(i, j, M(c, 3*m, 3*n+1)));
        tlist.push_back(Triplet(j, i, M(c, 3*m, 3*n+1)));

        i = 3*m_cells[c][m] + 1;
        j = 3*m_cells[c][n] + 2;
        tlist.push_back(Triplet(i, j, M(c, 3*m+1, 3*n+2)));
        tlist.push_back(Triplet(j, i, M(c, 3*m+1, 3*n+2)));

        i = 3*m_cells[c][m];
        j = 3*m_cells[c][n] + 2;
        tlist.push_back(Triplet(i, j, M(c, 3*m, 3*n+2)));
        tlist.push_back(Triplet(j, i, M(c, 3*m, 3*n+2)));
      }
    }
  }
}

void TetNonLinearFEM::costruct_cell_load_vector(Data2d & b)
{
  auto NN = m_nodes.size();
  auto NC = m_cells.size();

  b.init(NC, 12);

  //重力载荷
  auto f = m_d;
  f[0] *= m_rho*m_g;
  f[1] *= m_rho*m_g;
  f[2] *= m_rho*m_g;

  for(unsigned c = 0; c < NC; c++) 
  {
    double val = m_vols[c]/4.0;
    for(unsigned i = 0; i < 4; i++)
    {
      b(c, 3*i)   = val*f[0];
      b(c, 3*i+1) = val*f[1]; 
      b(c, 3*i+2) = val*f[2]; 
    }
  }
}

void TetNonLinearFEM::construct_load_vector(VectorXd & b)
{
  auto NN = m_nodes.size();
  auto NC = m_cells.size();

  for(int i = 0; i < 3*NN; i++) //这里开始忘记初始化，产生了很奇怪的错误
    b[i] = 0.0;

  if(m_gflag)
  {
    //重力载荷
    std::cout << "处理重力载荷！" << std::endl;
    auto f = m_d;
    f[0] *= m_rho*m_g;
    f[1] *= m_rho*m_g;
    f[2] *= m_rho*m_g;

    for(unsigned c = 0; c < NC; c++)
    {
      auto vol = m_vols[c]/4.0;
      for(unsigned i = 0; i < 3; i++)  
      {
        auto val = vol*f[i];
        b[3*m_cells[c][0] + i] += val; 
        b[3*m_cells[c][1] + i] += val; 
        b[3*m_cells[c][2] + i] += val; 
        b[3*m_cells[c][3] + i] += val; 
      }
    }
  }
}

void TetNonLinearFEM::construct_stiff_matrix(CSRMatrix & K)
{
  Data3d M;
  construct_cell_stiff_matrix(M);
  std::vector<Triplet> tlist;
  to_coo(M, tlist);
  K.setFromTriplets(tlist.begin(), tlist.end());
}

void TetNonLinearFEM::compute_non_linear_elasticity_model(int stype, int para)
{
  auto gdof = 3*m_nodes.size();

  /** 组装载荷向量 */
  std::cout << "组装载荷向量...\n";
  auto start = std::chrono::steady_clock::now();

  VectorXd f(gdof);
  construct_load_vector(f);

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "完成载荷向量组装， 花费时间: "
       << elapsed_seconds.count() << " s\n" << std::endl;
  /** 结束载荷向量组装*/

  /** 处理载荷边界条件 */
  std::cout<< "处理载荷边界条件......\n" <<std::endl;
  apply_load_boundary_condition(f);
  
  for(int j=0; j<para; j++)
  {
    std::cout << "更新次数："
         << j << " \n" << std::endl;

    /** 0. 组装刚度矩阵 */
    std::cout << "组装刚度矩阵...\n";
    start = std::chrono::steady_clock::now();

    CSRMatrix K(gdof, gdof);
    construct_stiff_matrix(K);
    
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "完成刚度矩阵组装，花费时间: "
         << elapsed_seconds.count() << " s\n" << std::endl;
    /** 结束刚度矩阵组装*/

    /** 分次施加载荷条件 */
    double m = double(j+1)/ double(para);
    VectorXd b = f*m;

    /** 处理 位移边界条件 */
    std::cout<< "处理位移边界条件......\n" <<std::endl;
    VectorXd x(gdof);
    for(int i = 0; i < gdof; i++)
    {
      x[i] = m_u.data()[i];
    }
    apply_disp_boundary_condition(x, K, b);

//    JacobiSVD<Eigen::MatrixXd> svd(K);
//    std::cout<<"gdof:"<<gdof<<std::endl;
//    std::cout<<"rank:"<<svd.rank()<<std::endl;

    /** 求解代数系统 */
    start = std::chrono::steady_clock::now();

    if(stype == 0)
    {
      std::cout << "直接法求解代数系统...\n";
      LUSolver solver;
      solver.analyzePattern(K);
      solver.factorize(K); /**< Compute the numerical factorization */
      x = solver.solve(b); /**< Use the factors to solve the linear system */
    }
    else if(stype == 1)
    {
      std::cout << "共轭梯度法求解代数系统...\n";
      CGSolver solver;
      solver.compute(K);
      x = solver.solve(b);
    }

    end = std::chrono::steady_clock::now();
    std::cout<< "残量范数 ||K0*x - b||_2 = " << (K*x - b).norm() << std::endl;
    elapsed_seconds = end - start;
    std::cout << "完成线性代数系统求解，花费时间: "
         << elapsed_seconds.count() << " s\n" << std::endl;

    /** 后处理 */
    for(unsigned i = 0; i < gdof; i++)
    {
      m_u.data()[i] = x[i];
    }

    //计算应变和应力
    compute_strain_and_stress();
  }
}


void TetNonLinearFEM::compute_non_linear_elasticity_model_1(int stype)
{
  auto gdof = 3*m_nodes.size();

  /** 组装载荷向量 */
  std::cout << "组装载荷向量...\n";
  auto start = std::chrono::steady_clock::now();

  VectorXd f(gdof);
  construct_load_vector(f);

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "完成载荷向量组装， 花费时间: "
       << elapsed_seconds.count() << " s\n" << std::endl;
  /** 结束载荷向量组装*/

  /** 处理载荷边界条件 */
  std::cout<< "处理载荷边界条件......\n" <<std::endl;
  apply_load_boundary_condition(f);
 
  int k = 0;
  while(k < 100)
  {
    k += 1;
    /** 0. 组装刚度矩阵 */
    std::cout << "组装刚度矩阵...\n";
    start = std::chrono::steady_clock::now();

    CSRMatrix K(gdof, gdof);
    construct_stiff_matrix(K);
    
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end - start;
    std::cout << "完成刚度矩阵组装，花费时间: "
         << elapsed_seconds.count() << " s\n" << std::endl;
    /** 结束刚度矩阵组装*/

    VectorXd b = f;

    /** 处理 位移边界条件 */
    std::cout<< "处理位移边界条件......\n" <<std::endl;
    VectorXd x(gdof);
    for(int i = 0; i < gdof; i++)
    {
      if(m_isddof.data()[i])
      {
        x[i] = m_u.data()[i];
      }
    }
    apply_disp_boundary_condition(x, K, b);


    /** 求解代数系统 */
    start = std::chrono::steady_clock::now();

    if(stype == 0)
    {
      std::cout << "直接法求解代数系统...\n";
      LUSolver solver;
      solver.analyzePattern(K);
      solver.factorize(K); /**< Compute the numerical factorization */
      x = solver.solve(f); /**< Use the factors to solve the linear system */
    }
    else if(stype == 1)
    {
      std::cout << "共轭梯度法求解代数系统...\n";
      CGSolver solver;
      solver.compute(K);
      x = solver.solve(f);
    }

    end = std::chrono::steady_clock::now();
    std::cout<< "残量范数 ||K0*x - b||_2 = " << (K*x - f).norm()/f.norm() << std::endl;
    elapsed_seconds = end - start;
    std::cout << "完成线性代数系统求解，花费时间: "
         << elapsed_seconds.count() << " s\n" << std::endl;

    double error = 0.0;
    /** 后处理 */
    for(unsigned i = 0; i < gdof; i++)
    {
      double val = x[i] - m_u.data()[i];
      error += val*val; 

      m_u.data()[i] = x[i];
    }

    std::cout << std::sqrt(error) << std::endl;

    //计算应变和应力
    compute_strain_and_stress();
  }
}

} // end of namespace FEM

} // end of namespace TOPT
