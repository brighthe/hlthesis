import numpy as np

from shape_gradient import TopLsfShapeGrad

# Cantilever 的默认参数
domain_width = 3
domain_hight = 2
nelx = 3
nely = 2
ew = domain_width / nelx
eh = domain_hight / nely
volReq = 0.5
stepLength = 2;
numReinit = 3
ts = TopLsfShapeGrad(domain_width = domain_width, domain_hight = domain_hight, \
                     nelx = nelx, nely = nely, \
                     volReq=volReq, stepLength=stepLength, numReinit=3)

fe_domain = [0, domain_width, 0, domain_hight]
fe_mesh = ts.generate_mesh(domain=fe_domain, nelx=nelx, nely=nely)
fe_NC = fe_mesh.number_of_cells()
fe_node = fe_mesh.entity('node') # 左下角按列增加
fe_cell = fe_mesh.entity('cell') # 左下角逆时针
fe_center_node = fe_mesh.entity_barycenter('cell')
fe_center_x = fe_center_node[:, 0]
fe_center_y = fe_center_node[:, 1]
fe_x = fe_node[:, 0]
fe_y = fe_node[:, 1]
#print("fe_x:", fe_x.shape, "\n", fe_x.round(4))
#print("fe_y:", fe_y.shape, "\n", fe_y.round(4))

#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#fe_mesh.add_plot(axes)
#fe_mesh.find_node(axes, showindex=True, markersize=6, fontsize=6, fontcolor='r')
#fe_mesh.find_cell(axes, showindex=True, markersize=8, fontsize=8, fontcolor='k')
#plt.show()

ls_domain = [-ew/2, domain_width+ew/2, -eh/2, domain_hight+eh/2]
ls_mesh = ts.generate_mesh(domain=ls_domain, nelx=nelx+1, nely=nely+1)
ls_NC = ls_mesh.number_of_cells()
ls_NN = ls_mesh.number_of_nodes()
ls_node = ls_mesh.entity('node') # 左下角按列增加
ls_cell = ls_mesh.entity('cell') # 左下角逆时针
ls_center_node = ls_mesh.entity_barycenter('cell')
ls_x = ls_node[:, 0]
ls_y = ls_node[:, 1]
#print("ls_x:", ls_x.shape, "\n", ls_x.round(4))
#print("ls_y:", ls_y.shape, "\n", ls_y.round(4))

ls_Phi = ts.init_lsf(mesh = ls_mesh)
print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi.round(4))

#ts.plot(x=ls_x, y=ls_y, z=ls_Phi, label='ls_Phi')

# 边界条件处理
boundary_condition = (ls_x - np.min(ls_x)) * (ls_x - np.max(ls_x)) * \
                     (ls_y - np.max(ls_y)) * (ls_y - np.min(ls_y)) \
                        <= 100 * np.finfo(float).eps
#print("boundary_condition:", boundary_condition)
ls_Phi[boundary_condition] = -1e-6
print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi.round(4))

# 导入 matlab 的数据
from scipy.io import loadmat
mat = loadmat('fe_phi_min.mat')
data = mat['FENd'][0, 0]
phi_data = data['Phi'].astype(np.float64)
fe_Phi = phi_data[:, 0]
print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))

#from scipy.io import loadmat
#mat = loadmat('fe_phi_max.mat')
#data = mat['FENd'][0, 0]
#phi_data = data['Phi'].astype(np.float64)
#fe_Phi = phi_data[:, 0]
#print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))

distances = (ls_x[None, :] - fe_center_x[:, None])**2 + \
            (ls_y[None, :] - fe_center_y[:, None])**2
print("distances:", distances.shape, "\n", distances)

# 有限元网格单元中心对应的水平集网格节点的索引
ele_lsgrid_id = np.argmin(distances, axis=1)
print("ele_lsgrid_id:", ele_lsgrid_id.shape, "\n", ele_lsgrid_id)
## 水平集函数值 Phi 从水平集节点投影到有限元节点
#from scipy.interpolate import griddata
#fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')
#print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))

#ts.plot_mesh(x0=ls_x, y0=ls_y, label0='ls_grid', x=fe_x, y=fe_y, z=fe_Phi, label='fe_Phi')


from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(fe_mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 定义荷载
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
tmp = vgdof - (nely + 1)
print("tmp:", tmp)
F[tmp, 0] = -1
print("F:", F.shape, "\n", F)

# 位移约束(supports) - short cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)
print("fixeddofs:", fixeddofs)

E0 = 1e-3
E1 = 1.0
nu = 0.3
lagV = 10;
lagCur = 0.1;
fea_Intercal = 10;

totalNum = 1
# 开始循环
for iterNum in range(totalNum):
    # 有限元分析
    print(f'Finite Element Analysis No.: {iterNum}')

    #fe_Phi = np.array([1, 2, 3, 4, 5, 0, 0, -0.5, -1, -2, -3, -4])
    #print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))
    U = ts.fe_analysis(mesh=fe_mesh, E0=E0, E1=E1, nu=nu, ew=ew, eh=eh, Phi=fe_Phi,
                       F=F, fixeddofs=fixeddofs)
    print("U:", U.shape, "\n", U)
    ux = U[:, 0]
    uy = U[:, 1]

    print("F[tmp]:", F[tmp])
    mean_compliances = F[tmp] * U.reshape(-1, 1)[tmp]
    print("mean_compliances:", mean_compliances)

    phi = ls_Phi.reshape(nelx+2, nely+2).T
    ls_curv = ts.calc_curvature(phi=phi, dx=ew, dy=eh)

    ls_Beta = ts.sensi_analysis(fe_mesh=fe_mesh, ls_mesh=ls_mesh, E1=E1, E0=E0, \
                                u=ux, v=uy, ew=ew, eh=eh, nu=nu, lag4Vol=lagV, lag4Curv=lagCur, \
                                ele_lsgrid_id=ele_lsgrid_id, phi=ls_Phi, curvature=ls_curv)

    # 得到速度场
    ls_Vn = ls_Beta / np.max(np.abs(ls_Beta))
    print("ls_Vn:", ls_Vn.shape, "\n", ls_Vn)

    # 水平集界面更新
    phi0 = ls_Phi.reshape(nelx+2, nely+2).T
    vn = ls_Vn.reshape(nelx+2, nely+2).T
    ls_Phi = ts.level_set_evolve(phi0=phi0, vn=vn, dx=ew, dy=eh, loop_num=fea_Intercal)

    # 水平集界面重置化
    if (iterNum == 0) or ((iterNum+1) % 5 == 0):
        phi0 = ls_Phi.reshape(nelx+2, nely+2).T
        ls_Phi = ts.reinitialize(phi0=phi0, dx=ew, dy=eh, loop_num=20)

    asd




#import matplotlib.pyplot as plt
#plt.figure(figsize=(8, 8))
## 绘制有限元网格
#for cell in fe_cell:
#    plt.fill(fe_node[cell, 0], fe_node[cell, 1], edgecolor='red', fill=False, linewidth=1.5)
#plt.plot(fe_node[:, 0], fe_node[:, 1], 's', color='red', label='FE Mesh Nodes')
## 绘制水平集网格
#for cell in ls_cell:
#    plt.fill(ls_node[cell, 0], ls_node[cell, 1], edgecolor='blue', fill=False, linewidth=1.5)
#plt.plot(ls_node[:, 0], ls_node[:, 1], 'o', color='blue', label='LS Mesh Nodes')
#plt.legend()
#plt.title('Combined Visualization of FE and LS Meshes')
#plt.axis('equal')
#plt.show()


## 定义 17 个孔洞的圆心位置
#cx = domain_width/200 * np.array([33.33,  100,  166.67,   0,    66.67, 133.33,
#                                  200,  33.33,  100,   166.67,   0,   66.67,
#                                 133.33, 200,  33.33,   100,   166.67], dtype=np.float64)
#cy = domain_hight/100 * np.array([0,  0,  0,   25,  25, 25,
#                                 25, 50, 50,  50,  75, 75,
#                                 75, 75, 100, 100, 100], dtype=np.float64)
#
## 计算每个水平集网格点与所有圆心的距离，创建初始界面
#tmpPhi = -np.sqrt((ls_x[:, None] - cx)**2 + (ls_y[:, None] - cy)**2) + domain_hight / 10
#print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))
#
## 取这些距离的最大负值来定义初始的水平集函数
#ls_Phi = -np.max(tmpPhi, axis=1)
#print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi.round(4))




#def output(filename, output, mesh, node_val=None, cell_val=None):
#    import os
#    if not os.path.exists(output):
#        os.makedirs(output)
#    if node_val is not None:
#        mesh.nodedata['node_val'] = node_val.flatten('F') # 按列增加
#        fname = os.path.join(output, f'{filename}.vtu')
#        mesh.to_vtk(fname=fname)
#    if cell_val is not None:
#        pass
#
#output(filename='ls_Phi', output='./visulaization/', mesh=ls_mesh, node_val=ls_Phi)
#output(filename='fe_Phi', output='./visulaization/', mesh=fe_mesh, node_val=fe_Phi)





#import matplotlib.pyplot as plt
#
#plt.figure(figsize=(8, 6))
#
## 绘制原始水平集网格点及Phi值
#plt.scatter(ls_x, ls_y, c=ls_Phi, cmap='viridis', label='LSgrid Phi', edgecolors='k', s=100)
#X, Y = np.meshgrid(ls_x, ls_y)
#print(X.shape)
#Z = ls_Phi.reshape(len(ls_y), len(ls_x))
#plt.contour(X, Y, Z, levels=[0], colors='red')  # 绘制零水平集边界
### 绘制有限元网格点
##plt.scatter(fe_x, fe_y, color='red', marker='s', label='FENd Points')
#
## 绘制有限元网格点及Phi值
##plt.scatter(fe_x, fe_y, c=fe_Phi, cmap='cool', label='FENd Phi', edgecolors='r', marker='s', s=100)
### 绘制水平集网格点
##plt.scatter(ls_x, ls_y, color='blue', marker='o', label='LSgrid Points')
#
#plt.colorbar(label='Phi values')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Visualization of Phi Mapping from LSgrid to FENd')
#plt.legend()
#plt.grid(True)
#plt.show()
