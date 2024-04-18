import numpy as np
import os

from wang import TopLSM

# Cantilever 的默认参数
domain_width = 2 # 设计区域的宽度
domain_hight = 1 # 设计区域的高度
nelx = 8 # x 方向的单元数
nely = 4 # y 方向的单元数
lagV = 50 # Lagrange multiplier for volume constraint
lagCur = 0 # Lagrange multiplier for perimeter constraint whose shape sensitivity is curvature
ts = TopLSM(domain_width = domain_width, domain_hight = domain_hight, \
                     nelx = nelx, nely = nely, lagV = lagV, lagCur = lagCur)
# 数据初始化
ew = domain_width / nelx
eh = domain_hight / nely
fe_domain = [0, domain_width, 0, domain_hight]
fe_mesh = ts.generate_mesh(domain=fe_domain, nelx=nelx, nely=nely)
fe_NC = fe_mesh.number_of_cells()
fe_node = fe_mesh.entity('node') # 左下角按列增加
print("fe_NC:", fe_NC)
fe_cell = fe_mesh.entity('cell') # 左下角逆时针
fe_center_node = fe_mesh.entity_barycenter('cell')
fe_center_x = fe_center_node[:, 0]
fe_center_y = fe_center_node[:, 1]
fe_x = fe_node[:, 0]
fe_y = fe_node[:, 1]

ls_domain = [-ew/2, domain_width+ew/2, -eh/2, domain_hight+eh/2]
ls_mesh = ts.generate_mesh(domain=ls_domain, nelx=nelx+1, nely=nely+1)
ls_NC = ls_mesh.number_of_cells()
ls_NN = ls_mesh.number_of_nodes()
ls_node = ls_mesh.entity('node') # 左下角按列增加
ls_cell = ls_mesh.entity('cell') # 左下角逆时针
ls_center_node = ls_mesh.entity_barycenter('cell')
ls_x = ls_node[:, 0]
ls_y = ls_node[:, 1]

distances = (ls_x[None, :] - fe_center_x[:, None])**2 + \
            (ls_y[None, :] - fe_center_y[:, None])**2
#print("distances:", distances.shape, "\n", distances)

# 有限元网格单元中心对应的水平集网格节点的索引
ele_lsgrid_id = np.argmin(distances, axis=1)
#print("ele_lsgrid_id:", ele_lsgrid_id.shape, "\n", ele_lsgrid_id)

# 初始化水平集函数
ls_Phi = ts.init_lsf(mesh = ls_mesh)
print("ls_Phi:", ls_Phi)

ls_mesh.nodedata['ls_phi0'] = ls_Phi.flatten('f')

# 边界条件处理
#boundary_condition = (ls_x - np.min(ls_x)) * (ls_x - np.max(ls_x)) * \
#                     (ls_y - np.max(ls_y)) * (ls_y - np.min(ls_y)) \
#                        <= 100 * np.finfo(float).eps
is_bd_node = ls_mesh.ds.boundary_node_flag()
#print("boundary_condition:", boundary_condition)
ls_Phi[is_bd_node] = -1e-6

ls_mesh.nodedata['ls_phi1'] = ls_Phi.flatten('f')

# 水平集函数值 Phi 从水平集节点投影到有限元节点
from scipy.interpolate import griddata
fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')
print("fe_Phi0:", fe_Phi.shape, "\n", fe_Phi.round(4))

fe_mesh.nodedata['fe_phi0'] = fe_Phi

ls_Phi = ts.reinitialize(phi0=ls_Phi.reshape(nelx+2, nely+2).T,
                            dx=ew, dy=eh, loop_num=50)
print("ls_NC:", ls_NN)
print("ls_Phi0:", ls_Phi.shape, "\n", ls_Phi.round(4))

ls_mesh.nodedata['ls_phi2'] = ls_Phi.flatten('f')


fname = os.path.join('./visulaization/', 'wang_ls.vtu')
ls_mesh.to_vtk(fname=fname)

fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')
print("fe_Phi1:", fe_Phi.shape, "\n", fe_Phi.round(4))

fe_mesh.nodedata['fe_Phi1'] = fe_Phi

fname = os.path.join('./visulaization/', 'wang_fe.vtu')
fe_mesh.to_vtk(fname=fname)

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
#print("tmp:", tmp)
F[tmp, 0] = -1
#print("F:", F.shape, "\n", F)

# 位移约束(supports) - short cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)
#print("fixeddofs:", fixeddofs)

E0 = 1e-3 # Young's modulus of void material
E1 = 1.0 # Young's modulus of elastic material
nu = 0.3 # 泊松比
fea_Intercal = 20; # 指定有限元单元的频率
reIter = 10 # 重初始化的时间步数


import os
totalNum = 50 # 总的迭代次数
mean_compliances = np.zeros(totalNum)
objective = np.zeros(totalNum)
# 开始循环
for iterNum in range(totalNum):
    # 有限元分析
    U = ts.fe_analysis(mesh=fe_mesh, E0=E0, E1=E1, nu=nu, ew=ew, eh=eh,
                       phi=fe_Phi, F=F, fixeddofs=fixeddofs)
    #print("U:", U.shape, "\n", U)
    ux = U[:, 0]
    uy = U[:, 1]

    mean_compliances[iterNum] = F[tmp] * U.reshape(-1, 1)[tmp]

    print(f'Iter: {iterNum}, Compliance.: {mean_compliances[iterNum]:.4f}')

    # 计算几何量
    ls_curv = ts.calc_curvature(phi=ls_Phi.reshape(nelx+2, nely+2).T, dx=ew, dy=eh)
    print("ls_curv:", ls_curv.shape, "\n", ls_curv)
    asd

    # 形状灵敏度分析
    ls_Beta = ts.sensi_analysis(fe_mesh=fe_mesh, ls_mesh=ls_mesh, E1=E1, E0=E0, \
                                u=ux, v=uy, ew=ew, eh=eh, nu=nu, lag4Vol=lagV, lag4Curv=lagCur, \
                                ele_lsgrid_id=ele_lsgrid_id, phi=ls_Phi, curvature=ls_curv)
    #print("ls_Beta:", ls_Beta.shape, "\n", ls_Beta)

    # 归一化法向速度场
    ls_Vn = ls_Beta / np.max(np.abs(ls_Beta))

    # 水平集界面更新
    ls_Phi = ts.level_set_evolve(phi0 = ls_Phi.reshape(nelx+2, nely+2).T,
                                 vn = ls_Vn.reshape(nelx+2, nely+2).T,
                                 dx = ew, dy = eh, loop_num = fea_Intercal)
    #print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi)

    # 水平集界面重置化
    if (iterNum == 0) or ((iterNum+1) % 5 == 0):
        print("iterNum:", iterNum)
        ls_Phi = ts.reinitialize(phi0=ls_Phi.reshape(nelx+2, nely+2).T, 
                                 dx=ew, dy=eh, loop_num=reIter)
        #print("re_ls_Phi:", ls_Phi.shape, "\n", ls_Phi)

    fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')
    #print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))

    fe_mesh.nodedata['u1'] = ux.flatten('F')
    fe_mesh.nodedata['u2'] = uy.flatten('F')
    fe_mesh.nodedata['fe_Phi'] = fe_Phi.flatten('F')
    fname = os.path.join('./visulaization/', f'wang{iterNum:010}.vtu')
    fe_mesh.to_vtk(fname=fname)

    #fe_node_Phi = fe_Phi.reshape(nelx+1, nely+1).T
    #ls_node_Phi = ls_Phi.reshape(nelx+2, nely+2).T

    #ax1.clear()
    #ax2.clear()

    #ax1.contourf(fe_node_x, fe_node_y, fe_node_Phi, \
    #             levels=bounds, cmap=cmap, norm=norm)
    #ax1.set_aspect('equal', adjustable='box')
    #ax1.grid(True)

    #surf = ax2.plot_surface(ls_node_x, ls_node_y, ls_node_Phi, \
    #                        cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #ax2.view_init(30, 37.5)
    #ax2.grid(True)
    #ax2.axis('equal')

    #plt.draw()
    #plt.pause(1e-3)

#plt.ioff()
#plt.show()
