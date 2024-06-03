import numpy as np

from shape_gradient import TopLsfShapeGrad

# Long Cantilever 的默认参数
domain_width = 2 # 设计区域的宽度
domain_hight = 1 # 设计区域的高度
nelx = 80 # x 方向的单元数
nely = 40 # y 方向的单元数
lagV = 50 # Lagrange multiplier for volume constraint
lagCur = 0 # Lagrange multiplier for perimeter constraint whose shape sensitivity is curvature
ts = TopLsfShapeGrad(domain_width = domain_width, domain_hight = domain_hight, \
                     nelx = nelx, nely = nely, lagV = lagV, lagCur = lagCur)
# 数据初始化
ew = domain_width / nelx
eh = domain_hight / nely
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

# 有限元网格单元中心对应的水平集网格节点的索引
ele_lsgrid_id = np.argmin(distances, axis=1)

# 初始化水平集函数
ls_Phi = ts.init_lsf(mesh = ls_mesh)
print("ls_Phi0:", ls_Phi.shape, "\n", ls_Phi)
asd
ls_mesh.nodedata['ls_Phi'] = ls_Phi
import os
fname = os.path.join('/visulaiztion', 'ls_Phi0.vtu')
ls_mesh.to_vtk(fname=fname)


# 边界条件处理
boundary_condition = (ls_x - np.min(ls_x)) * (ls_x - np.max(ls_x)) * \
                     (ls_y - np.max(ls_y)) * (ls_y - np.min(ls_y)) \
                        <= 100 * np.finfo(float).eps
ls_Phi[boundary_condition] = -1e-6

# 水平集函数值 Phi 从水平集节点投影到有限元节点
from scipy.interpolate import griddata
fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')

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


fe_node_x = fe_x.reshape(nelx+1, nely+1).T
fe_node_y = fe_y.reshape(nelx+1, nely+1).T
fe_node_Phi = fe_Phi.reshape(nelx+1, nely+1).T
ls_node_x = ls_x.reshape(nelx+2, nely+2).T
ls_node_y = ls_y.reshape(nelx+2, nely+2).T
ls_node_Phi = ls_Phi.reshape(nelx+2, nely+2).T

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors

fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122, projection='3d')
cmap = mcolors.ListedColormap(['white', 'black'])
bounds = [-np.inf, 0, np.inf]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# 初始化目标函数值列表
objective_function_values = []

plt.ion() # 开启交互模式

totalNum = 20 # 总的迭代次数
# 开始循环
for iterNum in range(totalNum):
    # 有限元分析
    U = ts.fe_analysis(mesh=fe_mesh, E0=E0, E1=E1, nu=nu, ew=ew, eh=eh,
                       phi=fe_Phi, F=F, fixeddofs=fixeddofs)
    ux = U[:, 0]
    uy = U[:, 1]
    print("U:", U.shape, "\n", U.reshape(-1, 1))

    mean_compliances = F[tmp] * U.reshape(-1, 1)[tmp]
    #print("mean_compliances:", mean_compliances)
    # 将每步的目标函数值添加到列表中
    objective_function_values.append(mean_compliances)

    # 计算几何量
    ls_curv = ts.calc_curvature(phi=ls_Phi.reshape(nelx+2, nely+2).T, dx=ew, dy=eh)

    # 形状灵敏度分析
    ls_Beta = ts.sensi_analysis(fe_mesh=fe_mesh, ls_mesh=ls_mesh, E1=E1, E0=E0, \
                                u=ux, v=uy, ew=ew, eh=eh, nu=nu, lag4Vol=lagV, lag4Curv=lagCur, \
                                ele_lsgrid_id=ele_lsgrid_id, phi=ls_Phi, curvature=ls_curv)

    # 归一化法向速度场
    ls_Vn = ls_Beta / np.max(np.abs(ls_Beta))

    # 水平集界面更新
    ls_Phi = ts.level_set_evolve(phi0 = ls_Phi.reshape(nelx+2, nely+2).T,
                                 vn = ls_Vn.reshape(nelx+2, nely+2).T,
                                 dx = ew, dy = eh, loop_num = fea_Intercal)

    # 水平集界面重置化
    if (iterNum == 0) or ((iterNum+1) % 5 == 0):
        ls_Phi = ts.reinitialize(phi0=ls_Phi.reshape(nelx+2, nely+2).T, 
                                 dx=ew, dy=eh, loop_num=reIter)

    fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')


    fe_node_Phi = fe_Phi.reshape(nelx+1, nely+1).T
    ls_node_Phi = ls_Phi.reshape(nelx+2, nely+2).T

    ax1.clear()
    ax2.clear()

    ax1.contourf(fe_node_x, fe_node_y, fe_node_Phi, \
                 levels=bounds, cmap=cmap, norm=norm)
    ax1.set_aspect('equal', adjustable='box')
    ax1.grid(True)

    surf = ax2.plot_surface(ls_node_x, ls_node_y, ls_node_Phi, \
                            cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax2.view_init(30, 37.5)
    ax2.grid(True)
    ax2.axis('equal')

    plt.draw()
    plt.pause(1e-3)

plt.ioff()
plt.show()

# 绘制目标函数的收敛曲线
plt.figure(figsize=(10, 6))
plt.plot(objective_function_values, label='Objective Function')
plt.xlabel('Iterations')
plt.ylabel('Objective function')
plt.title('Convergence of Objective Function')
plt.legend()
plt.grid(True)
plt.show()
