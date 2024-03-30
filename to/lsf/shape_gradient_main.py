import numpy as np
from numpy.lib.ufunclike import fix

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

ls_Phi = ts.init_lsf(mesh = ls_mesh)
#print("ls_Phi0:", ls_Phi.shape, "\n", ls_Phi.round(4))

#ts.plot(x=ls_x, y=ls_y, z=ls_Phi, label='ls_Phi')

# 边界条件处理
boundary_condition = (ls_x - np.min(ls_x)) * (ls_x - np.max(ls_x)) * \
                     (ls_y - np.max(ls_y)) * (ls_y - np.min(ls_y)) \
                        <= 100 * np.finfo(float).eps
#print("boundary_condition:", boundary_condition)
ls_Phi[boundary_condition] = -1e-6
#print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi.round(4))

# 水平集函数值 Phi 从水平集节点投影到有限元节点
from scipy.interpolate import griddata
fe_Phi = griddata((ls_x, ls_y), ls_Phi, (fe_x, fe_y), method='cubic')
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
F[vgdof - (nely+1), 0] = 1
#print("F:", F.shape, "\n", F)

# 位移约束(supports) - short cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)
print("fixeddofs:", fixeddofs)

totalNum = 1
# 开始循环
for iterNum in range(totalNum):
    # 有限元分析
    print(f'Finite Element Analysis No.: {iterNum}')

    E0 = 1e-3
    E1 = 1.0
    nu = 0.3

    fe_Phi = np.array([1, 2, 3, 4, 5, 0, 0, -0.5, -1, -2, -3, -4])
    print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))
    print(fe_cell)
    U = ts.fe_analysis(mesh=fe_mesh, E0=E0, E1=E1, nu=nu, ew=ew, eh=eh, Phi=fe_Phi,
                       F=F, fixeddofs=fixeddofs)
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

distances = (ls_x[None, :] - fe_center_x[:, None])**2 + \
            (ls_y[None, :] - fe_center_y[:, None])**2
print("distances:", distances.shape, "\n", distances)

# 有限元网格单元中心对应的水平集网格节点的索引
ele_lsgrid_id = np.argmin(distances, axis=1)
print("ele_lsgrid_id:", ele_lsgrid_id.shape, "\n", ele_lsgrid_id)

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

























# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))
print("struc:", struc.shape, "\n", struc)

# 水平集函数的初始化
#r = nely * 0.1 # 初始孔洞的半径
#lsf = ts.lsf_init(mesh = mesh, r=r)
lsf = ts.reinit(struc = struc)
print("lsf0:", lsf.shape, "\n", lsf.round(4))
strucFULL = (lsf < 0).astype(int)
struc = strucFULL[1:-1, 1:-1]
print("struc:", struc.shape, "\n", struc)
asd

# 初始化灵敏度
shapeSens = np.zeros((nely, nelx))

def stiff_matrix(nu, E0):
    """
    单元刚度矩阵.
    
    Parameters:
    - nu (flaot): Poisson 比.
    - E0 (int): 杨氏模量.

    Returns:
    - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.

    Note:
        考虑四个单元的四边形网格中的 0 号单元：
        0,1 - 6,7
        2,3 - 8,9
        拓扑 : cell2dof : 0,1,6,7,8,9,2,3
        FEALPy : cell2dof : 0,1,2,3,6,7,8,9
    """
    k = np.array([1/2 - nu/6,   1/8 + nu/8,   -1/4 - nu/12, -1/8 + 3 * nu/8,
                -1/4 + nu/12,  -1/8 - nu/8,    nu/6,         1/8 - 3 * nu/8])
    KE = E0 / (1 - nu**2) * np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
    ])

    return KE

def lame(nu, E0):
    """
    Lame 常数.

    Parameters:
    - nu (flaot): Poisson 比.
    - E0 (int): 杨氏模量.
    """
    lambda_ = E0 * nu / ((1 + nu) * (1 - nu))
    mu = E0 / (2 * (1 + nu))

    return lambda_, mu

E0 = 1.0
nu = 0.3
KE = stiff_matrix(nu=nu, E0=E0)
lambda_, mu = lame(nu=nu, E0=E0)

from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 定义荷载 - medium cantilever
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
nodal_loads_index = 2*( (nely+1)*nelx + int(np.ceil(nely/2)) ) + 1
F[nodal_loads_index, 0] = -1

# 位移约束(supports) - medium cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)

# Load Bearing 单元的索引 - short cantilever
loadBearingIndices = (-1, -1)

import time
time_stats = {
    'FE_computation': [],
    'shapeSens_computation': [],
    'plot': [],
    'augmented_lag': [],
    'design_update': [],
    'reinit_lsf': [],
    'iteration_time': []
}

# 绘制结果图
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# 设置初始的 augmented Lagrangian parameters
la = -0.01
La = 1000
alpha = 0.9
# 优化循环的最大迭代次数
num = 200
# 初始化 compliance objective value
objective = np.zeros(num)
for iterNum in range(num):
    iter_start_time = time.time()

    # 有限元计算全局位移和局部单元位移
    fe_start = time.time()

    U, Ue = ts.FE(mesh=mesh, struc=struc, KE=KE, F=F, fixeddofs=fixeddofs)

    fe_end = time.time()
    fe_time = fe_end - fe_start
    time_stats['FE_computation'].append(fe_time)

    # 添加来自每个荷载的灵敏度之前，将形状和拓扑灵敏度设置为 0
    shapeSens[:] = 0

    shapeSens_time = 0
    #topSens_time = 0
    for i in range(nLoads):
        # 计算每个单元的柔度的形状灵敏度
        shapeSens_start = time.time()

        stiff = 0.0001 # 0.0001: Stiffness of the void phase compared to the solid phase
                    # 可以为零，但较小的非零值可以提高算法的稳健性
        temp1 = -np.maximum(struc, stiff)
        temp2 = np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        shapeSens[:] = shapeSens[:] + np.einsum('ij, ij -> ij', temp1, temp2)

        shapeSens_end = time.time()
        shapeSens_time = shapeSens_end - shapeSens_start
        time_stats['shapeSens_computation'].append(shapeSens_time)

    # 存储当前迭代的 compliance objective
    objective[iterNum] = -np.sum(shapeSens)

    # 计算当前的 volume fraction
    volCurr = np.sum(struc) / (nelx*nely)

    # 打印当前迭代的结果
    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}')

    # 更新图像
    plot_start = time.time()

    # 白色区域：void 部分, 黑色区域：solid 部分
    cmap = mcolors.ListedColormap(['white', 'black'])
    bounds = [-np.inf, 0, np.inf]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(1)
    contourf = plt.contourf(lsf, levels=bounds, cmap=cmap, norm=norm)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-5)

    plot_end = time.time()
    plot_time = plot_end - plot_start
    time_stats['plot'].append(plot_time)


    # 五次迭代后执行收敛性检查
    start_num = 5 # Number of iterations at the start of the optimization 
                # for which the convergence criteria are not checked
    vol_tor = 0.005 # Tolerance for satisfaction of the volume constraint
    rel_tor = 0.01 # Relative tolerance on the objective values for termination of the algorithm
    if iterNum > start_num and (abs(volCurr-volReq) < vol_tor) and \
        np.all( np.abs(objective[iterNum] - objective[iterNum-start_num:iterNum]) \
               < rel_tor * np.abs(objective[iterNum]) ):
        break

    # 设置  augmented Lagrangian parameters
    lag_start = time.time()

    if iterNum == 0:
        pass
    else:
        # TODO 与理论不一致
        la = la - 1/La * (volCurr - volReq)
        La = alpha * La

    # Update the sensitivities with augmented Lagrangian terms
    # TODO 与理论不一致
    shapeSens = shapeSens - la + 1/La * (volCurr - volReq)

    lag_end = time.time()
    lag_time = lag_end - lag_start
    time_stats['augmented_lag'].append(lag_time)

    # 执行设计更新
    update_start = time.time()

    struc, lsf = ts.updateStep(lsf=lsf, shapeSens=shapeSens, stepLength=stepLength, loadBearingIndices=loadBearingIndices)

    update_end = time.time()
    update_time = update_end - update_start
    time_stats['design_update'].append(update_time)

    # 在特定的迭代步重置化水平集函数
    reinit_time = 0
    if (iterNum+1) % numReinit == 0:
        reinit_start = time.time()

        lsf = ts.reinit(struc)
        print("reinit_lsf:\n", lsf.shape, "\n", lsf.round(4))

        reinit_end = time.time()
        reinit_time = reinit_end - reinit_start
        time_stats['reinit_lsf'].append(reinit_time)

    iter_end_time = time.time()
    total_iter_time = iter_end_time - iter_start_time
    time_stats['iteration_time'].append(total_iter_time)

    reinit_time_str = f", Reinit_lsf: {reinit_time:.4f}" if reinit_time > 0 else ""
    print(f'Total: {total_iter_time:.4f}, FE: {fe_time:.4f}, '
          f'ShapeSens: {shapeSens_time:.4f}, '
          f'plot: {plot_time:.4f}, AugmentedLag: {lag_time:.4f}, '
          f'Update: {update_time:.4f}{reinit_time_str}, '
          f'iter_time: {total_iter_time - plot_time:.4f}')

plt.ioff()
plt.show()
