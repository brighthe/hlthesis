import numpy as np

from top_simp import TopSimp

# Cantilever beam with a fixed hole
nelx = 45
nely = 30
volfrac = 0.5
penal = 3.0
rmin = 1.5
ts = TopSimp(nelx=nelx, nely=nely, volfrac=volfrac, penal=penal, rmin=rmin)

# 初始化优化参数
nelx, nely, volfrac, penal, rmin = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, markersize=6, fontsize=6, fontcolor='r')
#mesh.find_cell(axes, showindex=True, markersize=6, fontsize=6, fontcolor='g')
#plt.show()

# 根据体积分数 volfrac 初始化设计变量场
x = np.full((nely, nelx), volfrac)

# 利用 passive 单元给结构打孔洞
passive = np.zeros((nely, nelx))
for ely in range(nely):
    for elx in range(nelx):
        if np.sqrt((ely+1 - nely/2.)**2 + (elx+1 - nelx/3.)**2) < nely / 3.:
            passive[ely, elx] = 1
            x[ely, elx] = 0.001

# 初始化每个单元的目标函数的灵敏度为 0
dc = np.zeros((nely, nelx))

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

E0 = 1.0
nu = 0.3
KE = stiff_matrix(nu=nu, E0=E0)

from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 定义荷载 - Cantilever beam with a fixed hole
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
F[vgdof-1, 0] = -1

# 位移约束(supports) - Cantilever beam with a fixed hole
fixeddofs = np.arange(0, 2*(nely+1), 1)

import time
time_stats = {
    'FE_computation': [],
    'c_dc_computation': [],
    'filter_dc': [],
    'update_x': [],
    'plot': [],
    'iteration_time': []
}

# 绘制材料密度图 - (白色 - 黑色): (0 - 1)
import matplotlib.pyplot as plt
plt.ion()
fig, ax = plt.subplots()
image = ax.imshow(-x, cmap='gray', vmin=-1, vmax=0)
ax.axis('off')

# 迭代次数
loop = 0
# 迭代中设计变量的最大变化
change = 1.0
# 优化循环，直到 change < 0.01 迭代终止
while change > 0.01:
    iter_start_time = time.time()

    loop += 1
    xold = np.copy(x)

    # 有限元计算全局位移和局部单元位移
    fe_start = time.time()

    U, Ue = ts.FE(mesh=mesh, x=x, penal=penal, KE=KE, F=F, fixeddofs=fixeddofs)
    print("U:", U.shape, "\n", U.round(4))
    print("Ue:", Ue.shape, "\n", Ue[620].round(4))

    fe_end = time.time()
    fe_time = fe_end - fe_start
    time_stats['FE_computation'].append(fe_time)

    # 初始化目标函数为 0
    c = 0
    # 初始化每个单元的目标函数的灵敏度为 0
    dc = np.zeros((nely, nelx))
    c_dc_time = 0
    for i in range(nLoads):
        # 计算目标函数
        c_dc_statr = time.time()

        temp1 = x ** penal
        temp2 = np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        c = c + np.einsum('ij, ij -> ', temp1, temp2)

        # 计算每个单元的目标函数的灵敏度
        temp3 = -penal * x ** (penal-1)
        dc[:] = dc[:] + np.einsum('ij, ij -> ij', temp3, temp2)

        c_dc_end = time.time()
        c_dc_time = c_dc_end - c_dc_statr
        time_stats['c_dc_computation'].append(c_dc_time)


    # 灵敏度过滤
    filter_dc_start = time.time()

    dc = ts.check(rmin=rmin, x=x, dc=dc)
    #print("dc:", dc.shape, "\n", dc.round(4))

    filter_dc_end = time.time()
    filter_dc_time = filter_dc_end - filter_dc_start
    time_stats['filter_dc'].append(filter_dc_time)

    # 使用 Optimality Criteria Method 更新设计
    update_x_start = time.time()

    x = ts.OC(volfrac=volfrac, x=x, dc=dc, passive=passive)
    print("x:", x.shape, "\n", x.round(4))

    update_x_end = time.time()
    update_x_time = update_x_end - update_x_start
    time_stats['update_x'].append(update_x_time)

    # 计算当前的 volume fraction
    volfrac = np.sum(x) / (nelx*nely)

    # 打印当前迭代的结果
    change = np.max(np.abs(x - xold))
    print(f' Iter.: {loop:4d} Objective.: {c:10.4f} Volfrac.: {volfrac:6.3f} change.: {change:6.3f}')

    # 更新图像
    plot_start = time.time()

    image.set_data(-x)
    plt.draw()
    plt.pause(1e-5)

    plot_end = time.time()
    plot_time = plot_end - plot_start
    time_stats['plot'].append(plot_time)

    iter_end_time = time.time()
    total_iter_time = iter_end_time - iter_start_time
    time_stats['iteration_time'].append(total_iter_time)

    print(f'Total: {total_iter_time:.4f}, FE: {fe_time:.4f}, c_dc: {c_dc_time:.4f},'
          f'Filter_dc: {filter_dc_time:.4f}, Update_x: {update_x_time:.4f}, plot: {plot_time:.4f}')

plt.ioff()
plt.show()
