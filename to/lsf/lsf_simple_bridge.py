import numpy as np

from top_lsf import TopLsf

# Simple Bridge 的默认参数
nelx = 60
nely = 20
volReq = 0.3
stepLength = 3;
topWeight = 2;
numReinit = 2
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq, stepLength=stepLength, topWeight=topWeight, numReinit=3)

# 初始化优化参数
nelx, nely, volReq = ts._nelx, ts._nely, ts._volReq
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

# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))

# 初始化水平集函数
lsf = ts.reinit(struc = struc)

# 初始化灵敏度
shapeSens = np.zeros((nely, nelx))
topSens = np.zeros((nely, nelx))

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

def trace_matrix(nu, E0):
    """
    应力张量的迹矩阵.

    Parameters:
    - nu (flaot): Poisson 比.
    - E0 (int): 杨氏模量.

    Returns:
    - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元迹矩阵.
    """
    k = np.array([1/3, 1/4, -1/3, 1/4, -1/6, -1/4, 1/6, -1/4])
    KTr = E0 / (1 - nu) * np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
    ])

    return KTr

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
KTr = trace_matrix(nu=nu, E0=E0)
lambda_, mu = lame(nu=nu, E0=E0)

from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 定义荷载 - simple bridge
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
F[2*(round(nelx/2)+1)*(nely+1)-1, 0] = -1

# 位移约束(supports) - simple bridge
fixeddofs = np.concatenate( [np.arange( 2*(nely+1)-2, 2*(nely+1) ),
                             np.arange( 2*(nelx+1)*(nely+1)-2, 2*(nelx+1)*(nely+1) )] )

# Load Bearing 单元的索引 - simple bridge
loadBearingIndices = ( -1, [0, nelx//2-1, nelx//2, -1] )

import time
time_stats = {
    'FE_computation': [],
    'shapeSens_computation': [],
    'topSens_computation': [],
    'plot': [],
    'augmented_lag': [],
    'sensitivity_smoothing': [],
    'design_update': [],
    'reinit_lsf': [],
    'iteration_time': []
}

# 绘制结果图 - 白色区域：void 部分, 黑色区域：solid 部分
import matplotlib.pyplot as plt
plt.ion()
fig, ax = plt.subplots()
image = ax.imshow(-struc, cmap='gray', vmin=-1, vmax=0)
ax.axis('off')

# 设置  初始的 augmented Lagrangian parameters
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
    topSens[:] = 0

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

        # 计算每个单元的柔度的拓扑灵敏度
        topSens_start = time.time()

        coef = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
        temp3 = (4 * mu) * \
            np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        temp4 = (lambda_ - mu) * \
            np.einsum('ij, jk, ki -> i', Ue[:, :, i], KTr, Ue[:, :, i].T).reshape(nelx, nely).T
        topSens[:] = topSens[:] + np.einsum('ij, ij -> ij', coef*struc, (temp3+temp4))

        topSens_end = time.time()
        topSens_time = topSens_end - topSens_start
        time_stats['topSens_computation'].append(topSens_time)

    # 存储当前迭代的 compliance objective
    objective[iterNum] = -np.sum(shapeSens)

    # 计算当前的 volume fraction
    volCurr = np.sum(struc) / (nelx*nely)

    # 打印当前迭代的结果
    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}')

    # 更新图像
    plot_start = time.time()

    image.set_data(-struc)
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
    topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )

    lag_end = time.time()
    lag_time = lag_end - lag_start
    time_stats['augmented_lag'].append(lag_time)

    # 灵敏度光滑化
    smoothing_start = time.time()

    shapeSens = ts.smooth_sens(sens=shapeSens)
    topSens = ts.smooth_sens(sens=topSens)

    smoothing_end = time.time()
    smoothing_time = smoothing_end - smoothing_start
    time_stats['sensitivity_smoothing'].append(smoothing_time)

    # 执行设计更新
    update_start = time.time()

    struc, lsf = ts.updateStep(lsf=lsf, shapeSens=shapeSens, topSens=topSens,
                               stepLength=stepLength, topWeight=topWeight,
                               loadBearingIndices=loadBearingIndices)
    print("lsf:\n", lsf.shape, "\n", lsf.round(4))
    print("struc:", struc.shape , "\n", struc.round(4))

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
    print(f'Total: {total_iter_time:.4f}, FE: {fe_time:.4f}, ShapeSens: {shapeSens_time:.4f},'
          f'topSens: {topSens_time:.4f}, plot: {plot_time:.4f}, AugmentedLag: {lag_time:.4f}, '
          f'Smoothing: {smoothing_time:.4f}, Update: {update_time:.4f}{reinit_time_str}, '
          f'iter_time: {total_iter_time - plot_time:.4f}')

plt.ioff()
plt.show()
