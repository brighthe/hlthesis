import numpy as np

from lsf_top_short_cantilever import TopLsf

# Short Cantilever
nelx = 32
nely = 20
volReq = 0.4
stepLength = 2;
topWeight = 2;
numReinit = 3
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq, stepLength=stepLength, topWeight=topWeight, numReinit=3)

# 初始化优化参数
nelx, nely, volReq = ts._nelx, ts._nely, ts._volReq
mesh = ts._mesh_top2

#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
#mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
#plt.show()

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

# 定义开始优化时的初始结构为 entirely solid
struc = np.ones((nely, nelx))

# 初始化水平集函数
lsf = ts.reinit(struc = struc)

# 初始化灵敏度
shapeSens = np.zeros((nely, nelx))
topSens = np.zeros((nely, nelx))

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
E0 = 1.0
nu = 0.3
integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
KE = integrator.stiff_matrix()
KTr = integrator.trace_matrix()
lambda_, mu = integrator.lame()

# 优化循环的最大迭代次数
num = 200
# 初始化 compliance objective value
objective = np.zeros(num)
for iterNum in range(num):
    # 计算全局位移和局部单元位移
    U, Ue = ts.FE(mesh=mesh, struc=struc)

    # 计算每个单元的柔度的形状灵敏度
    stiff = 0.0001 # 0.0001: Stiffness of the void phase compared to the solid phase
                # 可以为零，但较小的非零值可以提高算法的稳健性
    temp1 = -np.maximum(struc, stiff)
    temp2 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nelx, nely).T
    shapeSens[:] = np.einsum('ij, ij -> ij', temp1, temp2)

    # 计算每个单元的柔度的拓扑灵敏度
    coef = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
    temp3 = (4 * mu) * np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nelx, nely).T
    temp4 = (lambda_ - mu) * np.einsum('ij, jk, ki -> i', Ue, KTr, Ue.T).reshape(nelx, nely).T
    topSens[:] = np.einsum('ij, ij -> ij', coef*struc, (temp3+temp4))

    # 存储当前迭代的 compliance objective
    objective[iterNum] = -np.sum(shapeSens)

    # 计算当前的 volume fraction
    volCurr = np.sum(struc) / (nelx*nely)

    # 打印当前迭代的结果
    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}')

    # 绘制结果图
    import matplotlib.pyplot as plt
    plt.imshow(-struc, cmap='gray', vmin=-1, vmax=0)
    plt.axis('off')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-5)

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
    if iterNum == 0:
        la = -0.01
        La = 1000
        alpha = 0.9
    else:
        # TODO 与理论不一致
        la = la - 1/La * (volCurr - volReq)
        La = alpha * La

    # Update the sensitivities with augmented Lagrangian terms
    shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
    topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )

    # Smooth the sensitivities
    shapeSens = ts.smooth_sens(sens=shapeSens)
    topSens = ts.smooth_sens(sens=topSens)

    # 执行设计更新
    struc, lsf = ts.updateStep(lsf=lsf, shapeSens=shapeSens, topSens=topSens,
                               stepLength=stepLength, topWeight=topWeight)

    # Reinitialize the level set function at specified iterations
    if (iterNum+1) % numReinit == 0:
        lsf = ts.reinit(struc)

plt.ioff()
plt.show()
