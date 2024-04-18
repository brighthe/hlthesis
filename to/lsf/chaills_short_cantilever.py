import numpy as np
import os

from numpy._typing import _8Bit

from chaills import TopLsf

# Short Cantilever 的默认参数
nelx = 6
nely = 4
volReq = 0.4
stepLength = 2;
topWeight = 2;
numReinit = 3
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq, stepLength=stepLength, topWeight=topWeight, numReinit=3)

# 初始化优化参数
nelx, nely, volReq = ts._nelx, ts._nely, ts._volReq
#mesh = ts._mesh

struc_mesh = ts.generate_mesh(domain=[0, nelx, 0, nely], nelx=nelx, nely=nely)
lsf_mesh = ts.generate_mesh(domain=[0, nelx+2, 0, nely+2], nelx=nelx+2, nely=nely+2)

# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))

# 初始化水平集函数
lsf = ts.reinit(struc = struc)
print("node:", lsf_mesh.entity('node').shape)
lsf_mesh.celldata['lsf'] = lsf.flatten('F')
fname = os.path.join('./visulaization/', 'chaills_lsf.vtu')
lsf_mesh.to_vtk(fname=fname)

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
        0,1 - 4,5
        2,3 - 6,7
        拓扑 : cell2dof : 0,1,6,7,8,9,2,3
        FEALPy : cell2dof : 0,1,2,3,6,7,8,9
    """
    k0 = np.array([1/2 - nu/6,   1/8 + nu/8,   -1/4 - nu/12, -1/8 + 3 * nu/8,
                -1/4 + nu/12,  -1/8 - nu/8,    nu/6,         1/8 - 3 * nu/8])
    k_0 = 1/2 - nu/6
    k_1 = 1/8 + nu/8
    k_2 = nu/6
    k_3 = 1/8 - 3*nu/8
    k_4 = -1/4 - nu/12
    k_5 = -1/8 + 3*nu/8
    k_6 = -1/4 + nu/12
    k_7 = -1/8 - nu/8
    k = np.array([k_0, k_1, k_2, k_3, k_4, k_5, k_6, k_7])


    print("k0:", k0)
    print("k:", k)
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
    print("KE:", KE.round(4))

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
# 定义荷载 - short cantilever
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
#tmp = vgdof - 1
tmp = vgdof - 2*(nely+1) - 1
print("tmp:", tmp)
F[tmp, 0] = -1

# 位移约束(supports) - short cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)
print("fixeddofs:", fixeddofs)

# Load Bearing 单元的索引 - short cantilever
loadBearingIndices = (-1, -1)
#loadBearingIndices = (nely, -1)
print("loadBearingIndices:", loadBearingIndices)

# 设置初始的 augmented Lagrangian parameters
la = -0.01
La = 1000
alpha = 0.9
# 优化循环的最大迭代次数
num = 200
# 初始化 compliance objective value
objective = np.zeros(num)
for iterNum in range(num):

    # 有限元计算全局位移和局部单元位移
    U, Ue = ts.FE(mesh=mesh, struc=struc, KE=KE, F=F, fixeddofs=fixeddofs)

    # 添加来自每个荷载的灵敏度之前，将形状和拓扑灵敏度设置为 0
    shapeSens[:] = 0
    topSens[:] = 0

    for i in range(nLoads):
        # 计算每个单元的柔度的形状灵敏度
        stiff = 0.0001 # 0.0001: Stiffness of the void phase compared to the solid phase
                    # 可以为零，但较小的非零值可以提高算法的稳健性
        temp1 = -np.maximum(struc, stiff)
        temp2 = np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        shapeSens[:] = shapeSens[:] + np.einsum('ij, ij -> ij', temp1, temp2)

        # 计算每个单元的柔度的拓扑灵敏度
        coef = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
        temp3 = (4 * mu) * \
            np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        temp4 = (lambda_ - mu) * \
            np.einsum('ij, jk, ki -> i', Ue[:, :, i], KTr, Ue[:, :, i].T).reshape(nelx, nely).T
        topSens[:] = topSens[:] + np.einsum('ij, ij -> ij', coef*struc, (temp3+temp4))

    # 存储当前迭代的 compliance objective
    objective[iterNum] = -np.sum(shapeSens)

    # 计算当前的 volume fraction
    volCurr = np.sum(struc) / (nelx*nely)

    # 打印当前迭代的结果
    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}')

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
        pass
    else:
        # TODO 与理论不一致
        la = la - 1/La * (volCurr - volReq)
        La = alpha * La

    # Update the sensitivities with augmented Lagrangian terms
    # TODO 与理论不一致
    shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
    topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )

    # 灵敏度光滑化
    shapeSens = ts.smooth_sens(sens=shapeSens)
    topSens = ts.smooth_sens(sens=topSens)

    # 执行设计更新
    struc, lsf = ts.updateStep(lsf=lsf, shapeSens=shapeSens, topSens=topSens,
                               stepLength=stepLength, topWeight=topWeight,
                               loadBearingIndices=loadBearingIndices)

    #print("struc:", struc.shape, "\n", "struc")
    #print("lsf:", lsf.shape, "\n", "lsf")

    # 在特定的迭代步重置化水平集函数
    if (iterNum+1) % numReinit == 0:
        lsf = ts.reinit(struc)

    mesh.celldata['struc'] = struc.flatten('F')
    fname = os.path.join('./visulaization/', f'chaills{iterNum:010}.vtu')
    mesh.to_vtk(fname=fname)
