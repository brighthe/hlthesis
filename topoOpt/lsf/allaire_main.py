from operator import xor
import numpy as np
import os

# Cantilever 的默认参数
nelx = 80
nely = 40
xlength = 2  # 工作域长度
dx = xlength / nelx  # x 方向的步长
dy = dx  # y 方向的步长，与 x 方向相同
yheight = dy * nely  # 工作域高度

hx = 5  # x 方向孔洞数
hy = 3  # y 方向孔洞数
r = 0.5  # 孔洞大小（介于0和1之间）

lagV = 50 # Lagrange multiplier for volume constraint
lagP = 0 # Lagrange multiplier for perimeter constraint 
e2 = 4*dx**2 # coefficient for the regularization in front of the Laplacian

nodex = np.linspace(0, xlength, nelx + 1)  # x space for nodes
nodey = np.linspace(0, yheight, nely + 1)  # y space for nodes
print("nodey:", nodey.shape, "\n", nodey)
FEx = np.linspace(0, xlength, nelx)  # x space for finite elements
FEy = np.linspace(0, yheight, nely)  # y space for finite elements

RIiterinit = 50 # number of time steps in the reinitialization of the initial design
RIiter = 5 # number of time steps in the reinitialization of further designs
RIfreq = 5 # 重置化的频率，即两次重置化之间的 transport level set equation 的时间步

HJiter0 = 20 # original number of transport time steps 
entire = 20 # total number of optimization iterations
allow = 0.01 # fraction that the objective function can increase

from allaire import TopLSM
ts = TopLSM(nelx=nelx, nely=nely, xlength=xlength, yheight=yheight)

from fealpy.mesh import QuadrangleMesh
domain = [0, xlength, 0, yheight]
mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)

# 初始化
#phi0 = ts.init_lsf(mesh)
phi0 = ts.mesh0(hx, hy, r)
print("phi0:", phi0.shape, "\n", phi0.reshape(nelx+1, nely+1).T)

mesh.nodedata['phi0'] = phi0.flatten('F')
#mesh.nodedata['phi1'] = phi1.flatten('F')


# 重置化水平集函数
#phi00 = ts.mesh00(phi=phi0, RIiter=50)
#print("phi00:", phi00.shape, "\n", phi00)
phi0 = phi0.reshape((nelx+1, nely+1)).T
phi00 = ts.reinitialize(phi0=phi0, dx=dx, dy=dy, loop_num=50)
print("phi00:", phi00.shape, "\n", phi00)

mesh.nodedata['phi00'] = phi00.flatten('F')
#mesh.nodedata['phi11'] = phi11.flatten('F')

phi = phi00

# 基于 phi 定义单元密度
eps = 0.001 # "hole" (or ersatz) material density
fe_theta = ts.fe_density(phi.reshape(nelx+1, nely+1).T, eps)
print("fe_theta:", fe_theta.shape, "\n", fe_theta)

mesh.celldata['fe_theta'] = fe_theta.flatten('F')
fname = os.path.join('./visulaization/', "allaire.vtu")
mesh.to_vtk(fname = fname)

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

E1 = 1.0 # Young's modulus of elastic material
nu = 0.3 # 泊松比
KE = stiff_matrix(nu=nu, E0=E1)
from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
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

e3 = 1

from scipy.sparse import lil_matrix
# 有限差分矩阵
K1 = lil_matrix(((nelx + 1) * (nely + 1), (nelx + 1) * (nely + 1))) # 速度正则化的矩阵

num = 50
objective = np.zeros(num)
fea_eueu = np.zeros((nely, nelx))
for iterNum in range(num):
    U, Ue = ts.fe_analysis(mesh=mesh, fe_theta=fe_theta, KE=KE, F=F, fixeddofs=fixeddofs)
    #print("U:", U.shape, "\n", U)
    #print("Ue:", Ue.shape, "\n", Ue)

    # 计算每个单元的应变能密度
    for i in range(nLoads):
        temp = np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T).reshape(nelx, nely).T
        fea_eueu[:] = fea_eueu[:] + np.einsum('ij, ij -> ij', fe_theta, temp)

        # 创建一个比fea_eueu大一圈的数组，外围填充零
        padded = np.pad(fea_eueu, ((1, 1), (1, 1)), mode='edge')
        lvlAeueu = np.zeros((nely+1, nelx+1))
        # 累加四个角的应变能密度
        lvlAeueu += padded[:-1, :-1]   # 从左上角开始
        lvlAeueu += padded[:-1, 1:]  # 从右上角开始
        lvlAeueu += padded[1:, :-1]  # 从左下角开始
        lvlAeueu += padded[1:, 1:] # 从右下角开始
        # 计算平均值
        lvlAeueu /= 4
    #print("fea_eueu:", fea_eueu.shape, "\n", fea_eueu)
    #print("lvlAeueu:", lvlAeueu.shape, "\n", lvlAeueu)

    # 定义速度场
    v = lvlAeueu / (dx*dy) - lagV
    print("v:", v.shape, "\n", v)
    #v = ts.regularize(phi=phi, V=v, e2=e2)
    #print("v1:", v.shape, "\n", v)

    # 速度场的正则化

    # 计算目标函数
    totcomp = np.sum(fea_eueu)
    totvol = ts.volume(fe_theta)
    objective[iterNum] = lagV * totvol + totcomp + lagP * ts.perimeter(phi.reshape(nelx+1, nely+1).T)
    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Totvol.: {totvol:.3f}')

    dt = 0.5 * e3 * min(dx, dy) / np.max(np.abs(v)) ;
    phi = ts.solvelevelset(phi.reshape(nelx+1, nely+1).T, v, dt, HJiter0, lagP) ;
    #print("phi:", phi.shape, "\n", phi)

    if (iterNum+1) % 5 == 0:
        phi = ts.mesh00(phi=phi, RIiter=50)

    fe_theta = ts.fe_density(phi, eps)
    #print("fe_theta:", fe_theta.shape, "\n", fe_theta)

    mesh.celldata['fe_theta'] = fe_theta.flatten('F')
    mesh.nodedata['phi'] = phi.flatten('F')
    fname = os.path.join('./visulaization/', f'allarie_{iterNum:010}.vtu')
    mesh.to_vtk(fname=fname)








asd
'''
有限元矩阵
'''
KE = ts.lk() # 刚度矩阵
print("KE:", KE.shape, "\n", KE)
K = lil_matrix((2 * (nelx + 1) * (nely + 1), 2 * (nelx + 1) * (nely + 1))) # 全局刚度矩阵
F = lil_matrix((2 * (nely + 1) * (nelx + 1), 2)) # 力矩阵
U = lil_matrix((2 * (nely + 1) * (nelx + 1), 2)) # 位移向量矩阵

# 有限差分矩阵
K1 = lil_matrix(((nelx + 1) * (nely + 1), (nelx + 1) * (nely + 1))) # 速度正则化的矩阵
print("K1:", K1.shape, "\n", K1)

# 弹性问题的设置
forceMatrix = np.array([[1, 0.5, 1, -1]])
forceCount = forceMatrix.shape[0]  # 施加力的数量
for force in range(forceCount):
    # 计算力的作用节点索引
    node_index = int((forceMatrix[force, 0] * (nelx + 1) + forceMatrix[force, 1] * (nely + 1)) * 2)
    if forceMatrix[force, 2] == 0:  # 水平力
        F[node_index, force] = forceMatrix[force, 3]
    else:  # 垂直力
        F[node_index + 1, force] = forceMatrix[force, 3]
print("F:", F.shape, "\n", F)

# 固定边界条件
fixeddofs = np.arange(2 * (nely + 1))  # 固定x和y自由度的节点
alldofs = np.arange(2 * (nely + 1) * (nelx + 1))
freedofs = np.setdiff1d(alldofs, fixeddofs)  # 计算自由度

from fealpy.mesh import QuadrangleMesh
domain = [0, xlength, 0, yheight]
mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)
node = mesh.entity('node')
print("node:", node.shape, "\n", node)


# 重置化水平集函数
phi00 = ts.mesh00(phi0, rIiterinit)
print("phi00:", phi00.shape, "\n", phi00)

asd

# 检查拓扑梯度变量 'ntop' 是否已定义
try:
    ntop
except NameError:
    ntop = 0



def g(u1, u2, v1, v2):
    """
    数值通量函数
    这是 Hamilton-Jacobi 方程数值方案中的通量函数。
    """
    return np.sqrt(np.maximum(u1, 0)**2 + np.minimum(u2, 0)**2 + 
                   np.maximum(v1, 0)**2 + np.minimum(v2, 0)**2)


def mesh00(phi, RIiter):
    """
    重初始化水平集函数
    在初始化之后以及在求解运输方程的过程中周期性地重初始化水平集函数。
    这意味着我们求解方程 d(phi)/dt + sign(phi0)(|grad(phi)|-1) = 0，
    其中 phi(t=0, x) = phi_0(x)。
    我们这样做是为了确保水平集函数在零水平集附近（即我们的结构边界）不会变得过平或过陡。
    """
    cfl = 0.5  # CFL 条件
    dt0 = min(dx, dy) * cfl  # 根据 CFL 条件定义时间步长

    for _ in range(RIiter):
        phin = np.roll(phi, -1, axis=0)
        phis = np.roll(phi, 1, axis=0)
        phie = np.roll(phi, -1, axis=1)
        phiw = np.roll(phi, 1, axis=1)

        phinn = np.roll(phin, -1, axis=0)
        phiss = np.roll(phis, 1, axis=0)
        phiee = np.roll(phie, -1, axis=1)
        phiww = np.roll(phiw, 1, axis=1)

        dxm = (phi - phie) / dx
        dxp = (phiw - phi) / dx
        dym = (phi - phis) / dy
        dyp = (phin - phi) / dy

        dxmxm = (phi - 2 * phie + phiee) / (dx**2)
        dxpxp = (phiww - 2 * phiw + phi) / (dx**2)
        dxpxm = (phie - 2 * phi + phiw) / (dx**2)

        dymym = (phi - 2 * phis + phiss) / (dy**2)
        dypyp = (phinn - 2 * phin + phi) / (dy**2)
        dypym = (phin - 2 * phi + phis) / (dy**2)

        partA = dxm + .5 * dx * np.minimum(np.abs(dxmxm), np.abs(dxpxm))
        partB = dxp - .5 * dx * np.minimum(np.abs(dxpxp), np.abs(dxpxm))
        partC = dym + .5 * dy * np.minimum(np.abs(dymym), np.abs(dypym))
        partD = dyp - .5 * dy * np.minimum(np.abs(dypyp), np.abs(dypym))

        delp2 = g(partA, partB, partC, partD)
        delm2 = g(partB, partA, partD, partC)

        nabla = 0.5 * (dxm**2 + dxp**2 + dym**2 + dyp**2)
        sphi = phi / np.sqrt(phi**2 + np.sqrt(dx**2 + dy**2) * nabla / 10)
        sphip = np.maximum(sphi, 0)
        sphim = np.minimum(sphi, 0)

        phi = phi - dt0 * (sphip * delp2 + sphim * delm2 - sphi)

    return phi


import random

def FEdensity(phi, eps):
    """
    将水平集函数转换为每个矩形有限元的密度函数。
    phi 是水平集，eps 是填充孔洞的非常轻的替代材料的密度。
    """
    FEtheta = np.zeros((nely, nelx))

    # 消除 phi 中值为 0 的节点（我们的算法不能很好地处理这种情况）
    e = 0.000001
    indices = np.where(phi == 0)
    phi[indices] += (2 * np.round(np.random.rand(len(indices[0]))) - 1) * e

    # 遍历每个单元上四个节点处的 phi 值
    for i in range(nely):
        for j in range(nelx):
            phis = [phi[i, j], phi[i, j + 1], phi[i + 1, j + 1], phi[i + 1, j]]
            # 如果全为负，则单元有 full(1) 密度，如果全为正，则单元有 eps 密度
            if all(x < 0 for x in phis):
                FEtheta[i, j] = 1
            elif all(x > 0 for x in phis):
                FEtheta[i, j] = eps
            else:
                # 通过切割单元的两条主对角线将单元分割成四个三角形
                # 然后对周围的节点求和，得到第 5 个值，代表单元的中心值
                # We don't want it to be 0, so we bump it positive or negative 
                # according to equal probabilities.
                phis.append(sum(phis))
                if phis[-1] == 0:
                    phis[-1] += (2 * random.randint(0, 1) - 1) * e
                triMat = np.array([
                    [phis[4], phis[0], phis[1]],
                    [phis[4], phis[1], phis[2]],
                    [phis[4], phis[2], phis[3]],
                    [phis[4], phis[3], phis[0]]
                ])
                # 遍历每个三角形，添加密度到单元的密度值
                for k in range(4):
                    signs = np.sign(triMat[k, :])
                    sum_signs = np.sum(signs)
                    if sum_signs == -3:
                        FEtheta[i, j] += 0.25
                    elif sum_signs == 3:
                        FEtheta[i, j] += 0.25 * eps
                    else:
                        n = np.where(signs != np.sign(np.sum(triMat[k, :])))[0][0]
                        phicc = triMat[k, (n + 1) % 3]
                        phimid = triMat[k, n]
                        phic = triMat[k, (n - 1) % 3]
                        f1 = phimid / (phimid - phicc)
                        f2 = phimid / (phimid - phic)
                        if sum_signs == 1:
                            FEtheta[i, j] += 0.25 * ((1 - f1 * f2) * eps + f1 * f2)
                        else:
                            FEtheta[i, j] += 0.25 * ((1 - f1 * f2) + f1 * f2 * eps)

    return FEtheta

phi = phi00
# 基于 phi 定义有限元密度
FEtheta = FEdensity(phi, eps)
print("FEtheta:", FEtheta.shape, "\n", FEtheta)

def lk():
    """
    刚度矩阵
    这个矩阵通过符号操作得到，是使用有限元方法分析我们的弹性结构的关键。
    """
    E = 1.0  # 弹性模量
    nu = 0.3  # 泊松比
    k = [1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3*nu/8,
         -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3*nu/8]

    KE = E / (1 - nu**2) * np.array([
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

from scipy.sparse import lil_matrix
# 有限元矩阵
KE = lk() # 刚度矩阵

# 全局刚度矩阵
K = lil_matrix((2 * (nelx + 1) * (nely + 1), 2 * (nelx + 1) * (nely + 1)))
# 应用力矩阵
F = lil_matrix((2 * (nely + 1) * (nelx + 1), 2))
# 位移向量矩阵
U = lil_matrix((2 * (nely + 1) * (nelx + 1), 2))

# 有限差分矩阵，用于速度正则化
K1 = lil_matrix(((nelx + 1) * (nely + 1), (nelx + 1) * (nely + 1)))

# 弹性问题的设置
# 力的配置
forceMatrix = np.array([[1, 0.5, 1, -1]])  # 定义应用力
forceCount = forceMatrix.shape[0]  # 应用力的数量

# 设置力
for force in range(forceCount):
    # 计算力的作用节点索引
    node_index = int((forceMatrix[force, 0] * (nelx + 1) + forceMatrix[force, 1] * (nely + 1)) * 2)
    if forceMatrix[force, 2] == 0:  # 水平力
        F[node_index, force] = forceMatrix[force, 3]
    else:  # 垂直力
        F[node_index + 1, force] = forceMatrix[force, 3]
print("F:", F.shape, "\n", F)

# 固定边界条件
fixeddofs = np.arange(2 * (nely + 1))  # 固定x和y自由度的节点
alldofs = np.arange(2 * (nely + 1) * (nelx + 1))
freedofs = np.setdiff1d(alldofs, fixeddofs)  # 计算自由度

from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu

def FE(FEtheta, KE, K, F, U):
    """
    有限元分析和形状梯度计算
    使用密度函数以及关键的有限元向量进行分析，以获取弹性位移，
    该位移用于计算柔顺度和用于传输水平集方程的速度场。
    """
    # 分析
    for elx in range(1, nelx+1):
        for ely in range(1, nely+1):
            n1 = (nely+1) * (elx-1) + ely
            n2 = (nely+1) * elx + ely
            edof = np.array([2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]) - 1
            K[np.ix_(edof, edof)] += FEtheta[ely-1, elx-1] * KE

    # 求解
    sK = csc_matrix(K[np.ix_(freedofs, freedofs)])
    sKcho = splu(sK)

    for i in range(forceCount):
        force_vector = F[freedofs, i].toarray().flatten('F')  # 转换为 ndarray 并展平
        U[freedofs, i] = sKcho.solve(force_vector)

    U[fixeddofs, :] = 0

    # 计算速度场（弹性能量密度）
    FEAeueu = np.zeros((nely, nelx))
    for ely in range(nely):
        for elx in range(nelx):
            n1 = (nely+1) * (elx) + (ely + 1)
            n2 = (nely+1) * (elx + 1) + (ely + 1)
            edof = np.array([2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]) - 1
            for i in range(forceCount):
                Ue = U[edof, i]
                FEAeueu[ely, elx] += FEtheta[ely, elx] * Ue.T @ KE @ Ue

    # 弹性能量密度已经在每个有限元中计算出来，但我们需要将其在传输水平集方程中进一步使用，
    # 所以我们将每个节点周围的有限元的弹性能量密度平均化。
    lvlAeueu = np.zeros((nely+1, nelx+1))
    for i in range(nely+1):
        for j in range(nelx+1):
            bl = FEAeueu[max(0, i-1), max(0, j-1)]
            br = FEAeueu[min(nely-1, i), max(0, j-1)]
            tr = FEAeueu[min(nely-1, i), min(nelx-1, j)]
            tl = FEAeueu[max(0, i-1), min(nelx-1, j)]
            lvlAeueu[i, j] = (bl + br + tr + tl) / 4

    return lvlAeueu, FEAeueu

if ntop == 0:
    # 形状导数
    # 有限元分析的输出是弹性能量密度：
    # lvlAeueu 定义在节点上（用作传输水平集方程的速度），
    # FEAeueu 定义在元素上（用来计算符合性）。
    lvlAeueu, FEAeueu = FE(FEtheta, KE, K, F, U)
else:
    pass

print("lvlAeueu:", lvlAeueu.shape, "\n", lvlAeueu)
print("FEAeueu:", FEAeueu.shape, "\n", FEAeueu)
    
if ntop == 0:
    # 形状梯度
    V = lvlAeueu / (dx*dy) - lagV
    # 速度场的正则化

