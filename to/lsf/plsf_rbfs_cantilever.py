import numpy as np

from top_plsm_rbfs import TopPlsmRBFs

# Cantilever
nelx = 60
nely = 30
volfrac = 0.5
ts = TopPlsmRBFs(nelx=nelx, nely=nely, volfrac=volfrac)
# 初始化优化参数
nelx, nely, volfrac = ts._nelx, ts._nely, ts._volfrac
mesh = ts._mesh

#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, markersize=6, fontsize=6, fontcolor='r')
#mesh.find_cell(axes, showindex=True, markersize=6, fontsize=6, fontcolor='g')
#plt.show()

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

# 水平集函数的初始化
r = nely * 0.1 # 初始孔洞的半径
Phi = ts.lsf_init(mesh = mesh, r=r)
print("Phi:", Phi.shape, "\n", Phi.round(4))

def output(filename, output, node_val=None, cell_val=None):
    import os
    if not os.path.exists(output):
        os.makedirs(output)
    if node_val is not None:
        mesh.nodedata['node_val'] = node_val.flatten('F') # 按列增加
        fname = os.path.join(output, f'{filename}.vtu')
        mesh.to_vtk(fname=fname)
    if cell_val is not None:
        pass

output(filename='Phi0', output='./visulaization/', node_val=Phi)

# 计算 MQ 样条
A, G, pGpX, pGpY = ts.MQ_spline(mesh=mesh)
print("A:", A.shape, "\n", A.round(4))

# 径向基函数初始化
Alpha = ts.rbf_init(G = G, Phi = Phi)
print("Alpha:", Alpha.shape, "\n", Alpha.round(4))

def stiff_matrix(nu, E0):
    """
    具有 unit 杨氏模量的单元刚度矩阵.
    
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
    A11 = np.array([[12,  3, -6, -3],
                    [ 3, 12,  3,  0],
                    [-6,  3, 12, -3],
                    [-3,  0, -3, 12]])
    A12 = np.array([[-6, -3,  0,  3],
                    [-3, -6, -3, -6],
                    [ 0, -3, -6,  3],
                    [ 3, -6,  3, -6]])
    B11 = np.array([[-4,  3, -2,  9],
                    [ 3, -4, -9,  4],
                    [-2, -9, -4, -3],
                    [ 9,  4, -3, -4]])
    B12 = np.array([[ 2, -3,  4, -9],
                    [-3,  2,  9, -2],
                    [ 4,  9,  2,  3],
                    [-9, -2,  3,  2]])

    KE = E0 / (1-nu**2) / 24 * (np.block([[A11, A12], [A12.T, A11]]) +\
                            nu * np.block([[B11, B12], [B12.T, B11]]))

    return  KE

E0 = 1 # solid 材料的杨氏模量
Emin = 1e-9 # void 材料的杨氏模量
nu = 0.3 # 两种材料的泊松比
KE = stiff_matrix(nu = nu, E0 = E0)

# 初始化单元应变能量场(速度场)
eleComp = np.zeros(nelx*nely, )
# struc = np.ones((nely, nelx))

from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 定义荷载 - cantilever
nLoads = 1
F = np.zeros( (vgdof, nLoads) )
nodal_loads_index = 2*( (nely+1)*nelx + int(np.ceil(nely/2)) ) + 1
F[nodal_loads_index, 0] = -100

# 位移约束(supports) - cantilever
fixeddofs = np.arange(0, 2*(nely+1), 1)

# 迭代优化
nLoop = 200 # 优化的最大迭代次数
nRelax = 30
dt = 0.5 # 水平集演化的时间步长
delta = 10
# Lagrangian 更新方案的参数
mu = 20
gamma = 0.05
lag = 1
# 初始化目标函数值
comp = np.zeros(nLoop)
# 初始化总体积分数
vol = np.zeros(nLoop)

# 绘制结果图
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

eleNode = cell
for iT in range(nLoop):
    # 有限元分析
    # 构件精细网格 - (21, 21)
    s, t = np.meshgrid(np.arange(-1, 1.1, 0.1), np.arange(-1, 1.1, 0.1))
    #print("s:", s.shape, "\n", s.round(4))
    #print("t:", t.shape, "\n", t.round(4))
    c0 = ((1-s.flatten('F')) * (1-t.flatten('F')))[:, np.newaxis]
    c1 = ((1+s.flatten('F')) * (1-t.flatten('F')))[:, np.newaxis]
    c2 = ((1+s.flatten('F')) * (1+t.flatten('F')))[:, np.newaxis]
    c3 = ((1-s.flatten('F')) * (1+t.flatten('F')))[:, np.newaxis]
    tmpPhi = c0 / 4*Phi.flatten('F')[eleNode[:, 0]] + \
             c1 / 4*Phi.flatten('F')[eleNode[:, 1]] + \
             c2 / 4*Phi.flatten('F')[eleNode[:, 2]] + \
             c3 / 4*Phi.flatten('F')[eleNode[:, 3]]
    print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))

    '''
    这里数学上是取 tmpPhi >= 0，但由于不同编程语言中浮点数运算的精度存在差异，
    例如 matlab 中是 -5.55111512e-17，而 python 中是 5.55111512e-17，
    因此我们需要确定一个接近零的阈值
    '''
    # 计算 solid 部分的体积分数
    threshold = 1e-15
    eleVol = np.sum(tmpPhi >= threshold, axis=0) / np.size(s)
    print("eleVol:", eleVol.shape, "\n", eleVol.round(4))
    vol[iT] = np.sum(eleVol) / (nelx*nely)
    #print("vol:", vol.shape, "\n", vol.round(4))

    # 计算单元刚度矩阵中的等效杨氏模量
    E = Emin + eleVol * (E0 - Emin)
    print("E:", E.shape, "\n", E.round(4))

    # 有限元计算全局位移和局部单元位移
    U, Ue = ts.FE(mesh=mesh, E=E, KE=KE, F=F, fixeddofs=fixeddofs)
    print("U:", U.shape, "\n", U.round(4))
    print("Ue:", Ue.shape, "\n", Ue.round(4))

    eleComp[:] = 0
    for i in range(nLoads):
        # 计算单元应变能量场(速度场)
        temp1 = np.einsum('ij, jk, ki -> i', Ue[:, :, i], KE, Ue[:, :, i].T)
        #print("temp1:", temp1.shape, "\n", temp1.round(4))
        eleComp = np.einsum('c, c -> c', E, temp1)
        print("eleComp:", eleComp.shape, "\n", eleComp.round(4))

    # 计算目标函数值
    comp[iT] = np.sum(eleComp)
    #print("comp:", comp.shape, "\n", comp.round(4))

    # 打印当前迭代的结果
    print(f'Iter: {iT}, Obj.: {comp[iT]:.4f}, Vol.: {vol[iT]:.4f}')

    # 白色区域：void 部分, 黑色区域：solid 部分
    cmap = mcolors.ListedColormap(['white', 'black'])
    bounds = [-np.inf, 0, np.inf]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    fig = plt.figure(1)
    contourf = plt.contourf(Phi, levels=bounds, cmap=cmap, norm=norm)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.axis('off')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-5)

    # 收敛性检查
    if iT > nRelax and np.abs(vol[iT] - volfrac) / volfrac < 1e-3 and \
        np.all( np.abs(comp[iT] - comp[iT-8:iT]) / comp[iT] < 1e-3 ):
        break

    # Lagrange 乘子
    if iT < nRelax:
        lag = mu * (vol[iT] - vol[0] + (vol[0] - volfrac) * (iT + 1) / nRelax)
    else:
        lag = lag + gamma * (vol[iT] - volfrac)
        gamma = min(gamma + 0.05, 5)

    # 水平集函数演化
    gradPhi = np.sqrt( (pGpX @ Alpha) ** 2 + (pGpY @ Alpha) ** 2 )
    #print("gradPhi:", gradPhi.shape, "\n", gradPhi.round(4))

    indexDelta = np.abs(Phi) <= delta
    #indexDelta = np.abs(Phi.flatten('F')) <= delta
    #print("indexDelta:", indexDelta.shape, "\n", indexDelta)
    DeltaPhi = np.zeros_like(Phi)
    DeltaPhi[indexDelta] = 0.75 / delta * (1 - Phi[indexDelta] ** 2 / delta ** 2)
    #print("DeltaPhi:", DeltaPhi.shape, "\n", DeltaPhi.round(5))

    eleComp = eleComp.reshape(nelx, nely).T
    #print("eleComp:", eleComp.shape, "\n", eleComp.round(4))
    eleCompLR = np.hstack([eleComp[:, 0:1], eleComp]) + np.hstack([eleComp, eleComp[:, -1:]])
    #print("eleCompLR:", eleCompLR.shape, "\n", eleCompLR.round(4))

    nodeComp = ( np.vstack([eleCompLR, eleCompLR[-1, :]]) + \
                 np.vstack([eleCompLR[0, :], eleCompLR]) ) / 4
    #print("nodeComp:", nodeComp.shape, "\n", nodeComp.round(4))

    B = ( nodeComp.flatten('F') / np.median(nodeComp) - lag ) * \
            DeltaPhi.flatten('F') * delta / 0.75
    #print("B:", B.shape, "\n", B.round(4))

    B_augmented = np.concatenate([B, np.zeros(3)])
    X = np.linalg.solve(G, B_augmented)
    Alpha += dt * X
    #print("Alpha:", Alpha.shape, "\n", Alpha.round(4))

    condition = (eleVol.flatten('F') < 1) & (eleVol.flatten('F') > 0)
    unique_nodes = np.unique(eleNode[condition])
    mean_gradPhi = np.mean(gradPhi[unique_nodes])
    Alpha /= mean_gradPhi
    #print("Alpha:", Alpha.shape, "\n", Alpha.round(4))

    Phi = (G[:-3, :] @ Alpha).reshape(nelx+1, nely+1).T

    #mesh.nodedata['Phi'] = Phi.flatten('F') # 按列增加
    #import os
    #output = './visulaization/'
    #if not os.path.exists(output):
    #    os.makedirs(output)
    #fname = os.path.join(output, 'Phi.vtu')
    #mesh.to_vtk(fname=fname)


    #strucFULL = (Phi > 0).astype(int)
    #struc = strucFULL[1:-1, 1:-1]
    #print("Phi:", Phi.shape, "\n", Phi.round(4))


plt.ioff()
plt.show()
