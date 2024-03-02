import numpy as np

from plsf_top_rbfs_cantilever import TopRBFPlsm

# Cantilever
nelx = 60
nely = 30
volfrac = 0.5
ts = TopRBFPlsm(nelx=nelx, nely=nely, volfrac=volfrac)
# 初始化优化参数
nelx, nely, volfrac = ts._nelx, ts._nely, ts._volfrac
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

Phi = ts.lsf_init(mesh = mesh)
print("Phi:", Phi.shape, "\n", Phi.round(4))

A, G, pGpX, pGpY, Alpha = ts.rbf_init(mesh = mesh, Phi = Phi)
#print("pGpX:", pGpX.shape, "\n", pGpX.round(4))
#print("pGpY:", pGpY.shape, "\n", pGpY.round(4))
#print("Alpha:", Alpha.shape, "\n", Alpha.round(4))

def stiff_matrix(nu):
    """
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

    KE = 1 / (1-nu**2) / 24 * (np.block([[A11, A12], [A12.T, A11]]) +\
                            nu * np.block([[B11, B12], [B12.T, B11]]))

    return  KE
nu = 0.3
KE = stiff_matrix(nu = nu)

## 每个单元左上角节点的全局编号
#eleN1 = np.tile(np.arange(nely).reshape(nely, 1), (1, nelx)) +\
#        np.kron(np.arange(nelx), (nely + 1) * np.ones((nely, 1)))
#print("eleN1:", eleN1.shape, "\n", eleN1)
## 每个单元上四个节点的全局编号
#eleNode = np.tile(eleN1.flatten('F')[:, np.newaxis], (1, 4)) +\
#          np.tile(np.array([0, nely+1, nely+2, 1]), (nelx * nely, 1))
#eleNode = eleNode.astype(int)
#print("eleNode:", eleNode.shape, "\n", eleNode)
## 每个单元上八个自由度的全局编号
#edofMat = np.kron(eleNode, np.array([2, 2])) +\
#          np.tile(np.array([0, 1]), (nelx*nely, 4))
#print("edofMat:", edofMat.shape, "\n", edofMat)
## 计算索引向量
#iK = np.kron(edofMat, np.ones((8, 1))).flatten('C')
#iK = iK.astype(np.int64)
#print("iK:", iK.shape, "\n", iK)
#jK = np.kron(edofMat, np.ones((1, 8))).flatten('C')
#jK = jK.astype(np.int64)
#print("jK:", jK.shape, "\n", jK)


from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 边界条件定义 - Cantilever
F = np.zeros(vgdof) # 节点荷载
nodal_loads_index = 2 * ( (nely+1)*nelx + int(np.ceil(nely/2)) ) + 1
print("nodal_loads_index:", nodal_loads_index)
F[nodal_loads_index] = -100
print("F:", F.shape, "\n", F)
fixeddofs = np.arange(0, 2*(nely+1), 1) # 位移约束
print("fixeddofs:", fixeddofs.shape, "\n", fixeddofs)

#cell2node = vspace[0].cell_to_dof()
#print("cell2node:", cell2node.shape, "\n", cell2node)
eleNode = cell

# 迭代优化
nLoop = 32 # 优化的最大迭代次数
nRelax = 30
dt = 0.5 # 水平集演化的时间步长
delta = 10
# 增广 Lagrangian 更新方案的参数
mu = 20
gamma = 0.05
# 目标函数值
comp = np.zeros(nLoop)
# 总体积分数
vol = np.zeros(nLoop)

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
    #print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))

    #from scipy.io import loadmat
    #data = loadmat('tmpPhi.mat')
    #tmpPhi2 = data['tmpPhi']
    #print("tmpPhi2:", tmpPhi2.shape, "\n", tmpPhi2[:, 0])
    #error = np.sum(np.abs(tmpPhi[:, 0] - tmpPhi2[:, 0]))
    #print("error:", error)

    #index1 = tmpPhi[:, 0] < 0
    #print("index1:", index1.shape, "\n", index1)
    #test1 = np.sum(tmpPhi[:, 0] < 0, axis=0)
    #print("test1:", test1.shape, "\n", test1.round(4))
    #index2 = tmpPhi2[:, 0] < 0
    #print("index2:", index2.shape, "\n", index2)
    #test2 = np.sum(tmpPhi2[:, 0] < 0, axis=0)
    #print("test2:", test2.shape, "\n", test2.round(4))

    '''
    这里数学上是取 tmpPhi >= 0，但由于不同编程语言中浮点数运算的精度存在差异，
    例如 matlab 中是 -5.55111512e-17，而 python 中是 5.55111512e-17，
    因此我们需要确定一个接近零的阈值
    '''
    # Solid 部分的体积分数
    threshold = 1e-15
    eleVol = np.sum(tmpPhi >= threshold, axis=0) / np.size(s)
    #print("eleVol:", eleVol.shape, "\n", eleVol.round(4))
    vol[iT] = np.sum(eleVol) / (nelx*nely)
    #print("vol:", vol.shape, "\n", vol.round(4))

    E0 = 1
    Emin = 1e-9
    # 计算等效杨氏模量
    E = Emin + eleVol * (E0 - Emin)
    #print("E:", E.shape, "\n", E.round(4))
    #sK = np.outer(KE.flatten('F'), E).flatten('F')
    #print("sK:", sK.shape, "\n", sK.round(4))

    #from scipy.sparse import coo_matrix
    #K = coo_matrix((sK, (iK, jK)))
    ## CSR 格式对于算术运算和矩阵向量操作更高效
    #K = K.tocsr()
    #K = (K + K.T) / 2
    #print("K:", K.shape, "\n", K.toarray().round(4))

    U, Ue = ts.FE(mesh=mesh, eleVol=eleVol, KE=KE, F=F, fixeddofs=fixeddofs)
    #print("U:", U.shape, "\n", U.round(4))
    #print("Ue:", Ue.shape, "\n", Ue.round(4))

    temp1 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T)
    #print("temp1:", temp1.shape, "\n", temp1.round(4))
    eleComp = np.einsum('c, c -> c', E, temp1)
    #print("eleComp:", eleComp.shape, "\n", eleComp.round(4))
    # 计算目标函数值
    comp[iT] = np.sum(eleComp)
    #print("comp:", comp.shape, "\n", comp.round(4))

    # 打印当前迭代的结果
    print(f'Iter: {iT}, Obj.: {comp[iT]:.4f}, Vol.: {vol[iT]:.4f}')

    # 绘制结果图
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    # 黑色区域：solid 部分，白色区域：void 部分
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

    #x = np.linspace(0, 30, Phi.shape[1])
    #y = np.linspace(0, 60, Phi.shape[0])
    #X, Y = np.meshgrid(x, y)
    #fig = plt.figure(2)
    #ax = fig.add_subplot(111, projection='3d')
    #surf = ax.plot_surface(X, Y, Phi, cmap='viridis')
    #ax.set_xlim(0, 30)
    #ax.set_ylim(0, 60)
    #ax.set_zlim(-12, 12)
    #ax.view_init(elev=20, azim=-35)

    #plt.figure(3)
    #plt.subplot(2, 1, 1)
    #plt.plot(comp[0:iT], '-')
    #plt.title('Compliance')

    #plt.subplot(2, 1, 2)
    #plt.plot(vol[0:iT], '-')
    #plt.title('Volume fraction')
    #plt.tight_layout() # 调整布局以防止标题相互重叠

    # 收敛性检查
    if iT > nRelax and np.abs(vol[iT] - volfrac) / volfrac < 1e-3 and \
        np.all( np.abs(comp[iT] - comp[iT-8:iT]) / comp[iT] < 1e-3 ):
        break

    # Lagrange 乘子
    if iT < nRelax:
        lag = mu * (vol[iT] - vol[0] + (vol[0] - volfrac) * (iT + 1) / nRelax)
    else:
        lag = lag + gamma * (vol[iT-1] - volfrac)
        gamma = min(gamma + 0.05, 5)
        print("lag:", lag)
        print("gamma:", gamma)
    #print("iT:", iT)
    #print("nRelax:", nRelax)
    #print("test:", vol[iT-1] - vol[0])
    #print("test2:", (vol[0] - volfrac) * (iT + 1) / nRelax)
    #print("lag:", lag)

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
    #print("Phi:", Phi.shape, "\n", Phi.round(4))


plt.ioff()
plt.show()




#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
#mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
#plt.show()
