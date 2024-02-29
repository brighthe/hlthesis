import numpy as np

from plsf_top_rbfs_cantilever import TopRBFPlsm

# Cantilever
nelx = 6
nely = 3
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

from beam_operator_integrator import BeamOperatorIntegrator
E0 = 1.0
nu = 0.3
integrator = BeamOperatorIntegrator(nu=nu, E0=E0)
KE = integrator.stiff_matrix()
# 每个单元左上角节点的全局编号
eleN1 = np.tile(np.arange(nely).reshape(nely, 1), (1, nelx)) +\
        np.kron(np.arange(nelx), (nely + 1) * np.ones((nely, 1)))
print("eleN1:", eleN1.shape, "\n", eleN1)
# 每个单元上四个节点的全局编号
eleNode = np.tile(eleN1.flatten('F')[:, np.newaxis], (1, 4)) +\
          np.tile(np.array([0, nely+1, nely+2, 1]), (nelx * nely, 1))
eleNode = eleNode.astype(int)
print("eleNode:", eleNode.shape, "\n", eleNode)
# 每个单元上八个自由度的全局编号
edofMat = np.kron(eleNode, np.array([2, 2])) +\
          np.tile(np.array([0, 1]), (nelx*nely, 4))
print("edofMat:", edofMat.shape, "\n", edofMat)


from fealpy.functionspace import LagrangeFESpace as Space
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
# 边界条件定义 - Cantilever
F = np.zeros(vgdof) # 节点荷载
F[2*((nely+1)*nelx+int(np.ceil(nely/2)))] = -100
print("F:", F.shape, "\n", F)
fixeddofs = np.arange(0, 2*(nely+1), 1) # 位移约束
print("fixeddofs:", fixeddofs.shape, "\n", fixeddofs)

# 迭代优化
nLoop = 1 # 优化的最大迭代次数
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
    print("s:", s.shape, "\n", s.round(4))
    print("t:", t.shape, "\n", t.round(4))
    c0 = ((1-s.flatten('F')) * (1-t.flatten('F')))[:, np.newaxis]
    c1 = ((1+s.flatten('F')) * (1-t.flatten('F')))[:, np.newaxis]
    c2 = ((1+s.flatten('F')) * (1+t.flatten('F')))[:, np.newaxis]
    c3 = ((1-s.flatten('F')) * (1+t.flatten('F')))[:, np.newaxis]
    tmpPhi = c0 / 4*Phi.flatten('F')[eleNode[:, 0]] + \
             c1 / 4*Phi.flatten('F')[eleNode[:, 1]] + \
             c2 / 4*Phi.flatten('F')[eleNode[:, 2]] + \
             c3 / 4*Phi.flatten('F')[eleNode[:, 3]]
    print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))
            
    eleVol = np.sum(tmpPhi>=0,axis=1)/np.size(s)



import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
plt.show()
