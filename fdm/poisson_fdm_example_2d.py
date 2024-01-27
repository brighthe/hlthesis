import numpy as np 
import matplotlib.pyplot as plt 

from poisson_2d import SinSinData

from fealpy.mesh.uniform_mesh_2d import UniformMesh2d

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags, csr_matrix, spdiags

# PDE 模型
pde = SinSinData()
domain = pde.domain()

# 定义网格
nx = 5
ny = 5
hx = (domain[1] - domain[0])/nx
hy = (domain[3] - domain[2])/ny
mesh = UniformMesh2d((0, nx, 0, ny), h=(hx, hy), origin=(domain[0], domain[2]))

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_edge(axes, showindex=True, fontsize=12, fontcolor='g') 
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')

# 组装 Laplace 算子 ∆u 对应的有限差分离散矩阵
cx = 1 / (hx**2)
cy = 1 / (hy**2)
n0 = mesh.ds.nx + 1
n1 = mesh.ds.ny + 1
NN = mesh.number_of_nodes()
k = np.arange(NN).reshape(n0, n1)
A = diags([2*(cx+cy)], [0], shape=(NN, NN), format='csr')

val = np.broadcast_to(-cx, (NN-n1, ))
I = k[1:, :].flat
J = k[0:-1, :].flat
A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
val = np.broadcast_to(-cy, (NN-n0, ))
I = k[:, 1:].flat
J = k[:, 0:-1].flat
A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
print("A1:", A.shape, "\n", A.toarray())

# 右端项
node = mesh.entity('node')
F = pde.source(node)
print("F1:", F.shape, "\n", F.round(4))

# 处理边界条件
uh = mesh.function(etype='node').reshape(-1)
print("uh:", uh.shape, "\n", uh)
isBdNode = pde.is_dirichlet_boundary(node)
isBdNode = isBdNode.reshape(-1)
print("isBdNode:", isBdNode.shape, "\n", isBdNode)
uh[isBdNode] = pde.dirichlet(node[isBdNode])

F = F.flat-A@uh
F[isBdNode] = uh[isBdNode]
print("F2:", F.shape, "\n", F.round(4))

bdIdx = np.zeros(A.shape[0], dtype=np.int_)
bdIdx[isBdNode] = 1
D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
A = D0@A@D0 + D1
print("A2:", A.shape, "\n", A.toarray())

# 求解方程组
uh[:] = spsolve(A, F)
print("uh:", uh.shape, "\n", uh)

# 计算误差和测试收敛阶
maxit = 5
et = ['$||u - u_h||_{\infty}$', '$||u - u_h||_{0}$', '$||u - u_h||_{1}$']
em = np.zeros((len(et), maxit), dtype=np.float64)

for i in range(maxit):
    hx = mesh.h[0]
    hy = mesh.h[1]
    cx = 1 / (hx**2)
    cy = 1 / (hy**2)
    n0 = mesh.ds.nx + 1
    n1 = mesh.ds.ny + 1
    NN = mesh.number_of_nodes()
    k = np.arange(NN).reshape(n0, n1)
    A = diags([2*(cx+cy)], [0], shape=(NN, NN), format='csr')
    val = np.broadcast_to(-cx, (NN-n1, ))
    I = k[1:, :].flat
    J = k[0:-1, :].flat
    A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
    A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
    val = np.broadcast_to(-cy, (NN-n0, ))
    I = k[:, 1:].flat
    J = k[:, 0:-1].flat
    A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
    A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
    node = mesh.entity('node')
    F = pde.source(node)
    uh = mesh.function(etype='node').reshape(-1)
    isBdNode = pde.is_dirichlet_boundary(node)
    isBdNode = isBdNode.reshape(-1)
    uh[isBdNode] = pde.dirichlet(node[isBdNode])
    F = F.flat - A@uh
    F[isBdNode] = uh[isBdNode]
    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[isBdNode] = 1
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
    A = D0@A@D0 + D1
    uh[:] = spsolve(A, F)

    nx = mesh.nx
    ny = mesh.ny
    node = mesh.entity('node')
    uI = pde.solution(node)
    e = uI - uh
    em[0, i] = np.max(np.abs(e))
    em[1, i] = np.sqrt(hx * hy * np.sum(e**2))
    em[2, i] = np.sqrt(1 / ((nx - 1) * (ny - 1)) * np.sum(e ** 2))

    if i < maxit:
        mesh.uniform_refine()

print("em:\n", em)
print("em_ratio:\n", em[:, 0:-1]/em[:, 1:])
