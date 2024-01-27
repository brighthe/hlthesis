import numpy as np
import matplotlib.pyplot as plt 

from poisson_1d import Sin4PiData

from fealpy.mesh.uniform_mesh_1d import UniformMesh1d

from scipy.sparse.linalg import spsolve
from scipy.sparse import diags, csr_matrix, spdiags

# PDE 模型
pde = Sin4PiData()
domain = pde.domain()

# 定义网格
nx = 10
hx = (domain[1] - domain[0])/nx
mesh = UniformMesh1d([0, nx], h=hx, origin=domain[0])

fig = plt.figure(1)
axes = fig.gca() 
mesh.add_plot(axes, nodecolor='r') # 在 mesh_base 里面的 plot.py 定义的
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r') # 在 mesh_base 里面的 plot.py 定义的
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b') # 在 mesh_base 里面的 plot.py 定义的

# 组装 Laplace 算子 ∆u 对应的有限差分离散矩阵
cx = 1 / (hx**2)
NN = mesh.number_of_nodes()
k = np.arange(NN)
A = diags([2*cx], [0], shape=(NN, NN), format='csr')
val = np.broadcast_to(-cx, (NN-1, ))
I = k[1:]
J = k[0:-1]
A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
print("A1:", A.shape, "\n", A.toarray())

# 右端项
node = mesh.node
F = pde.source(node)
print("F1:", F.shape, "\n", F.round(4))

# 处理边界条件
uh = mesh.function(etype='node')
index = pde.is_dirichlet_boundary(node)
uh[index] = pde.dirichlet(node[index])

F = F-A@uh
F[index] = uh[index]
print("F2:", F.shape, "\n", F.round(4))

bdIdx = np.zeros(A.shape[0], dtype=np.int_)
bdIdx[index] = 1
D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
A = D0@A@D0 + D1
print("A2:", A.shape, "\n", A.toarray())

# 求解方程组
uh[:] = spsolve(A, F)
print("uh:", uh.shape, "\n", uh)

# 计算误差
et = ['$||u - u_h||_{\infty}$', '$||u - u_h||_{0}$', '$||u - u_h||_{1}$']
eu = np.zeros(len(et), dtype=np.float64) 
uI = pde.solution(node)
e = uI - uh
eu[0] = np.max(np.abs(e))
eu[1] = np.sqrt(hx * np.sum(e**2))
de = e[1:] - e[:-1]
eu[2] = np.sqrt(np.sum(de**2)/hx + eu[1]**2)
et = np.array(et)
print(np.vstack((et, eu)))

# 测试收敛阶
maxit = 5
em = np.zeros((len(et), maxit), dtype=np.float64)
egradm = np.zeros((len(et), maxit), dtype=np.float64) 
hx_list = []

for i in range(maxit):
    hx = mesh.h
    hx_list.append(hx)
    cx = 1 / (hx**2)
    NN = mesh.number_of_nodes()
    k = np.arange(NN)
    A = diags([2*cx], [0], shape=(NN, NN), format='csr')
    val = np.broadcast_to(-cx, (NN-1, ))
    I = k[1:]
    J = k[0:-1]
    A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
    A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
    node = mesh.node
    F = pde.source(node)
    uh = mesh.function(etype='node')
    index = pde.is_dirichlet_boundary(node)
    uh[index] = pde.dirichlet(node[index])
    F = F-A@uh
    F[index] = uh[index]
    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[index] = 1
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
    A = D0@A@D0 + D1
    uh[:] = spsolve(A, F)

    grad_uh = np.gradient(uh[:], hx, edge_order=1)

    uI = pde.solution(node)
    e = uI - uh
    em[0, i] = np.max(np.abs(e))
    em[1, i] = np.sqrt(hx * np.sum(e**2))
    de = e[1:] - e[:-1]
    em[2, i] = np.sqrt(np.sum(de**2)/hx + em[1, i]**2)

    grad_uI = pde.gradient(node)
    grad_e = grad_uI - grad_uh
    egradm[0, i] = np.max(np.abs(grad_e))
    egradm[1, i] = np.sqrt(hx * np.sum(grad_e**2))
    grad_de = grad_e[1:] - grad_e[:-1]
    egradm[2, i] = np.sqrt(np.sum(grad_de**2)/hx + egradm[1, i]**2)

    if i < maxit:
        mesh.uniform_refine()

print("em:\n", em)
print("em_ratio:\n", em[:, 0:-1]/em[:, 1:])
print("egradm:\n", egradm)
print("egradm_ratio:\n", egradm[:, 0:-1]/egradm[:, 1:])

# 绘制收敛阶
linetype = ['r-*', 'g-o', 'y-+']

fig, axes = plt.subplots(2, 1, figsize=(10, 8))
#fig = plt.figure(5)
#plt.xlabel('h')
#plt.ylabel('error')
#axes = fig.gca()
c = np.polyfit(np.log(hx_list), np.log(em[0]), 1)
axes[0].loglog(hx_list, em[0], linetype[0], label='$||u-u_h||_{\infty} = O(h^{%0.4f})$'%(c[0]))
c = np.polyfit(np.log(hx_list), np.log(em[1]), 1)
axes[0].loglog(hx_list, em[0], linetype[1], label='$||u-u_h||_{0} = O(h^{%0.4f})$'%(c[0]))
c = np.polyfit(np.log(hx_list), np.log(em[2]), 1)
axes[0].loglog(hx_list, em[0], linetype[2], label='$||u-u_h||_{1} = O(h^{%0.4f})$'%(c[0]))
axes[0].set_xlabel('h')
axes[0].set_ylabel('error')
axes[0].legend()

c = np.polyfit(np.log(hx_list), np.log(egradm[0]), 1)
axes[1].loglog(hx_list, egradm[0], linetype[0], label='$||\\nabla u-\\nabla u_h||_{\infty} = O(h^{%0.4f})$'%(c[0]))
c = np.polyfit(np.log(hx_list), np.log(egradm[1]), 1)
axes[1].loglog(hx_list, egradm[0], linetype[1], label='$||\\nabla u- \\nabla u_h||_{0} = O(h^{%0.4f})$'%(c[0]))
c = np.polyfit(np.log(hx_list), np.log(egradm[2]), 1)
axes[1].loglog(hx_list, egradm[0], linetype[2], label='$||\\nabla u- \\nabla u_h||_{1} = O(h^{%0.4f})$'%(c[0]))
axes[1].set_xlabel('h')
axes[1].set_ylabel('error')
axes[1].legend()

plt.tight_layout()
plt.show()
