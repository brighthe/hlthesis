import numpy as np
import matplotlib.pyplot as plt

from parabolic_2d import SinSin4PiExpData

from scipy.sparse import diags, csr_matrix, spdiags
from scipy.sparse.linalg import spsolve

from fealpy.mesh.uniform_mesh_2d import UniformMesh2d

# PDE 模型
pde = SinSin4PiExpData()

# 空间离散
domain = pde.domain()
nx = 40
ny = 40
hx = (domain[1] - domain[0])/nx
hy = (domain[3] - domain[2])/ny
mesh = UniformMesh2d([0, nx, 0, ny], h=(hx, hy), origin=(domain[0], domain[2]))

# 时间离散
duration = pde.duration()
nt = 6400
tau = (duration[1] - duration[0])/nt 
print("网比r_x:", tau/(hx**2)) 
print("网比r_y:", tau/(hy**2))

# 准备初值
node = mesh.node
uh0 = pde.init_solution(node)
print("初始解:", uh0.shape)

def advance_backward(n):
    """
    @brief 时间步进格式为向后欧拉方法

    @param[in] n int, 表示第 n 个时间步（当前时间步） 
    """
    t = duration[0] + n*tau
    if n == 0:
        return uh0, t
    else:
        rx = tau/mesh.h[0]**2 
        ry = tau/mesh.h[1]**2 
        if rx + ry > 0.5:
            raise ValueError(f"The rx+ry: {rx + ry} should be smaller than 0.5")

        NN = mesh.number_of_nodes()
        n0 = mesh.nx + 1
        n1 = mesh.ny + 1
        k = np.arange(NN).reshape(n0, n1)

        A = diags([1+2*rx+2*ry], 0, shape=(NN, NN), format='csr')

        val = np.broadcast_to(-rx, (NN-n1, ))
        I = k[1:, :].flat
        J = k[0:-1, :].flat
        A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
        A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)

        val = np.broadcast_to(-ry, (NN-n0, ))
        I = k[:, 1:].flat
        J = k[:, 0:-1].flat
        A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
        A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)

        source = lambda p: pde.source(p, t + tau)

        node = mesh.node
        f = source(node)
        f *= tau
        f += uh0
        gD = lambda p: pde.dirichlet(p, t + tau)

        uh = mesh.function('node').reshape(-1)
        f = f.reshape(-1)
        node = mesh.entity('node')
        isBdNode = pde.is_dirichlet_boundary(node)
        isBdNode = isBdNode.reshape(-1)
        uh[isBdNode] = gD(node[isBdNode])

        f -= A @ uh
        f[isBdNode] = uh[isBdNode]

        bdIdx = np.zeros(A.shape[0], dtype=np.int_)
        bdIdx[isBdNode] = 1
        D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
        D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
        A = D0@A@D0 + D1

        uh0.flat = spsolve(A, f)

        uI = lambda p: pde.solution(p, t + tau)
        
        node = mesh.node
        u_exact = uI(node)
        e = u_exact - uh0
        emax = np.max(np.abs(e))
        print(f"the max error is {emax}")

        return uh0, t

# 绘制结果
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')
box = [0, 1, 0, 1, -1, 1] # 图像显示的范围 0 <= x <= 1, 0 <= y <= 1, -1 <= uh <= 1
mesh.show_animation(fig, axes, box, advance_backward, 
                    fname='parabolic_2d_ab.mp4', plot_type='surface', frames=nt + 1)
plt.show()


