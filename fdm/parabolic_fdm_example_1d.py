import numpy as np
import matplotlib.pyplot as plt

from parabolic_1d import Sin4PiExpData

from scipy.sparse import diags, csr_matrix, spdiags
from scipy.sparse.linalg import spsolve

from fealpy.mesh.uniform_mesh_1d import UniformMesh1d

# PDE 模型
pde = Sin4PiExpData()

# 空间离散
domain = pde.domain()
nx = 40
hx = (domain[1] - domain[0])/nx
mesh = UniformMesh1d([0, nx], h=hx, origin=domain[0])

# 时间离散
duration = pde.duration()
nt = 3200
tau = (duration[1] - duration[0])/nt 
print("网比 r:", tau/(hx**2))

# 准备初值
node = mesh.node
uh0 = pde.init_solution(node)

def advance_backward(n):
    """
    @brief 时间步进格式为向后欧拉方法

    @param[in] n int, 表示第 n 个时间步（当前时间步） 
    """
    t = duration[0] + n*tau
    if n == 0:
        return uh0, t
    else:
        r = tau/mesh.h**2 

        NN = mesh.number_of_nodes()
        k = np.arange(NN)
        A = diags([1+2*r], 0, shape=(NN, NN), format='csr')

        val = np.broadcast_to(-r, (NN-1, ))
        I = k[1:]
        J = k[0:-1]
        A += csr_matrix((val, (I, J)), shape=(NN, NN), dtype=mesh.ftype)
        A += csr_matrix((val, (J, I)), shape=(NN, NN), dtype=mesh.ftype)
        source = lambda p: pde.source(p, t + tau)

        node = mesh.node
        f = source(node)
        f *= tau
        f += uh0
        gD = lambda p: pde.dirichlet(p, t + tau)

        uh = mesh.function('node')
        index = pde.is_dirichlet_boundary(node).reshape(-1)
        uh[index] = gD(node[index])

        f -= A @ uh
        f[index] = uh[index]

        bdIdx = np.zeros(A.shape[0], dtype=np.int_)
        bdIdx[index] = 1
        D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
        D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
        A = D0@A@D0 + D1

        uh0[:] = spsolve(A, f)

        uI = lambda p: pde.solution(p, t + tau)
        
        node = mesh.node
        u_exact = uI(node)
        e = u_exact - uh
        emax = np.max(np.abs(e))

        print(f"the max error is {emax}")

        return uh0, t

# 绘制结果
fig, axes = plt.subplots()
box = [0, 1, -1.5, 1.5] # 图像显示的范围 0 <= x <= 1, -1.5 <= y <= 1.5
mesh.show_animation(fig, axes, box, advance_backward, fname='advance_backward.mp4', 
                    frames=nt+1, lw=2, interval=50, linestyle='--', color='blue')
plt.show()

