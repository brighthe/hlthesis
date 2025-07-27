import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fealpy.mesh import TriangleMesh

mesh = TriangleMesh.from_one_triangle()
mesh.uniform_refine(5)

nodes = mesh.entity('node')
cells = mesh.entity('cell')
NC = mesh.number_of_cells()


# 创建图形
fig = plt.figure(figsize=(20, 5))

phival = 1 - nodes[:, 0] - nodes[:, 1]
# 绘制第一个基函数
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_trisurf(nodes[:, 0], nodes[:, 1], phival, triangles=cells,
                cmap='viridis', edgecolor=None, alpha=0.8)
ax1.set_title('$\\phi_0$')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('$\\phi_0$')

# 绘制第二个基函数
phival = nodes[:, 0]
ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_trisurf(nodes[:, 0], nodes[:, 1], phival, triangles=cells,
                cmap='viridis', edgecolor=None, alpha=0.8)
ax2.set_title('$\\phi_1$')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('$\\phi_1$')

# 绘制第三个基函数
phival = nodes[:, 1]
ax3 = fig.add_subplot(133, projection='3d')
ax3.plot_trisurf(nodes[:, 0], nodes[:, 1], phival, triangles=cells,
                cmap='viridis', edgecolor=None, alpha=0.8)
ax3.set_title('$\\phi_2$')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax3.set_zlabel('$\\phi_2$')

plt.tight_layout()
plt.savefig('fem_basis_function.pdf', dpi=300)
plt.show()
