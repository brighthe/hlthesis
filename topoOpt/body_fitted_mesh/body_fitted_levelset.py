import numpy as np
import matplotlib.pyplot as plt
from fealpy.mesh import TriangleMesh
from domain_2d import BoxWithZeroLevelSetDomain

nx = 80
ny = 50
domain = [-0.5*nx/100, 0.5*nx/100, -0.5*ny/100, 0.5*ny/100]
mesh = TriangleMesh.from_box(box=domain, nx=nx, ny=ny)
node = mesh.entity('node')

dN = np.sin(6*np.pi*2*node[:, 0]/(domain[1]-domain[0])) * \
     np.cos(6*np.pi*2*node[:, 1]/(domain[1]-domain[0])) + 0.5
xn = node[:, 0].reshape(nx+1, ny+1).T
yn = node[:, 1].reshape(nx+1, ny+1).T
dN = dN.reshape(nx+1, ny+1).T.round(4)

# 节点密度 dN 的零水平集
contour = plt.contour(xn, yn, dN, levels=[0])
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# 获取等高线数据，即具体零水平集值
zero_level_set = [path.vertices for path in contour.collections[0].get_paths()]

# 绘制零水平集数据
plt.figure()
for level_set in zero_level_set:
    plt.plot(level_set[:, 0], level_set[:, 1], 'r-')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Zero Level Set of dN')
plt.show()

hmin = 0.025
hmax = 0.2
domain = BoxWithZeroLevelSetDomain(box=domain, zero_level_set=zero_level_set, \
                                   hmin=hmin, hmax=hmax)

maxiter = 50
mesh = TriangleMesh.from_domain_distmesh(domain, maxit=maxiter)

# 绘制网格
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.show()
