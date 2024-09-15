import numpy as np
import matplotlib.pyplot as plt
from fealpy.mesh import TriangleMesh
from domain_2d import BoxWithCircleHolesDomain
from scipy.spatial import ConvexHull

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

# 零水平集近似为圆，计算每个零水平集的质心和半径
circles = []
for level_set in zero_level_set:
    hull = ConvexHull(level_set)
    center = np.mean(level_set[hull.vertices], axis=0)
    radius = np.max(np.sqrt(np.sum((level_set[hull.vertices] - center) ** 2, axis=1)))
    circles.append((center[0], center[1], radius))

# 绘制零水平集近似的圆
plt.figure()
for circle in circles:
    circle_patch = plt.Circle((circle[0], circle[1]), circle[2], color='r', fill=False)
    plt.gca().add_patch(circle_patch)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.25, 0.25)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

hmin = 0.0125
hmax = 0.1
domain = BoxWithCircleHolesDomain(box=domain, circles=circles, hmin=hmin, hmax=hmax)

maxiter = 200
mesh = TriangleMesh.from_domain_distmesh(domain, maxit=maxiter)

# 绘制生成的网格
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.show()
