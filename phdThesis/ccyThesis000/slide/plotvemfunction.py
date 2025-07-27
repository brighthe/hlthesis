import numpy as np

from fealpy.mesh import TriangleMesh
import gmsh
import numpy as np
import matplotlib.pyplot as plt

from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_mesh(h):
    # 初始化GMSH
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2*h)  # 最大网格尺寸
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)    # 最小网格尺寸


    # 创建一个模型
    gmsh.model.add("pentagon")

    # 设置五边形的参数
    center = [0, 0]  # 五边形中心
    radius = 1.0  # 半径
    num_sides = 5  # 五边形边数

    # 生成五边形的顶点
    angles = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
    points = []
    for i in range(num_sides):
        x = center[0] + radius * np.cos(angles[i])
        y = center[1] + radius * np.sin(angles[i])
        points.append([x, y])
        
    # 创建五边形的顶点并存储点编号
    point_tags = []
    for i, (x, y) in enumerate(points):
        point_tags.append(gmsh.model.geo.addPoint(x, y, 0))

    # 创建五边形的边
    for i in range(num_sides):
        gmsh.model.geo.addLine(point_tags[i], point_tags[(i+1) % num_sides])

    # 创建五边形区域
    gmsh.model.geo.addCurveLoop([i + 1 for i in range(num_sides)], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    # 同步几何信息到GMSH
    gmsh.model.geo.synchronize()

    # 生成网格
    gmsh.model.mesh.generate(2)

    # 导出网格为文件
    #gmsh.write("pentagon.msh")

    # 可视化网格
    #gmsh.fltk.run()

    # get mesh information
    node = gmsh.model.mesh.get_nodes()[1].reshape(-1, 3)[:, :2]
    NN = node.shape[0]

    nid2tag = gmsh.model.mesh.get_nodes()[0] 
    tagmax = int(np.max(nid2tag))
    tag2nid = np.zeros(tagmax+1, dtype = np.int_)
    tag2nid[nid2tag] = np.arange(NN)

    cell = gmsh.model.mesh.get_elements(2, -1)[2][0].reshape(-1, 3)
    cell = tag2nid[cell]

    # Construct FEALPy mesh
    mesh = TriangleMesh(node, cell)

    #mesh.to_vtk(fname='fiveedges.vtu')

    # 关闭GMSH
    gmsh.finalize()
    return mesh

class BasisFunctionOfVEM():
    def __init__(self):
        pass

    def get_mesh(self, h=0.1):
        return get_mesh(h)

    def source(self, p):
        return np.zeros(p[:, 0].shape)

    def dirichlete(self, p):
        flag = p[:, 0] >= np.cos(2*np.pi/5)
        re = np.zeros_like(p[:, 0])
        re[flag] = 1 - np.abs(p[flag, 1] / np.sin(2*np.pi/5))
        return re


from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC

pde = BasisFunctionOfVEM()

mesh = pde.get_mesh(0.02)
space = LagrangeFiniteElementSpace(mesh)
gdof = space.number_of_global_dofs()

A = space.stiff_matrix()
F = np.zeros(gdof, dtype=np.float_)

bc = DirichletBC(space, pde.dirichlete)
A, F = bc.apply(A, F)

uh = space.function()
uh[:] = spsolve(A, F)



def plot_3d_triangular_mesh(nodes, elements):
    """
    使用 matplotlib 的 plot_trisurf 绘制三维三角形网格
    :param nodes: 包含所有节点的坐标，格式为 [(x1, y1, z1), (x2, y2, z2), ...]
    :param elements: 包含所有三角形单元，格式为 [(n1, n2, n3), (n4, n5, n6), ...]
                     每个三角形由三个节点索引组成
    """
    nodes = np.array(nodes)
    elements = np.array(elements)

    # 创建3D绘图
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制三角形网格
    ax.plot_trisurf(nodes[:, 0], nodes[:, 1], nodes[:, 2], triangles=elements,
                    cmap='viridis', edgecolor=None, alpha=0.8)
    
    # 设置坐标轴和标题
    ax.set_title("Basis function", fontsize=16)
    ax.set_xlabel("X Coordinate", fontsize=12)
    ax.set_ylabel("Y Coordinate", fontsize=12)
    ax.set_zlabel("Z Coordinate", fontsize=12)
    
    plt.show()


node = mesh.entity('node')
cell = mesh.entity('cell')
node = np.concatenate([node, uh.reshape(-1, 1)], axis=1)

plot_3d_triangular_mesh(node, cell)













