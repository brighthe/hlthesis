import numpy as np
import gmsh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC
from fealpy.mesh import TriangleMesh
from scipy.sparse.linalg import spsolve

edge_num = 3
cell_num = 6

def get_mesh(h):
    # 初始化GMSH
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2*h)  # 最大网格尺寸
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)    # 最小网格尺寸


    # 创建一个模型
    gmsh.model.add("pentagon")

    # 设置五边形的参数
    center = [-1, 0]  # 五边形中心
    radius = 1.0  # 半径
    num_sides = 3  # 五边形边数

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


def get_all_mesh(h):
    # 初始化GMSH
    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2*h)  # 最大网格尺寸
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)    # 最小网格尺寸

    # 创建一个模型
    gmsh.model.add("three_pentagons")

    # 设置五边形的参数
    radius = 1.0  # 半径
    num_sides = 6  # 五边形边数
    centers = np.array([[-1, 0]])  # 三个五边形中心的坐标

    # 生成五边形的顶点
    points = []
    point_tags = []
    for center in centers:
        angles = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
        for i in range(num_sides):
            x = center[0] + radius * np.cos(angles[i])
            y = center[1] + radius * np.sin(angles[i])
            points.append([x, y])
            point_tags.append(gmsh.model.occ.addPoint(x, y, 0))

    R = np.array([[np.cos(2*np.pi / 3), -np.sin(2*np.pi / 3)],
                  [np.sin(2*np.pi / 3), np.cos(2*np.pi / 3)]])

    # 创建五边形的边
    num_points = len(point_tags)
    line_tags = []
    for i in range(0, num_points, num_sides):
        for j in range(num_sides):
            l0 = gmsh.model.occ.addLine(point_tags[i + j], point_tags[i + (j + 1) % num_sides])
            line_tags.append(l0)

    p0 = R@ points[2] 
    p0 = gmsh.model.occ.addPoint(p0[0], p0[1], 0)
    p1 = R @ points[3]
    p1 = gmsh.model.occ.addPoint(p1[0], p1[1], 0)
    p2 = R @ points[4]
    p2 = gmsh.model.occ.addPoint(p2[0], p2[1], 0)
    p3 = R @ points[5]
    p3 = gmsh.model.occ.addPoint(p3[0], p3[1], 0)

    l0 = gmsh.model.occ.addLine(p0, p1)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(6, p0)
    l3 = gmsh.model.occ.addLine(p2, p3)
    l4 = gmsh.model.occ.addLine(1, p3)
    line_tags += [l0, l1, l2, l3, l4]

    p4 = R@R@ points[2] 
    p4 = gmsh.model.occ.addPoint(p4[0], p4[1], 0)
    p5 = R@R @ points[3]
    p5 = gmsh.model.occ.addPoint(p5[0], p5[1], 0)
    p6 = R@R @ points[4]
    p6 = gmsh.model.occ.addPoint(p6[0], p6[1], 0)
    p7 = R@R@points[5]
    p7 = gmsh.model.occ.addPoint(p7[0], p7[1], 0)

    l0 = gmsh.model.occ.addLine(p4, p3)
    l1 = gmsh.model.occ.addLine(p4, p5)
    l2 = gmsh.model.occ.addLine(p5, p6)
    l3 = gmsh.model.occ.addLine(p6, p7)
    line_tags += [l0, l1, l2, l3]

    a = gmsh.model.occ.addRectangle(-8, -8, 0, 16, 16, 0)

    gmsh.model.occ.fragment([(2, a)], [(1, l) for l in line_tags])
        
    # 同步几何信息到GMSH
    gmsh.model.occ.synchronize()

    # 生成网格
    gmsh.model.mesh.generate(2)

    #gmsh.fltk.run()

    # 获取网格节点和元素
    node = gmsh.model.mesh.get_nodes()[1].reshape(-1, 3)[:, :2]
    NN = node.shape[0]

    nid2tag = gmsh.model.mesh.get_nodes()[0]
    tagmax = int(np.max(nid2tag))
    tag2nid = np.zeros(tagmax + 1, dtype=np.int_)
    tag2nid[nid2tag] = np.arange(NN)

    ctags = []
    cell = []
    dimtags = gmsh.model.getEntities(2)
    for dim, tag in dimtags:
        cell0 = gmsh.model.mesh.get_elements(dim, tag)[2][0].reshape(-1, 3)
        cell.append(cell0)
        ctags.append(np.ones(cell0.shape[0], dtype=np.int_)*tag)

    cell = np.concatenate(cell)
    ctags = np.concatenate(ctags)
    cell = tag2nid[cell]

    # Construct FEALPy mesh
    mesh = TriangleMesh(node, cell)
    mesh.celldata['celltype'] = ctags

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
        flag = p[:, 0]+1 >= np.cos(2 * np.pi / 3)+1e-13
        re = np.zeros_like(p[:, 0])
        re[flag] = 1 - np.abs(p[flag, 1] / np.sin(2 * np.pi / 3))
        return re


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


# 使用FEALPy
#pde = BasisFunctionOfVEM()
#mesh = pde.get_mesh(0.02)

# 获取网格信息
#node = mesh.entity('node')
#cell = mesh.entity('cell')

# 可视化网格
#plot_3d_triangular_mesh(node, cell)



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

node = mesh.entity('node')
cell = mesh.entity('cell')
node = np.concatenate([node, uh.reshape(-1, 1)], axis=1)
mesh.nodedata['uh'] = uh
mesh.node = node
print(mesh.node.shape)
mesh.to_vtk(fname='basisfunction0.vtu')


R = np.array([[np.cos(2*np.pi / 6), -np.sin(2*np.pi / 6)],
              [np.sin(2*np.pi / 6), np.cos(2*np.pi / 6)]])

for i in range(1, 6):
    node[:, :2] = node[:, :2]@R.T
    mesh.node = node
    mesh.to_vtk(fname='basisfunction%d.vtu' % i)

mesh1 = get_all_mesh(0.02)
node1 = mesh1.entity('node')
cell1 = mesh1.entity('cell')
mesh1.node = node1
mesh1.nodedata['c'] = np.zeros(node1.shape[0])
mesh1.to_vtk(fname='basisfunctionall.vtu')

plot_3d_triangular_mesh(node, cell)























