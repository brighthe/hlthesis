import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData

pde = ClassicalLsfData()
#mesh = pde.polygon_mesh()
mesh = pde.polygon_mesh_2(n=10)
#mesh = pde.triangle_mesh(domain=[0, 1, 0, 1], nx=10, ny=10)

#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=12, fontsize=20, fontcolor='r')
#mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=12, fontsize=20, fontcolor='b')
#mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=12, fontsize=20, fontcolor='g')
#plt.show()

NN = mesh.number_of_nodes()
print("NN:", NN)
NE = mesh.number_of_edges()
print("NE:", NE)
NC = mesh.number_of_cells()
print("NC:", NC)

solver = Mimetic(mesh)
M_c = solver.M_c()
print("M_c:", M_c.shape, "\n", M_c)
M_f = solver.M_f()
print("M_f:", M_f.shape, "\n", M_f)
div_operator = solver.div_operator()
print("div_operator:", div_operator.shape, "\n", div_operator)
gradh = solver.gard_operator()
#print("gradh:", gradh.shape, "\n", gradh)
# Initialize the level set function $phi0$ and velocity field $u$ on the mesh nodes
node = mesh.entity('node')
print("node:", node)
cells = mesh.entity('cell')
print("cells:", cells)
edge = mesh.entity('edge')
print("edge:", edge)

#edge_tangent = mesh.edge_unit_tangent()
#print("edge_tangent:", edge_tangent.shape, "\n", edge_tangent)
edge_norm = mesh.edge_unit_normal()
print("edge_norm:", edge_norm.shape, "\n", edge_norm)
cell_measure = mesh.entity_measure('cell') 
print("cell_measure:", cell_measure.shape, "\n", cell_measure)

#u_nodes = pde.velocity_field(node)
#print("u_nodes:", u_nodes.shape, "\n", u_nodes)
#u_edges = (u_nodes[edge[:, 0]] + u_nodes[edge[:, 1]]) / 2
#print("u_edges:", u_edges.shape, "\n", u_edges)
#u0 = np.einsum('ij, ij -> i', u_edges, edge_norm)
#print("u0:", u0.shape, "\n", u0)
#
#cell2edge = mesh.ds.cell_to_edge()

A10 = solver.u_M_f(velocity=pde.velocity_field)
#print("uhc_M_values:\n", uhc_M_values)
#A10 = np.zeros((NC, NE))
#for i, (cell_edges, uhc_M) in enumerate(zip(cell2edge, uhc_M_values)):
#    A10[i, cell_edges] = uhc_M
print("A10:", A10.shape, "\n", A10)

#asd
#for i in range(NC):
#    indexi, indexj = np.meshgrid(cell2edge[i], cell2edge[i])
#    print("indexi:\n", indexi)
#    print("indexj:\n", indexj)
#
#asd
#
#
#
## 初始化一个列表来存储每个单元内边上的离散速度函数值
#u_cell_edges = []
#
## 遍历 cells，提取每个单元内边上的离散速度值
#for cell in cells:
#    u_cell = u0[cell]  # 提取这个单元内所有边上的速度值
#    u_cell_edges.append(u_cell)  # 添加到列表中
#
## 打印每个单元的速度函数离散值，查看结果
#for i, u_cell in enumerate(u_cell_edges):
#    print(f"Cell {i}: {u_cell}")
#asd
#
#
#
#asd
#u0_Mf = M_f @ u0
#print("u0_Mf", u0_Mf.shape, "\n", u0_Mf)

#u_cell_edges = np.zeros((NC, NE))
#
## 遍历 cells，填充 u_cell_edges
#for i, cell in enumerate(cells):
#    for edge in cell:
#        u_cell_edges[i, edge] = u0_Mf[edge]
#
#print("u_cell_edges:", u_cell_edges.shape, "\n", u_cell_edges)


phi0 = mesh.integral(pde.circle, q=5, celltype=True) / mesh.entity_measure('cell')
print("phi0:", phi0.shape, "\n", phi0)


#from fealpy.functionspace import LagrangeFESpace
#space = LagrangeFESpace(mesh, p=1)
## Compute phi and the gradient of phi at quadrature points
#qf = mesh.integrator(3)
#bcs, _ = qf.get_quadrature_points_and_weights()
#phi_quad = space.value(uh=phi0, bc=bcs)
#print("phi_quad:", phi_quad.shape, "\n", phi_quad)
#grad_phi_quad = space.grad_value(uh=phi0, bc=bcs)
#print("grad_phi_quad:", grad_phi_quad.shape, "\n", grad_phi_quad)

## Compute the magnitude of the gradient at quadrature points
#magnitude = np.linalg.norm(grad_phi_quad, axis=-1)
#
## Identify points at the interface
#at_interface_mask = np.abs(phi_quad) <= 1e-3
#
## Compute the difference between the magnitude and 1 at the interface
#diff = np.abs(magnitude[at_interface_mask]) - 1
#
#diff_avg = np.mean(diff) if np.any(at_interface_mask) else 0
#diff_max = np.max(diff) if np.any(at_interface_mask) else 0
#print(f"Average diff: {diff_avg:.4f}, Max diff: {diff_max:.4f}")



#b = solver.source_neumann(fun=pde.circle)


cell_centers = mesh.entity_barycenter(etype=2) # (NC, GD)
print("cell_centers:", cell_centers.shape, "\n", cell_centers)
x = cell_centers[:, 0]
y = cell_centers[:, 1]
print("x:", x.shape)
x, y = np.meshgrid(x, y)
print("x:", x.shape)


from fealpy.timeintegratoralg import UniformTimeLine
T = 1
nt = 100
# Generate the uniform timeline
timeline = UniformTimeLine(0, T, nt)
dt = timeline.dt

phi0_values = []  # 用于存储每一步的 phi0 值

for i in range(nt):
    t1 = timeline.next_time_level()
    print("t1=", t1)

    A01 = -div_operator.T @ M_c
    #print("A01", A01.shape, "\n", A01)
    A = np.bmat([ [M_f, A01], [dt*A10, M_c] ])
    #print("A", A.shape, "\n", A)

    b1 = cell_measure * phi0
    gamma = np.zeros(NE)
    b = np.hstack((-gamma, -b1))
    #print("b:", b.shape, "\n", b)

    x = np.linalg.solve(A, b)
    phi0[:] = x[-NC:]
    print("phi0:", phi0.shape, "\n", phi0)

    phi0_values.append(phi0.copy())


    #import os
    #mesh.edgedata['phi'] = phi0
    #mesh.celldata['velocity'] = u0
    #output_dir = './visualization/'
    #file_name_prefix = 'lsf'
    #timestep = i+1
    #fname = os.path.join(output_dir, f'{file_name_prefix}_{timestep:010}.vtu')
    #mesh.to_vtk(fname=fname)


    # Move to the next time level
    timeline.advance()

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import griddata

fig, ax = plt.subplots()

x = np.linspace(cell_centers[:, 0].min(), cell_centers[:, 0].max(), 100)
y = np.linspace(cell_centers[:, 1].min(), cell_centers[:, 1].max(), 100)
X, Y = np.meshgrid(x, y)

def update(frame):
    plt.cla()
    scatter = ax.scatter(cell_centers[:, 0], cell_centers[:, 1], c=phi0_values[frame], cmap='viridis')
    Z = griddata(cell_centers, phi0_values[frame], (X, Y), method='cubic')
    contour = ax.contour(X, Y, Z, levels=[0], colors='red')
    title = ax.set_title(f'Time: {frame*dt:.2f}')
    return scatter, contour,

# 创建动画
ani = FuncAnimation(fig, update, frames=nt, blit=False)

plt.show()
