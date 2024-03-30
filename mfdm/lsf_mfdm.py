import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData

pde = ClassicalLsfData()
mesh = pde.polygon_mesh()

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=12, fontsize=20, fontcolor='r')
mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=12, fontsize=20, fontcolor='b')
mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=12, fontsize=20, fontcolor='g')
plt.show()

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
cell = mesh.entity('cell')
print("cell:", cell)
edge = mesh.entity('edge')

edge_center = mesh.entity_barycenter(etype=1)
print("edge_center:", edge_center.shape, "\n", edge_center)
u0 = pde.velocity_field(edge_center)
print("u0:", u0.shape, "\n", u0)

cell_values = []
for i in cell:
    cell_values.append(u0[i])

u0_cell = cell_values
for i, values in enumerate(u0_cell):
    print(f"u0_cell {i}:")
    print(values.round(4))

A10 = -div_operator.T @ M_c
print("A10", A10.shape)
A = np.bmat([ [M_f, A10], [A10, M_c] ])


asd

p1 = node[edge[:, 0], :]
p2 = node[edge[:, 1], :]
phi0_edge = np.vstack((pde.circle(p1), pde.circle(p2))).T
print("phi_edge:", phi0_edge.shape, "\n", phi0_edge)
asd

edge = mesh.entity('edge')
edge_tangent = node[edge[:, 1], :] - node[edge[:, 0], :]
print("edge_tangent:", edge_tangent.shape, "\n", edge_tangent)
edge_length = mesh.entity_measure('edge')
print("edge_length:", edge_length.shape, "\n", edge_length)
edge_center = mesh.entity_barycenter(etype=1)
print("edge_center:", edge_center.shape, "\n", edge_center)

u0 = pde.velocity_field(edge_center)
print("u0:", u0.shape, "\n", u0)
u = u0 * edge_center
print("u:", u.shape, "\n", u)

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=12, fontsize=20, fontcolor='r')
mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=12, fontsize=20, fontcolor='b')
mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=12, fontsize=20, fontcolor='g')
plt.show()


# Generate the uniform timeline
from fealpy.timeintegratoralg import UniformTimeLine
T = 2
nt = 200
timeline = UniformTimeLine(0, T, nt)
dt = timeline.dt
# Time iteration
for i in range(nt):
    t1 = timeline.next_time_level()
    print("t1=", t1)

    A = diags([1 + dt * np.dot(u, gradh)], 0, shape=(NN, NN), format='csr')
    phi0[:] = spsolve(A, phi0)



