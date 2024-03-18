import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData


pde = ClassicalLsfData()
mesh = pde.polygon_mesh()

NN = mesh.number_of_nodes()
print("NN:", NN)
NE = mesh.number_of_edges()
print("NE:", NE)

solver = Mimetic(mesh)
gradh = solver.gard_operator()
print("gradh:", gradh.shape, "\n", gradh)

print("--------------", solver.M_f().shape)
# Initialize the level set function $phi0$ and velocity field $u$ on the mesh nodes
node = mesh.entity('node')
edge = mesh.entity('edge')
print("edge:", edge.shape, "\n", edge)

p1 = node[edge[:, 0], :]
p2 = node[edge[:, 1], :]
phi0_edge = np.vstack((pde.circle(p1), pde.circle(p2))).T
print("phi_edge:", phi0_edge.shape, "\n", phi0_edge)

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



