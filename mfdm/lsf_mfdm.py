import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData


pde = ClassicalLsfData()
mesh = pde.polygon_mesh()

solver = Mimetic(mesh)
gradh = solver.gard_operator()
print("gradh:", gradh.shape, "\n", gradh)


# Initialize the level set function $phi0$ and velocity field $u$ on the mesh nodes
node = mesh.entity('node')
phi0 = pde.circle(node)
print("phi0:", phi0.shape, "\n", phi0)

edge = mesh.entity('edge')
edge_tangent = node[edge[:, 1], :] - node[edge[:, 0], :]
print("edge_tangent:", edge_tangent.shape, "\n", edge_tangent)
edge_length = mesh.entity_measure('edge')
print("edge_length:", edge_length.shape, "\n", edge_length)
edge_center = mesh.entity_barycenter(etype=1)
print("edge_center:", edge_center.shape, "\n", edge_center)


GD = mesh.geo_dimension()
qf = mesh.integrator(q=5, etype='cell')
bcs, ws = qf.get_quadrature_points_and_weights()
#ps = mesh.bc_to_point(bcs)
u0 = pde.velocity_field(edge_center)
print("u0:", u0.shape, "\n", u0)
u = mesh.integral(u=pde.velocity_field*edge_tangent, q=5, celltype=True) / edge_length
print("u:", u.shape, "\n", u)
#u = space.interpolate(pde.velocity_field, dim=2)

NN = mesh.number_of_nodes()

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


fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=12, fontsize=20, fontcolor='r')
mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=12, fontsize=20, fontcolor='b')
mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=12, fontsize=20, fontcolor='g')
plt.show()


