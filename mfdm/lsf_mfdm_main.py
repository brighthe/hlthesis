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


NN = mesh.number_of_nodes()
print("NN:", NN)
NE = mesh.number_of_edges()
print("NE:", NE)
NC = mesh.number_of_cells()
print("NC:", NC)

solver = Mimetic(mesh)
M_c = solver.M_c()
#print("M_c:", M_c.shape, "\n", M_c)
M_f = solver.M_f()
#print("M_f:", M_f.shape, "\n", M_f)
div_operator = solver.div_operator()
#print("div_operator:", div_operator.shape, "\n", div_operator)
gradh = solver.gard_operator()
#print("gradh:", gradh.shape, "\n", gradh)

node = mesh.entity('node')
#print("node:", node)
cell = mesh.entity('cell')
#print("cell:", cell)
edge = mesh.entity('edge')
#print("edge:", edge)

edge_norm = mesh.edge_unit_normal()
#print("edge_norm:", edge_norm.shape, "\n", edge_norm)
cell_measure = mesh.entity_measure('cell') 
#print("cell_measure:", cell_measure.shape, "\n", cell_measure)

A10 = solver.u_M_f(velocity=pde.velocity_field)
#print("A10:", A10.shape, "\n", A10)

phi0 = mesh.integral(pde.circle, q=5, celltype=True) / mesh.entity_measure('cell')
print("phi0:", phi0.shape, "\n", phi0)

cell_centers = mesh.entity_barycenter(etype=2) # (NC, GD)
#print("cell_centers:", cell_centers.shape, "\n", cell_centers)
x = cell_centers[:, 0]
y = cell_centers[:, 1]
#print("x:", x.shape)
x, y = np.meshgrid(x, y)
#print("x:", x.shape)

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
    A = np.bmat([ [M_f, -A01], [-dt*A10, M_c] ])

    b1 = cell_measure * phi0
    gamma = np.zeros(NE)
    b = np.hstack((-gamma, b1))

    x = np.linalg.solve(A, b)
    phi0[:] = x[-NC:]
    #print("phi0:", phi0.shape, "\n", phi0)

    phi0_values.append(phi0.copy())

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
