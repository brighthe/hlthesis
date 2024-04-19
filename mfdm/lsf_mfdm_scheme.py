import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData

pde = ClassicalLsfData()
ns=5
mesh = pde.polygon_mesh_2(n=ns)

maxit = 5
errorType = ['$|| u - u_h||_{\\Omega,0}$', 
        '$||\\nabla u - \\nabla u_h||_{\\Omega, 0}$']
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
NDof = np.zeros(maxit, dtype=np.int_)

phi0_values = []  # 用于存储每一步的 phi0 值

from fealpy.timeintegratoralg import UniformTimeLine
T = 1
nt = 10000
timeline = UniformTimeLine(0, T, nt)
dt = timeline.dt
for iter in range(maxit):

    NN = mesh.number_of_nodes()
    print("NN:", NN)
    NE = mesh.number_of_edges()
    print("NE:", NE)
    NC = mesh.number_of_cells()
    print("NC:", NC)

    phi0 = mesh.integral(pde.circle, q=5, celltype=True) / mesh.entity_measure('cell')

    solver = Mimetic(mesh)
    M_c = solver.M_c()
    M_f = solver.M_f()
    div_operator = solver.div_operator()
    gradh = solver.gard_operator()

    node = mesh.entity('node')
    cell = mesh.entity('cell')
    edge = mesh.entity('edge')

    edge_norm = mesh.edge_unit_normal()
    cell_measure = mesh.entity_measure('cell') 

    A10 = solver.u_M_f(velocity=pde.velocity_field)

    A01 = -div_operator.T @ M_c
    A = np.bmat([ [M_f, -A01], [-dt*A10, M_c] ])

    b1 = cell_measure * phi0
    gamma = np.zeros(NE)
    b = np.hstack((-gamma, b1))

    x = np.linalg.solve(A, b)
    phi0[:] = x[-NC:]

    phi0_values.append(phi0.copy())

    if iter < maxit-1:
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    NDof[iter] = NN

phi0_differences = []  # 用于存储每次网格细化后的 phi0 差异

for i in range(1, len(phi0_values)):
    difference = np.linalg.norm(phi0_values[i] - phi0_values[i-1], ord=2)  # 计算二范数
    phi0_differences.append(difference)
print("phi0_values:", phi0_values)

import matplotlib.pyplot as plt
from fealpy.tools.show import showmultirate

showmultirate(plt, 2, NDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()
