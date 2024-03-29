import matplotlib.pyplot as plt

import numpy as np

from poisson import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic


pde = SinSinData()
#pde = SinSin5Data()
#pde = CosCos5Data()
mesh = pde.polygon_mesh_2(n=40)
NC = mesh.number_of_cells()
print("NC:", NC)
NE = mesh.number_of_edges()
print("NE:", NE)
EDdof = mesh.ds.boundary_edge_index()

solver = Mimetic(mesh)
div_operator = solver.div_operator()
print("div_operator:", div_operator.shape, "\n", div_operator)
M_c = solver.M_c()
print("M_c:", M_c.shape, "\n", M_c)
M_f = solver.M_f()
print("M_f:", M_f.shape, "\n", M_f)

b = solver.source(fun=pde.source, gddof=EDdof, D=pde.Dirichlet)
print("b:", b.shape, "\n", b)

A10 = -M_c @ div_operator
A = np.bmat([[M_f, A10.T], [A10, np.zeros((NC, NC))]])
#print("A:", A.shape, "\n", A)

# 单元的积分平均 - (NC, )
p = mesh.integral(pde.solution, q=5, celltype=True) / mesh.entity_measure('cell')

#Ddof = mesh.ds.boundary_cell_flag()
#A, b = solver.boundary_treatment(A,b, Dirichlet, Ddof, so=p)

x = np.linalg.solve(A, b)
#print("x:", x.shape, "\n", x)
ph = x[-NC:]
print("ph:", ph.shape, "\n", ph)

error = p - ph
print(np.max(np.abs(error)))

#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=12, fontsize=20, fontcolor='r')
#mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=12, fontsize=20, fontcolor='b')
#mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=12, fontsize=20, fontcolor='g')
#plt.show()

