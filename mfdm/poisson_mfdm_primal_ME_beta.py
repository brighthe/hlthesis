import numpy as np

from fealpy.decorator import cartesian
from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic

pde = SinSinData()
ns = 2
#mesh = pde.polygon_mesh()
mesh = pde.polygon_mesh_2(n=ns)
import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes, cellcolor='w')
mesh.find_node(axes, showindex=True, 
               color='k', marker='o', markersize=8, fontsize=10, fontcolor='r')
mesh.find_cell(axes, showindex=True, 
               color='k', marker='o', markersize=8, fontsize=10, fontcolor='g')
mesh.find_edge(axes, showindex=True, 
               color='k', marker='o', markersize=8, fontsize=10, fontcolor='b')
plt.show()

maxit = 1
nDof = np.zeros(maxit, dtype=np.int_)

for iter in range(maxit):
    NC = mesh.number_of_cells()
    print("NC:", NC)
    NE = mesh.number_of_edges()
    print("NE:", NE)
    NN = mesh.number_of_nodes()
    print("NN:", NN)

    solver = Mimetic(mesh)

    print("edge:\n", mesh.entity('edge'))
    ME, LHS, error = solver.gme()


    if iter < maxit-1:
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    nDof[iter] = NN

