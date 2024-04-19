import matplotlib.pyplot as plt

import numpy as np

from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic


pde = SinSinData()
#pde = SinSin5Data()
#pde = CosCos5Data()
ns = 5
mesh = pde.polygon_mesh()
#mesh = pde.polygon_mesh_2(n=ns)

maxit = 1
errorType = ['$|| p - p_h||_{\\Omega,0}$']
errorMatrix = np.zeros((1, maxit), dtype=np.float64)
nDof = np.zeros(maxit, dtype=np.int_)

for iter in range(maxit):
    NC = mesh.number_of_cells()
    print("NC:", NC)
    NE = mesh.number_of_edges()
    print("NE:", NE)
    NN = mesh.number_of_nodes()
    print("NN:", NN)
    eDdof = mesh.ds.boundary_edge_index()

    solver = Mimetic(mesh)
    div_operator = solver.div_operator()
    #print("div_operator:", div_operator.shape, "\n", div_operator)

    M_c = solver.M_c()
    #print("M_c:", M_c.shape, "\n", M_c)
    M_f = solver.M_f()
    #print("M_f:", M_f.shape, "\n", M_f)

    b = solver.source(fun=pde.source, gddof=eDdof, D=pde.Dirichlet)
    #print("b:", b.shape, "\n", b)

    A10 = -M_c @ div_operator
    A = np.bmat([[M_f, A10.T], [A10, np.zeros((NC, NC))]])
    #print("A:", A.shape, "\n", A)

    # 单元的积分平均 - (NC, )
    p = mesh.integral(pde.solution, q=5, celltype=True) / mesh.entity_measure('cell')

    x = np.linalg.solve(A, b)

    ph = x[-NC:]
    #print("ph:", ph.shape, "\n", ph)

    errorMatrix[0, iter] = np.max(np.abs(p - ph))
    print("errorMatrix:", errorMatrix)

    if iter < maxit-1:
        #mesh.uniform_refine()
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    nDof[iter] = NN

import matplotlib.pyplot as plt
from fealpy.tools.show import showmultirate

showmultirate(plt, 2, nDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()

