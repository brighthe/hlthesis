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
    nDdof = mesh.entity('edge')[eDdof][:, 0]

    solver = Mimetic(mesh)

    m_v = solver.gmv()
    #print("M_v:", m_v.shape, "\n", m_v)

    m_e = solver.gme()
    #print("m_e:", m_e.shape, "\n", m_e)

    grad_h = solver.grad_operator()
    #print("grad_h:", grad_h.shape, "\n", grad_h)

    A = grad_h.T @ m_e @ grad_h
    print("A:", A.shape, "\n", A)

    rhs = solver.source_primal(fun=pde.source, gddof=nDdof, D=pde.Dirichlet)
    b = m_v @ rhs
    print("b:", b.shape, "\n", b)

    p = np.linalg.solve(A, b)
    print("p:", p.shape, "\n", p)

    node = mesh.entity('node')
    ph = pde.solution(node)
    print("ph:", ph.shape, "\n", ph)

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
