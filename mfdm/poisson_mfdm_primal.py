import numpy as np

from fealpy.decorator import cartesian
from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic

pde = SinSinData()
#pde = SinSin5Data()
#pde = CosCos5Data()
ns = 2
mesh = pde.polygon_mesh_2(n=ns)

maxit = 5
errorType = ['$|| p - p_h ||_{\\Omega,0}$']
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
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

    MV = solver.gmv()
    print("MV:", MV.shape, "\n", MV)

    ME = solver.gme()
    print("ME:", ME.shape, "\n", ME)

    grad_h = solver.grad_operator()

    A = grad_h.T @ ME @ grad_h

    node = mesh.entity('node')
    p = pde.solution(node)

    node = mesh.entity('node')
    rhs = pde.source(node)
    ph = np.zeros(NN)

    ph[nDdof] = pde.solution(node[nDdof])
    b = MV @ rhs
    b = b - A @ ph
    print("b:", b.shape, "\n", b)

    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[nDdof.flat] = 1
    from scipy.sparse import spdiags
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    A = D0 @ A @ D0 + D1
    print("A:", A.shape, "\n", A)


    b[nDdof] = pde.Dirichlet(node[nDdof])
    print("b:", b.shape, "\n", b.round(3))

    #f = pde.source(node)
    #b = MV@f

    #eDdof = mesh.ds.boundary_edge_index()
    #nDdof = mesh.entity('edge')[eDdof][:, 0]
    #b[nDdof] = pde.solution(node[nDdof])

    #bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    #bdIdx[nDdof.flat] = 1
    #from scipy.sparse import spdiags
    #D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    #D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    #A = D0 @ A + D1


    ph = np.linalg.solve(A, b)
    print("ph:", ph.shape, "\n", ph.round(3))


    errorMatrix[0, iter] = np.sum(np.abs(p - ph)) / NN
    errorMatrix[1, iter] = np.sum(np.abs(p[nDdof] - ph[nDdof])) / NN

    if iter < maxit-1:
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    nDof[iter] = NN

print("errorMatrix:\n", errorMatrix)
import matplotlib.pyplot as plt
from fealpy.tools.show import showmultirate

showmultirate(plt, 2, nDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()
