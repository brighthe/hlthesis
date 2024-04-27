import numpy as np

from fealpy.decorator import cartesian
from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic

pde = SinSinData()
#pde = SinSin5Data()
#pde = CosCos5Data()
ns = 2
#mesh = pde.polygon_mesh()
mesh = pde.polygon_mesh_2(n=ns)
import matplotlib.pyplot as plt

#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#plt.show()

@cartesian
def fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    val = x + y

    return val


maxit = 3
errorType = ['$|| p - p_h||_{\\Omega,0}$']
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
    print("nDdof:", nDdof)

    solver = Mimetic(mesh)

    MV, t2 = solver.gmv()
    print("t2:", t2.shape, "\n", t2)
    #print("MV:", MV.shape, "\n", MV.round(3))

    ME = solver.gme()

    grad_h = solver.grad_operator()

    t1 = mesh.integral(fun, q=5, celltype=True)
    print("t1:", t1.shape, "\n", t1)
    print("error:", np.sum(np.abs(t1 - t2)))

    A = grad_h.T @ ME @ grad_h
    print("A:", A.shape, "\n", A.round(3))

    node = mesh.entity('node')
    p = pde.solution(node)
    print("p:", p.shape, "\n", p.round(3))

    node = mesh.entity('node')
    rhs = pde.source(node)
    ph = np.zeros(NN)
    ph[nDdof] = pde.solution(node[nDdof])
    b = MV @ rhs
    b = b - A @ ph
    print("b:", b.shape, "\n", b.round(3))

    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[nDdof.flat] = 1
    from scipy.sparse import spdiags
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    A = D0 @ A @ D0 + D1
    print("A:", A.shape, "\n", A.round(3))


    b[nDdof] = pde.Dirichlet(node[nDdof])
    print("b:", b.shape, "\n", b.round(3))

    ph = np.linalg.solve(A, b)
    print("ph:", ph.shape, "\n", ph.round(3))


    errorMatrix[0, iter] = np.max(np.abs(p - ph))
    errorMatrix[1, iter] = np.max(np.abs(p[nDdof] - ph[nDdof]))
    print("errorMatrix:\n", errorMatrix)

    if iter < maxit-1:
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    nDof[iter] = NN

import matplotlib.pyplot as plt
from fealpy.tools.show import showmultirate

showmultirate(plt, 2, nDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()
