import numpy as np

from fealpy.mesh.polygon_mesh import PolygonMesh
from poisson_model import SinSinData, ExpSinData, CosCosData
from mimetic_solver import Mimetic

#pde = SinSinData()
#pde = ExpSinData()
pde = CosCosData()
ns = 4
mesh = PolygonMesh.from_unit_square(nx=ns, ny=ns)

import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#plt.show()

maxit = 4
errorType = ['$|| p - p_h ||_{\\Omega, 0}$',
             '$|| p - p_h||_{\\Omega, \\infty}$']
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
nDof = np.zeros(maxit, dtype=np.int_)

for iter in range(maxit):
    print("The {}-th computation:".format(iter))
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

    ME = solver.gme()

    grad_h = solver.grad_operator()

    A = grad_h.T @ ME @ grad_h

    node = mesh.entity('node')
    p = pde.solution(node)

    node = mesh.entity('node')
    rhs = pde.source(node)
    ph = np.zeros(NN)

    ph[nDdof] = pde.solution(node[nDdof])
    b = MV @ rhs
    #print("A:", np.sum(np.abs(A)))
    #print("b:", np.sum(np.abs(b)))

    b = b - A @ ph

    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[nDdof.flat] = 1
    from scipy.sparse import spdiags
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0]).toarray()
    A = D0 @ A @ D0 + D1

    b[nDdof] = pde.Dirichlet(node[nDdof])

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

    #print("A:", np.sum(np.abs(A)))
    #print("b:", np.sum(np.abs(b)))
    ph = np.linalg.solve(A, b)
    #print("ph:", np.sum(np.abs(ph)))

    # l2 误差
    errorMatrix[0, iter] = np.sqrt( np.sum( np.abs(p - ph)**2 ) * 1/NN )
    # l_infty 误差
    errorMatrix[1, iter] = np.max( np.abs(p - ph))

    if iter < maxit-1:
        ns = ns*2
        mesh = PolygonMesh.from_unit_square(nx=ns, ny=ns)

    for i, errType in enumerate(errorType):
        print(errType)
        print(errorMatrix[i])
        print('------')
        nDof[iter] = NN

print("errorMatrix:\n", errorMatrix)

from fealpy.tools.show import showmultirate

showmultirate(plt, 2, nDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()
