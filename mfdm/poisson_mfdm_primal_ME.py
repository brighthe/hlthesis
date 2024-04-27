import numpy as np

from fealpy.decorator import cartesian
from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic

pde = SinSinData()
ns = 2
mesh = pde.polygon_mesh()
#mesh = pde.polygon_mesh_2(n=ns)

@cartesian
def fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    val = x + y

    return val

@cartesian
def grad_fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    val = np.zeros_like(p)
    val[:, 0] = 1
    val[:, 1] = 1

    return val

maxit = 1
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
    node = mesh.entity('node')
    edge2node = mesh.ds.edge_to_node()
    edge_measure = mesh.entity_measure(etype=1)
    eDdof = mesh.ds.boundary_edge_index()
    nDdof = mesh.entity('edge')[eDdof][:, 0]

    solver = Mimetic(mesh)

    #MV, t2 = solver.gmv()

    RHS = mesh.integral(grad_fun, q=5, celltype=True)
    print("RHS:", RHS.shape, "\n", RHS.round(3))

    ME, LHS = solver.gme()
    print("ME:", ME.shape, "\n", ME.round(3))
    print("LHS:", LHS.shape, "\n", LHS.round(3))

    if iter < maxit-1:
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    nDof[iter] = NN

