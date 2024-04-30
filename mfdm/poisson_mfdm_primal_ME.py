import numpy as np

from fealpy.decorator import cartesian
from poisson_model import SinSinData, SinSin5Data, CosCos5Data
from mimetic_solver import Mimetic

pde = SinSinData()
ns = 2
#mesh = pde.polygon_mesh()
mesh = pde.polygon_mesh_2(n=ns)

@cartesian
def fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    val = x + y
    #val = np.sin(np.pi*x)*np.sin(np.pi*y)

    return val

@cartesian
def grad_fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    #val = np.array([1.0, 1.0])
    val = np.zeros_like(p)
    val[..., 0] = 1
    val[..., 1] = 1
    #pi = np.pi
    #val[..., 0] = pi*np.cos(pi*x)*np.sin(pi*y)
    #val[..., 1] = pi*np.sin(pi*x)*np.cos(pi*y)

    return val

maxit = 5
errorType = ['$|| RHS - LHS ||_{\\Omega,0}$']
errorMatrix = np.zeros((3, maxit), dtype=np.float64)
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
    cell2node = mesh.ds.cell_to_node()
    edge_measure = mesh.entity_measure(etype=1)
    eDdof = mesh.ds.boundary_edge_index()
    nDdof = mesh.entity('edge')[eDdof][:, 0]

    solver = Mimetic(mesh)

    # 验证 tau_e 和 nabla_h 的方向是否一致
    grad_h = solver.grad_operator()
    edge_tagnet = node[edge2node[:, 1]] - node[edge2node[:, 0]] # (NE, GD)
    edge_unit_tagnet = edge_tagnet / edge_measure[:, np.newaxis]
    errorMatrix[0, iter] = np.sum(np.abs(grad_h@node - edge_unit_tagnet))

    RHS = mesh.integral(grad_fun, q=5, celltype=True)
    print("RHS:", RHS.shape, "\n", RHS.round(3))

    ME, LHS, error = solver.gme()
    print("LHS:", LHS.shape, "\n", LHS.round(3))

    errorMatrix[1, iter] = np.sum(np.abs(error))

    errorMatrix[2, iter] = np.sum(np.abs(RHS - LHS))

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
