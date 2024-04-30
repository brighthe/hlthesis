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

@cartesian
def fun(p, index=None):
    x = p[..., 0]
    y = p[..., 1]
    #val = x + y
    val = np.sin(np.pi*x)*np.sin(np.pi*y)

    return val

maxit = 5
errorType = ['$|| RHS - LHS ||_{\\Omega,0}$']
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

    MV, LHS, error = solver.gmv()
    #print("LHS:", LHS.shape, "\n", LHS)
    errorMatrix[0, iter] = np.sum(np.abs(error))
    #print("MV:", MV.shape, "\n", MV.round(3))

    RHS = mesh.integral(fun, q=5, celltype=True)
    #print("RHS:", RHS.shape, "\n", RHS)

    errorMatrix[1, iter] = np.sum(np.abs(RHS - LHS))

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
