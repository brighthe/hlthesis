#!/usr/bin/env python3
# 

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
# solver
from scipy.sparse.linalg import spsolve

from matplotlib import rc
rc('text', usetex=True)
from mpl_toolkits.mplot3d import Axes3D

from pde import TripleLaplacePDE, DoubleLaplacePDE, LaplacePDE 
from cm_conforming_fespace import CmConformingFiniteElementSpace2d
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC 
from fealpy.tools.show import showmultirate
from fealpy.tools.show import show_error_table
from fealpy.mesh import TriangleMesh
from fealpy.decorator import barycentric
import sys

degree = int(sys.argv[1])
maxit = 4

x = sp.symbols("x")
y = sp.symbols("y")
#u = sp.sin(x)**3*sp.sin(y)**3 #example
u = (sp.sin(4*sp.pi*x)*sp.sin(4*sp.pi*y))**5
#u = x**3*y**4

pde = TripleLaplacePDE(u)

mesh = TriangleMesh.from_box([0, 1, 0, 1], 4, 4)
#mesh = TriangleMesh.from_polygon_gmsh([[0, 0], [1, 0], [1, 1], [0, 1]], 0.2)

errorType = [
             '$|| u - u_h||_{\Omega,0}$',
             '$||\\nabla u - \\nabla u_h||_{\Omega, 0}$'
             ]
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
NDof = np.zeros(maxit, dtype=np.int_)

for i in range(maxit):
    print("The {}-th computation:".format(i))
    space = CmConformingFiniteElementSpace2d(mesh, p=degree, m=2)
    NDof[i] = space.dof.number_of_global_dofs()
    bc1 = DirichletBC(space, [pde.dirichlet, pde.dirichlet, pde.dirichlet]) 

    uh = space.function() # uh 即是一个有限元函数，也是一个数组
    A = space.grad_m_matrix(3) # (\nabla uh, \nabla v_h)
    M = space.mass_matrix()
    F = space.source_vector(pde.source) # (f, vh)

    A, F = bc1.apply(A, F, uh) # 处理 Dirichlet 条件
    uh[:] = spsolve(A, F)

    @barycentric
    def ug2val(p):
        return uh.grad_m_value(p, 3)

    errorMatrix[0, i] = mesh.error(pde.solution, uh)
    errorMatrix[1, i] = mesh.error(pde.grad_3, ug2val)

    import copy
    m = copy.deepcopy(mesh)
    node = m.entity('node')
    cell = m.entity('cell')
    bcs = np.eye(3)
    val = uh.grad_m_value(bcs, 2)[..., 0]
    nval = np.zeros(m.number_of_nodes(), dtype=np.float_)
    nval[cell] = val.T

    NN = mesh.number_of_nodes()
    nval = uh[1:NN*6:6]

    node = np.c_[node, nval[:, None]]
    m.node = node
    m.to_vtk(fname='out'+str(i)+'.vtu')

    m = copy.deepcopy(mesh)
    node = m.entity('node')
    cell = m.entity('cell')
    node = np.c_[node, pde.gradient(node)[:, None, 0]]
    m.node = node
    m.to_vtk(fname='tru'+str(i)+'.vtu')

    if i < maxit-1:
        mesh.uniform_refine()

showmultirate(plt, 2, NDof, errorMatrix,  errorType, 
        propsize=20)

show_error_table(NDof, errorType, errorMatrix)

plt.show()
