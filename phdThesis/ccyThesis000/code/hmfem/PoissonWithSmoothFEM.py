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
import sys

degree = int(sys.argv[1])
maxit = 4

x = sp.symbols("x")
y = sp.symbols("y")
#u = sp.sin(x)**3*sp.sin(y)**3 #example
u = (sp.sin(sp.pi*y)*sp.sin(sp.pi*x))**2
#u = x**3*y**4
pde = LaplacePDE(u)

mesh = TriangleMesh.from_box([0, 1, 0, 1], 4, 4)

errorType = ['$|| u - u_h||_{\Omega,0}$',
             '$||\\nabla u - \\nabla u_h||_{\Omega, 0}$'
             ]
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
NDof = np.zeros(maxit, dtype=np.int_)

for i in range(maxit):
    print("The {}-th computation:".format(i))
    space = CmConformingFiniteElementSpace2d(mesh, p=degree, m=1)
    #space = LagrangeFiniteElementSpace(mesh, p=degree)
    NDof[i] = space.dof.number_of_global_dofs()
    #bc1 = DirichletBC(space, [pde.dirichlet, pde.dirichlet]) 
    bc1 = DirichletBC(space, [pde.dirichlet]) 

    uh = space.function() # uh 即是一个有限元函数，也是一个数组
    A = space.grad_m_matrix(1) # (\nabla uh, \nabla v_h)
    #A = space.stiff_matrix()
    F = space.source_vector(pde.source) # (f, vh)

    A, F = bc1.apply(A, F, uh) # 处理 Dirichlet 条件
    uh[:] = spsolve(A, F)

    errorMatrix[0, i] = mesh.error(pde.solution, uh)
    errorMatrix[1, i] = mesh.error(pde.gradient, uh.grad_value)

    if i < maxit-1:
        mesh.uniform_refine()

showmultirate(plt, 2, NDof, errorMatrix,  errorType, 
        propsize=20)

show_error_table(NDof, errorType, errorMatrix)

plt.show()
