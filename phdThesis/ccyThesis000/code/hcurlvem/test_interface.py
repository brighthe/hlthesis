#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.vector import CoordSys3D, Del, curl

from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh, MeshFactory
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d

#from pde import PDE0, PDE1, PDE2, PDE3, boxmesh2d, plusmesh
from pde0 import PDE0, PDE1

from mumps import DMumpsContext
from scipy.sparse.linalg import minres, gmres, cg

def Solve(A, b):
    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_centralized_sparse(A)

    ctx.set_rhs(b)
    ctx.run(job=6)
    ctx.destroy() # Cleanup
    '''
    x, _ = minres(A, b, x0=b, tol=1e-10)
    #x, _ = gmres(A, b, tol=1e-10)
    '''
    return b

def compute(mesh, pde):
    space = HCurlVirtualElementSpace2d(mesh)
    gdof = space.dof.number_of_global_dofs()
    NDof[i] = gdof

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(alpha=pde.beta)
    C = space.curl_matrix(beta=pde.alpha)
    b = space.source_vector(pde.source)
    B = C-M 

    B, b = bc.apply(B, b)
    #uh = spsolve(B, b)
    uh = Solve(B, b)
    # 计算误差

    err0 = space.L2_error(pde.solution, uh, celltype=True)
    err1 = space.curl_error(pde.curl_solution, uh, celltype=True)

    #eps = pde.eps
    #ff = lambda x, y : np.sum(pde.solution(x)*pde.solution(x), axis=-1)
    #norm = space.integralalg.integral(ff, celltype=True)
    #cellbar = mesh.entity_barycenter('cell')
    #flag = (cellbar[:, 0]<eps) & (cellbar[:, 0]>0)

    #err2 = np.sqrt(np.sum(err0[flag]))/np.sqrt(np.sum(norm[flag]))
    #print(err2)

    err0 = np.sqrt(np.sum(err0))
    err1 = np.sqrt(np.sum(err1))
    return err0, err1

s = float(sys.argv[1])-0.5
#pde = PDE0(s, eps=1e-7)
pde = PDE1()

maxit = 7
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.int_)
for i in range(maxit):
    print('第', i, '次计算')
    mesh = pde.get_mesh(nx = 2**(i+2), ny = 2**(i+2))
    errorMatrix[0, i], errorMatrix[1, i] = compute(mesh, pde)

    if i == 1000:
        fname = 'mesh_' + sys.argv[2] + '.svg'
        fig = plt.figure()
        axes = fig.gca()
        mesh.add_plot(axes)
        plt.savefig(fname, dpi=400)
t1 = time.time()

showmultirate(plt, 2, NDof, errorMatrix,  errorType, propsize=20)
show_error_table(NDof, errorType, errorMatrix)
#plt.savefig(sys.argv[2]+ '_' + sys.argv[1] +'.png', dpi=400)
plt.show()
