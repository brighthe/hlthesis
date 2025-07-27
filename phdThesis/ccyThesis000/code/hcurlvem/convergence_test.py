#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)

from fealpy.mesh import TriangleMesh, PolygonMesh
from fealpy.mesh.halfedge_mesh import HalfEdgeMesh2d
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d
#from HCurlVirtualElementSpace2d0 import HCurlVirtualElementSpace2d
from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve 

from pde0 import PDE1, PDE2

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

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(alpha=pde.beta)
    C = space.curl_matrix(beta=pde.alpha)
    b = space.source_vector(pde.source)
    B = C-M 
    #B = M

    B, b = bc.apply(B, b)
    F = b.copy()
    uh = Solve(B, b)
    uI = space.interpolation(pde.solution)

    # 计算误差
    err0 = space.L2_error(pde.solution, uh, celltype=True)
    err1 = space.curl_error(pde.curl_solution, uh, celltype=True)

    cellmea = mesh.entity_measure('cell')
    #NDof[i] = np.max(np.sqrt(cellmea))
    NDof[i] = gdof
    cellbar = mesh.entity_barycenter('cell')
    val = pde.interface(cellbar)
    e2c = mesh.ds.edge_to_cell()
    flag = val[e2c[:, 0]]*val[e2c[:, 1]] < 0
    iscutcell = np.zeros_like(cellmea, dtype=np.bool_)
    iscutcell[e2c[flag]] = True

    err00 = np.sqrt(np.sum(err0[iscutcell]/np.sum(cellmea[iscutcell])))
    err11 = np.sqrt(np.sum(err1[iscutcell]/np.sum(cellmea[iscutcell])))

    err0 = np.sqrt(np.sum(err0))
    err1 = np.sqrt(np.sum(err1))
    print("error : ", err0)
    return err00, err11

interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
#pde = PDE2(interface)
pde = PDE1(r0=np.pi/5)

maxit = 5
errorType = ['$|| u - u_h||_{0}$',
             '$||\mathrm{rot} u - \mathrm{rot} u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float64)
#NDof = np.zeros(maxit, dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.int_)
for i in range(maxit):
    print('第', i, '次计算')
    mesh = pde.get_mesh(nx = 2**(i+3), ny = 2**(i+3))
    errorMatrix[0, i], errorMatrix[1, i] = compute(mesh, pde)

    if i == 111:
        fname = 'mesh_' + sys.argv[2] + '.svg'
        fig = plt.figure()
        axes = fig.gca()
        #mesh.add_plot(axes, aspect=1)
        plt.savefig(fname, dpi=400)

showmultirate(plt, 2, NDof, errorMatrix,  errorType, propsize=20)
show_error_table(NDof, errorType, errorMatrix)
#plt.savefig(sys.argv[2]+ '_' + sys.argv[1] +'.png', dpi=400)
plt.show()
