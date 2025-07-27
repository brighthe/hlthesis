#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt

from fealpy.mesh import TriangleMesh, PolygonMesh
from fealpy.mesh.halfedge_mesh import HalfEdgeMesh2d
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d
from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve 

from pde0 import PDE1, PDE2
from fealpy.pde.MaxwellPDE_2d import MaxwellPDE2d

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

def compute(mesh, pde, k):
    space = HCurlVirtualElementSpace2d(mesh)
    gdof = space.dof.number_of_global_dofs()
    NDof[i] = gdof

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(alpha=pde.beta)
    C = space.curl_matrix(beta=pde.alpha)
    b = space.source_vector(pde.source)
    B = C-M 
    #B = M

    B, b = bc.apply(B, b)
    uh = Solve(B, b)
    uI = space.interpolation(pde.solution)
    # 计算误差

    err0 = space.L2_error(pde.solution, uh, celltype=True)
    err1 = space.curl_error(pde.curl_solution, uh, celltype=True)

    el = mesh.entity_measure('edge')
    et = uh[:, None]*mesh.edge_tangent()/el[:, None]
    e2c = mesh.ds.edge_to_cell()

    c2e, loc = mesh.ds.cell_to_edge()
    c2e = np.split(c2e, loc[1:-1])
    uhI = np.zeros((mesh.ds.NC, 2), dtype=np.float_)
    for ii in range(len(uhI)):
        uhI[ii] = np.average(et[c2e[ii]], axis=0)

    cm = mesh.entity_measure("cell")
    cbary = space.cellbarycenter
    cuh = space.get_curl_value_on_cell(uh.real)
    puh = space.projection_to_smspace(uh.real)
    mesh.celldata['cuh'] = cuh
    mesh.celldata['puh'] = puh
    mesh.celldata['puherror'] = err0/cm
    mesh.celldata['cuherror'] = err1/cm
    mesh.celldata['u'] = pde.solution(cbary)
    mesh.celldata['uhI'] = uhI
    mesh.celldata['cu'] = pde.curl_solution(cbary)
    mesh.celldata['uerror'] = np.abs(pde.solution(cbary)-puh)

    #err0 = space.L2_error(pde.solution, uI, celltype=True)
    #err1 = space.curl_error(pde.curl_solution, uI, celltype=True)
    cuh = space.get_curl_value_on_cell(uI.real)
    puh = space.projection_to_smspace(uI.real)
    mesh.celldata['cuI'] = cuh
    mesh.celldata['puI'] = puh
    #mesh.celldata['puIerror'] = err0/cm
    #mesh.celldata['cuIerror'] = err1/cm
    mesh.to_vtk(fname = 'convergnece_out_'+str(k)+'.vtu')
    k+=1

    err0 = np.sqrt(np.sum(err0))
    err1 = np.sqrt(np.sum(err1))
    print(err0)
    return err0, err1

def make_mesh(eps):
    node = np.array([[0, 0], [1, 0], [1, 1-eps], [1, 1], [1-eps, 1], [0, 1]],
            dtype=np.float_)
    cell = np.array([2, 3, 4, 0, 1, 2, 4, 5], dtype=np.int_)
    cellloc = np.array([0, 3, 8], dtype=np.int_)
    mesh = PolygonMesh(node, cell, cellloc)
    mesh = HalfEdgeMesh2d.from_mesh(mesh)
    return mesh


pde = PDE2()

maxit = 9
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.int_)
for i in range(maxit):
    print('第', i, '次计算')
    mesh = make_mesh(1/2**(i+3))
    errorMatrix[0, i], errorMatrix[1, i] = compute(mesh, pde, i)

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
