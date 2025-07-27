

#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.vector import CoordSys3D, Del, curl

from fealpy.functionspace import FirstNedelecFiniteElementSpace2d
from SecondNedelecFiniteElementSpace2d import SecondNedelecFiniteElementSpace2d 
from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d

#from pde import PDE0, PDE1, PDE2, PDE3, boxmesh2d, plusmesh
from pde0 import PDE0, PDE1, PDE6

#from mumps import DMumpsContext
from scipy.sparse.linalg import minres, gmres, cg

def Solve(A, b):
    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_centralized_sparse(A)

    ctx.set_rhs(b)
    ctx.run(job=6)
    ctx.destroy() # Cleanup
    #x, _ = minres(A, b, x0=b, tol=1e-10)
    #x, _ = gmres(A, b, tol=1e-10)
    return b

def compute(mesh, pde):
    space = FirstNedelecFiniteElementSpace2d(mesh)
    #space = SecondNedelecFiniteElementSpace2d(mesh, 1)
    gdof = space.dof.number_of_global_dofs()
    NDof[i] = gdof

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(c=pde.beta)
    C = space.curl_matrix(c=pde.alpha)
    b = space.source_vector(pde.source)

    bc = DirichletBC(space, pde.dirichlet) 
    B = C-M 
    B, b = bc.apply(B, b)
    uh = space.function()
    uh[:] = spsolve(B, b)
    #uh = Solve(B, b)
    # 计算误差

    err0 = space.error(pde.solution, uh, celltype=True)
    err1 = space.curl_error(pde.curl_solution, uh.rot_value, celltype=True)

    cellmea = mesh.entity_measure('cell')
    #NDof[i] = np.max(np.sqrt(cellmea))
    NDof[i] = gdof
    cellbar = mesh.entity_barycenter('cell')
    val = pde.interface(cellbar)
    e2c = mesh.ds.edge_to_cell()
    flag = val[e2c[:, 0]]*val[e2c[:, 1]] < 0
    iscutcell = np.zeros_like(cellmea, dtype=np.bool_)
    iscutcell[e2c[flag]] = True
    iscutcell[[0, 1]] = False
    print(np.where(cellbar[iscutcell][:, 0] > 1/2**(i+3)))
    print(np.where(cellbar[iscutcell][:, 0] < 0))
    print(2**(i+3), iscutcell.sum())

    err00 = np.sqrt(np.sum(err0[iscutcell]/np.sum(cellmea[iscutcell])))
    err11 = np.sqrt(np.sum(err1[iscutcell]/np.sum(cellmea[iscutcell])))

    mesh.celldata['uh'] = uh(np.array([[1/3, 1/3, 1/3]]))[0]
    mesh.celldata['u'] = pde.solution(mesh.entity_barycenter('cell'))
    mesh.to_vtk(fname='out.vtu')

    #eps = pde.eps
    #ff = lambda x, y : np.sum(pde.solution(x)*pde.solution(x), axis=-1)
    #norm = space.integralalg.integral(ff, celltype=True)
    #cellbar = mesh.entity_barycenter('cell')
    #flag = (cellbar[:, 0]<eps) & (cellbar[:, 0]>0)

    #err2 = np.sqrt(np.sum(err0[flag]))/np.sqrt(np.sum(norm[flag]))
    #print(err2)

    err0 = space.error(pde.solution, uh, celltype=True)
    err1 = space.curl_error(pde.curl_solution, uh.rot_value, celltype=True)
    err0 = np.sqrt(np.sum(err0))
    err1 = np.sqrt(np.sum(err1))
    return err00, err11

s = float(sys.argv[1])-0.5
pde = PDE0(s, eps=1e-4)
#pde = PDE6()
print(pde.beta)

maxit = 5
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.int_)
for i in range(maxit):
    print('第', i, '次计算')
    mesh = pde.get_mesh(nx = 2**(i+3), ny = 2**(i+3))
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
