#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.sparse import csc_matrix, bmat

from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d
from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve 

from pde0 import PDE5, PDE4, PDE3

from mumps import DMumpsContext

fname = 'two_circle'
def Solve(A, b):

    N = A.shape[0]
    A = bmat([[A.real, -A.imag], [A.imag, A.real]])
    b = np.r_[b.real, b.imag]

    print("自由度个数 : ", b.shape[0])
    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_centralized_sparse(A)

    ctx.set_rhs(b)
    ctx.run(job=6)
    ctx.destroy() # Cleanup
    return b[:N]+b[N:]*(1j)

def compute(mesh, pde):
    space = HCurlVirtualElementSpace2d(mesh)
    gdof = space.dof.number_of_global_dofs()

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(alpha=pde.beta, dtype=np.complex128)
    C = space.curl_matrix(beta=pde.alpha, dtype=np.complex128)
    br = space.source_vector(lambda x : pde.source(x).real, dtype=np.float_)
    bi = space.source_vector(lambda x : pde.source(x).imag, dtype=np.float_)
    b = br+bi*(1j)

    B = C-M 

    B, b = bc.apply(B, b)
    uh = Solve(B, b)
    #uh = spsolve(B, b)

    #err0 = space.L2_error(pde.solution, uh.real)
    #err1 = space.curl_error(pde.curl_solution, uh.real)
    #print('err0 : ', err0)
    #print('err1 : ', err1)

    cuh = space.get_curl_value_on_cell(uh.real)
    puh = space.projection_to_smspace(uh.real)
    mesh.celldata['cuh'] = cuh
    mesh.celldata['puh'] = puh
    mesh.to_vtk(fname = fname+'_vem_'+str(mesh.number_of_cells())+'.vtu')

    # 画图
    #fig = plt.figure()
    #axes = fig.gca()
    #mesh.add_plot(plt, cellcolor=puh[:, 0], cmap=plt.cm.jet, linewidths=0, 
    #        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    #plt.title('x component of the real part of $u_h^8$')
    #plt.savefig(fname+'_vem_puhx.svg')

    #fig = plt.figure()
    #axes = fig.gca()
    #mesh.add_plot(plt, cellcolor=puh[:, 1], cmap=plt.cm.jet, linewidths=0, 
    #        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    #plt.title('y component of the real part of $u_h^8$')
    #plt.savefig(fname+'_vem_puhy.svg')

    #fig = plt.figure()
    #axes = fig.gca()
    #mesh.add_plot(plt, cellcolor=np.sqrt(puh[:, 0]**2 + puh[:, 1]**2), cmap=plt.cm.jet, linewidths=0, 
    #        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    #plt.title('magnitude the real part of $u_h^8$')
    #plt.savefig(fname+'_vem_puhxy.svg')

    #fig = plt.figure()
    #axes = fig.gca()
    #mesh.add_plot(plt, cellcolor=cuh, cmap=plt.cm.jet, linewidths=0, 
    #        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    #plt.title('rot of the real part of $u_h^8$')
    #plt.savefig(fname+'_vem_cuh.svg')
    #plt.show()
    return space, uh

def compute_error(spaceD, uhD, space, uh, cidx):
    cuhDr = spaceD.get_curl_value_on_cell(uhD.real)
    puhDr = spaceD.projection_to_smspace(uhD.real)
    cuhDi = spaceD.get_curl_value_on_cell(uhD.imag)
    puhDi = spaceD.projection_to_smspace(uhD.imag)

    cuhr = space.get_curl_value_on_cell(uh.real)
    puhr = space.projection_to_smspace(uh.real)
    cuhi = space.get_curl_value_on_cell(uh.imag)
    puhi = space.projection_to_smspace(uh.imag)
    
    cm = spaceD.mesh.entity_measure('cell')
    err0 = np.sqrt(np.sum(cm[:, None]*((puhr[cidx]-puhDr)**2 + (puhi[cidx]-puhDi)**2)))
    err1 = np.sqrt(np.sum(cm*((cuhr[cidx]-cuhDr)**2 + (cuhi[cidx]-cuhDi)**2)))
    return err0, err1

maxit = 5
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.float_)

#pde = PDE5(om = 100, eps=0.01)
#pde = PDE4(level=5)
pde = PDE3(level=6, om = 100, eps=0.01, interfacetype='five_band')
meshD = pde.mesher.mesh0
#meshD = pde.get_mesh(2**10, 2**8)
spaceD, uhD = compute(meshD, pde) 
for i in range(maxit):
    print('第', i, '次计算')
    mesh, cidx = pde.get_mesh(i+1)
    space, uh = compute(mesh, pde)
    NDof[i] = np.max(space.mesh.entity_measure('edge'))

    errorMatrix[0, i], errorMatrix[1, i] = compute_error(spaceD, uhD, space, uh, cidx)

    if i == 1000:
        fname = 'mesh_double_band.svg'
        fig = plt.figure()
        axes = fig.gca()
        mesh.add_plot(axes, aspect=1)
        plt.savefig(fname, dpi=400)


#mesh, _ = pde.get_mesh(2)
#fname = 'mesh_double_band.svg'
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes, aspect=1)
#plt.savefig(fname, dpi=400)

showmultirate(plt, maxit-3, NDof, errorMatrix,  errorType, propsize=20)
show_error_table(NDof, errorType, errorMatrix)
plt.savefig('error.svg', dpi=400)
plt.show()
