#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from scipy.sparse import csc_matrix, bmat

from fealpy.mesh import TriangleMesh, PolygonMesh
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包
from fealpy.mesh.halfedge_mesh import HalfEdgeMesh2d

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d
from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve, DoubleBandY 

from pde0 import PDE5, PDE4, PDE3

from mumps import DMumpsContext
from scipy.sparse.linalg import minres, gmres, cg, lgmres

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
    """
    #b, x = minres(A, b, show=True, tol=1e-7)
    #b, x = gmres(A, b, tol=1e-2)
    #b, x = lgmres(A, b, atol=1e-5)
    #return b[:N]+b[N:]*(1j)
    x, _ = lgmres(A, b, atol=1e-6)
    #x, _ = gmres(A, b, atol=1e-5, callback=ff, callback_type='x')
    #x, _ = minres(A, b, show=True, tol=1e-10)
    print(np.max(x))
    print(np.max(np.abs(b - A@x)))
    """
    return b[:N]+b[N:]*(1j)

import psutil
import os
fname = sys.argv[1]

#pde = PDE5(om = 1000, eps=0.01)
#pde = PDE4(level=5)
pde = PDE3(level=2, om = 1000, eps=0.01, interfacetype=fname)

mesh = pde.get_mesh_0(2**6, 2**4)

#interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
#mesh = HalfEdgeMesh2d.from_interface_cut_box(interface, [0, 4, 0, 1], 111, 40)
fname = 'poly_interface.svg'
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes, aspect=1, linewidths=0.15, cellcolor=[128/256, 230/256, 115/256])
plt.savefig(fname)

space = HCurlVirtualElementSpace2d(mesh)
gdof = space.dof.number_of_global_dofs()

bc = DirichletBC(space, pde.dirichlet) 

M = space.mass_matrix(alpha=pde.beta, dtype=np.complex128)
print("Assmbly Curl Matrix")
C = space.curl_matrix(beta=pde.alpha, dtype=np.complex128)
print("Assmbly the real part of Source Vector")
br = space.source_vector(lambda x : pde.source(x).real, dtype=np.float_)
print("Assmbly the imag part of Source Vector")
bi = space.source_vector(lambda x : pde.source(x).imag, dtype=np.float_)
b = br+bi*(1j)
B = C-M 

B, b = bc.apply(B, b)
print("Solve the linear equations")
uh = Solve(B, b)

print("Write to VTK")
fname = 'five'
cuh = space.get_curl_value_on_cell(uh.real)
puh = space.projection_to_smspace(uh.real)
mesh.celldata['cuh'] = cuh
mesh.celldata['puh'] = puh
mesh.to_vtk(fname = fname+'_vem_'+str(mesh.number_of_cells())+'.vtu')
#mesh.to_vtk(fname = 'out.vtu')

