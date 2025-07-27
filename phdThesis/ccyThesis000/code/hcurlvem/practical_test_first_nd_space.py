#!/usr/bin/env python3
# 

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.vector import CoordSys3D, Del, curl

from fealpy.geometry import CircleCurve, FoldCurve, Polygon
from interface_mesh_2d import HalfEdgeMesh2dWithInterface
from fealpy.mesh.interface_mesh_generator import interfacemesh2d 

from fealpy.pde.MaxwellPDE_2d import SinData as PDE
from fealpy.functionspace import FirstNedelecFiniteElementSpace2d
from SecondNedelecFiniteElementSpace2d import SecondNedelecFiniteElementSpace2d 
from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh
from fealpy.decorator import cartesian, barycentric
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包
from fealpy.pde.MaxwellPDE_2d import MaxwellPDE2d

from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve, DoubleBandY, Polygon, BandY 

from scipy.sparse import csc_matrix, bmat
from scipy.sparse.linalg import spsolve, cg
from pde0 import PDE5
from mumps import DMumpsContext
def Solve(A, b):
    N = A.shape[0]
    A = bmat([[A.real, -A.imag], [A.imag, A.real]])
    b = np.r_[b.real, b.imag]

    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_centralized_sparse(A)

    ctx.set_rhs(b)
    ctx.run(job=6)
    ctx.destroy() # Cleanup
    return b[:N]+b[N:]*(1j)
"""
def Solve(A, b):
    x, _ = lgmres(A, b, atol=1e-6)
    print(np.max(x))
    print(np.max(np.abs(b - A@x)))
    return x
"""

class PDE2(MaxwellPDE2d):
    def __init__(self):
        self.interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))

        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        f = sym.sin(10*C.y)*C.i + 0.01*sym.sin(10*C.x)*C.j
        super(PDE2, self).__init__(f)

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = self.interface(p) < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = self.interface(p) < 0
        val = 1000*np.ones_like(x)
        val[flag] = val[flag]*1
        return val

    @cartesian
    def source(self, p):
        x = p[..., 0, None]
        y = p[..., 1, None]

        alpha = self.alpha(p)[..., None]
        beta = self.beta(p)[..., None]

        ccFx = self.curlcurlFx(x, y)
        ccFy = self.curlcurlFy(x, y)
        if type(ccFx) is not np.ndarray:
            ccFx = np.ones(x.shape, dtype=np.float_)*ccFx
        if type(ccFy) is not np.ndarray:
            ccFy = np.ones(x.shape, dtype=np.float_)*ccFy
        ccf = np.c_[ccFx, ccFy] 
        return alpha*ccf - beta*self.solution(p)

    def get_mesh(self, nx=1, ny=1):
        mesh = TriangleMesh.interfacemesh_generator([0, 4, 0, 1], nx, ny,
                self.interface)
        return mesh


class PDE3():
    def __init__(self, om = 5, eps = 0.01, sig0 = 1, sig1=0.1, level = 7,
            interfacetype='double_circle'):
        self.om = om
        self.eps = eps
        self.sig0 = sig0
        self.sig1 = sig1

        if interfacetype=='double_circle':
            self.interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
        elif interfacetype=='double_band':
            self.interface = DoubleBandY()
        elif interfacetype=='five_band':
            self.interface = BandY([0.11, 0.13, 0.24, 0.26, 0.49, 0.51, 0.74, 0.76, 0.87, 0.89])
        elif interfacetype=='poly':
            points = np.array([[311+2/3, 50], [290+0.625, 47+2/3], [225, 68+1/39],
             [206.25, 19+1/3], [198.05, 66+2/3], [178.125, 15+1/3], [175+1/3, 62+1/3],
             [100+1/3, 69-1/4], [243.75, 82+1/3], [287.5, 52+1/3]])/100
            points[:, 1] = 1-points[:, 1]
            self.interface = Polygon(points)

    @cartesian
    def source(self, p):
        om = self.om
        eps = self.eps

        x = p[..., 0]
        val = np.zeros_like(p, dtype=np.complex128)
        val[..., 1] = -om*(1j)*np.exp(-((x-3)**2)/(eps**2))
        return val

    @cartesian
    def alpha(self, p):
        val = np.ones_like(p[..., 0])
        return val 

    @cartesian
    def beta(self, p):
        flag = self.interface(p)<0

        om = self.om
        eps = self.eps
        sig0 = self.sig0
        sig1 = self.sig1

        val = np.zeros_like(p[..., 0], dtype=np.complex128)
        val[flag] = np.ones_like(p[..., 0])[flag]*(om**2*(eps + sig0/om * (1j)))
        val[~flag] = np.ones_like(p[..., 0])[~flag]*(om**2*(eps + sig1/om * (1j)))
        return val 

    @cartesian
    def dirichlet(self, p, t):
        # p : (NE, 2)
        val = np.zeros_like(p[..., 0])
        return val 

    def get_deepest_mesh(self):
        return self.mesher.mesh0

    def get_mesh(self, nx=2**10, ny=2**8):
        """
        level 要小于 init 里的 level
        """
        #mesh = interfacemesh2d([0, 4, 0, 1], self.interface, 2**5, 2**3)
        mesh = TriangleMesh.interfacemesh_generator([0, 4, 0, 1], nx, ny,
                self.interface)
        return mesh 

fname = sys.argv[1]
if fname == "two_circle":
    #pde = PDE3(level=1, om = 1000, eps=0.01, interfacetype='double_circle')
    pde = PDE2()
    mesh = pde.get_mesh(2**9, 2**7)
elif fname == "two_band":
    pde = PDE3(level=1, om = 1000, eps=0.01, interfacetype='double_band')
    mesh = pde.get_mesh(2**13, 2**11)
elif fname == "five_band":
    pde = PDE3(om = 1000, eps=0.01)
    pde = PDE3(level=1, om = 1000, eps=0.01, interfacetype='five_band')
    mesh = pde.get_mesh(2**12, 2**10)
elif fname == "poly":
    pde = PDE3(level=1, om = 1000, eps=0.01, interfacetype='poly')
    points = np.array([[311+2/3, 50], [290+0.625, 47+2/3], [225, 68+1/39],
     [206.25, 19+1/3], [198.05, 66+2/3], [178.125, 15+1/3], [175+1/3, 62+1/3],
     [100+1/3, 69-1/4], [243.75, 82+1/3], [287.5, 52+1/3]])/100
    points[:, 1] = 1-points[:, 1]
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.0001)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0005)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    model = gmsh.model.occ
    box = model.add_rectangle(0, 0, 0, 4, 1)
    p=[]
    l=[]
    for i in range(len(points)):
        p.append(model.addPoint(points[i, 0], points[i, 1], 0, 1.0))
    for i in range(len(p)):
        l.append(model.addLine(p[i], p[(i+1)%len(p)]))
    curve_loop = model.addCurveLoop(l)
    model.fragment([(1, ll) for ll in l], [(2, box)])
    model.synchronize()
    gmsh.model.mesh.generate(2)
    #gmsh.fltk.initialize()
    #gmsh.fltk.run()

    node = gmsh.model.mesh.get_nodes()[1].reshape(-1, 3)[:, :2]
    NN = node.shape[0]
    nid2tag = gmsh.model.mesh.get_nodes()[0]
    tag2nid = np.zeros(NN+1000, dtype = np.int_)
    tag2nid[nid2tag] = np.arange(NN)
    cell = gmsh.model.mesh.get_elements(2, -1)[2][0].reshape(-1, 3)
    cell = tag2nid[cell]
    mesh = TriangleMesh(node, cell)

space = FirstNedelecFiniteElementSpace2d(mesh)
#space = SecondNedelecFiniteElementSpace2d(mesh, 1)
gdof = space.dof.number_of_global_dofs()

bc = DirichletBC(space, pde.dirichlet) 

M = space.mass_matrix(c=pde.beta)
C = space.curl_matrix(c=pde.alpha)

@cartesian
def sr(x) : return pde.source(x).real
br = space.source_vector(sr)
@cartesian
def si(x) : return pde.source(x).imag
bi = space.source_vector(si)

b = br + bi*(1j)
B = C-M 
uh = space.function()
B, b = bc.apply(B, b)
uh[:] = Solve(B, b).real

bc = np.array([[1/3, 1/3, 1/3]])
val = space.value(uh, bc) #(NQ, NC, 2)
cval = space.curl_value(uh, bc) #(NQ, NC)

cm = mesh.entity_measure("cell")
cbary = mesh.entity_barycenter('cell')

mesh.celldata['cellval'] = np.average(val, axis=0)
mesh.celldata['cellcval'] = np.average(cval, axis=0)

NN = mesh.number_of_nodes() 
cell = mesh.entity('cell')

bc = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float_)
val = space.value(uh, bc) #(NQ, NC, 2)
cval = space.curl_value(uh, bc) #(NQ, NC)

NN2C = np.zeros(NN, dtype=np.int_)
nval = np.zeros((NN, 2), dtype=np.float_)
ncval = np.zeros(NN, dtype=np.float_)

np.add.at(NN2C, cell, 1)
np.add.at(nval, cell.T, val)
np.add.at(ncval, cell.T, cval)

nval = nval/NN2C[:, None]
ncval = ncval/NN2C

mesh.nodedata['val'] = nval
mesh.nodedata['cval'] = ncval

mesh.to_vtk(fname = fname+'_fem.vtu')

# 画图
if 0:
    puh = val[0]
    cuh = cval[0]
    print(puh.shape)
    fig = plt.figure()
    axes = fig.gca()
    mesh.add_plot(plt, cellcolor=puh[:, 0], cmap=plt.cm.jet, linewidths=0, 
            showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    plt.title('x component of the real part of $u_h^8$')
    plt.savefig(fname+'_fem_puhx.svg')

    fig = plt.figure()
    axes = fig.gca()
    mesh.add_plot(plt, cellcolor=puh[:, 1], cmap=plt.cm.jet, linewidths=0, 
            showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    plt.title('y component of the real part of $u_h^8$')
    plt.savefig(fname+'_fem_puhy.svg')

    fig = plt.figure()
    axes = fig.gca()
    mesh.add_plot(plt, cellcolor=np.sqrt(puh[:, 0]**2 + puh[:, 1]**2), cmap=plt.cm.jet, linewidths=0, 
            showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    plt.title('magnitude the real part of $u_h^8$')
    plt.savefig(fname+'_fem_puhxy.svg')

    fig = plt.figure()
    axes = fig.gca()
    mesh.add_plot(plt, cellcolor=cuh, cmap=plt.cm.jet, linewidths=0, 
            showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
    plt.title('rot of the real part of $u_h^8$')
    plt.savefig(fname+'_fem_cuh.svg')


