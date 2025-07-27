#!/usr/bin/env python3
# 

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.vector import CoordSys3D, Del, curl

from fealpy.geometry import CircleCurve, FoldCurve
from interface_mesh_2d import HalfEdgeMesh2dWithInterface
from fealpy.mesh.interface_mesh_generator import interfacemesh2d 

from fealpy.pde.MaxwellPDE2d import SinData as PDE
from fealpy.functionspace import FirstNedelecFiniteElementSpace2d
from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh, MeshFactory
from fealpy.decorator import cartesian, barycentric
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包
from fealpy.pde.MaxwellPDE2d import MaxwellPDE2d


from scipy.sparse.linalg import spsolve, cg

class PDE1(MaxwellPDE2d):
    def __init__(self):
        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        super(PDE1, self).__init__(f)

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*2
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
        box = [-1, 1, -1, 1]
        mesh = MeshFactory.boxmesh2d(box, nx=nx, ny=ny, meshtype='tri')
        return mesh 

class PDE2():
    """
    @brief 0: 内部; 1: 外部
    """
    def __init__(self, k1=20, r0=np.pi/5, r1=1, mu0=0.1, mu1=0.1):
        self.interface = CircleCurve([0.0, 0.0], r0)

        C = CoordSys3D('C')
        k0 = k1*(r1**2-r0**2)
        u0 = -mu0*k0*(r0**2 - C.x**2 - C.y**2)*C.y*C.i - mu0*k0*(r0**2 - 
                C.x**2 - C.y**2)*C.x*C.j
        u1 = -mu1*k1*(r1**2 - C.x**2 - C.y**2)*(r0**2 - 
                C.x**2 - C.y**2)*C.y*C.i -mu1*k1*(r1**2 - 
                        C.x**2 - C.y**2)*(r0**2 - C.x**2 - C.y**2)*C.x*C.j 

        self.pde0 = MaxwellPDE2d(u0)
        self.pde1 = MaxwellPDE2d(u1)

    def solution(self, p, *args):
        flag = self.interface(p)<0
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag] = self.pde0.solution(p[flag])
        val[~flag] = self.pde1.solution(p[~flag])
        return val

    def rotrotsolution(self, p, *args):
        flag = self.interface(p)<0
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag] = self.pde0.curl_curl_solution(p[flag])
        val[~flag] = self.pde1.curl_curl_solution(p[~flag])
        return val

    @cartesian
    def alpha(self, p):
        flag = self.interface(p)<0
        val = 1*np.ones(p.shape[:-1], dtype=p.dtype)
        val[flag] = 1
        return val 

    @cartesian
    def beta(self, p):
        flag = self.interface(p)<0
        val = 10*np.ones(p.shape[:-1], dtype=p.dtype)
        val[flag] = 1
        return val 

    @cartesian
    def curl_solution(self, p):
        flag = self.interface(p)<0
        val = np.zeros(p.shape[:-1], dtype=p.dtype)
        val[flag] = self.pde0.curl_solution(p[flag])
        val[~flag] = self.pde1.curl_solution(p[~flag])
        return val

    @cartesian
    def source(self, p):
        alpha = self.alpha(p)[..., None]
        beta = self.beta(p)[..., None]
        val = alpha*self.rotrotsolution(p)-beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx, ny):
        mesh = interfacemesh2d([-1.1, 1, -1.1, 1], self.interface, nx,
                meshtype='tri')
        node = mesh.entity('node')
        cell = mesh.entity('cell')[0].reshape(-1, 3)
        return TriangleMesh(node, cell)

pde = PDE2()
maxit = 7
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float)
NDof = np.zeros(maxit, dtype=np.int_)
for i in range(maxit):
    print("The {}-th computation:".format(i))
    mesh = pde.get_mesh(nx = 2**(i+2), ny = 2**(i+2))

    space = FirstNedelecFiniteElementSpace2d(mesh)
    gdof = space.dof.number_of_global_dofs()
    NDof[i] = gdof

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix(c=pde.beta)
    C = space.curl_matrix(c=pde.alpha)
    b = space.source_vector(pde.source)
    B = C-M 

    uh = space.function()
    B, b = bc.apply(B, b)
    uh[:] = spsolve(B, b)
    # 计算误差

    errorMatrix[0, i] = space.integralalg.L2_error(pde.solution, uh)
    #errorMatrix[1, i] = space.integralalg.error(pde.curl_solution, uh)
    print(errorMatrix)

showmultirate(plt, 2, NDof, errorMatrix,  errorType, propsize=20)
show_error_table(NDof, errorType, errorMatrix)
plt.show()

