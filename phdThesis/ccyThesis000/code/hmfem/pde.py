import numpy as np
import sympy as sp

from fealpy.mesh import TriangleMesh
from fealpy.decorator import cartesian

class LaplacePDE():
    def __init__(self, u):
        x = sp.symbols("x")
        y = sp.symbols("y")

        ux = sp.diff(u, x)
        uy = sp.diff(u, y)
        uxx = sp.diff(ux, x)
        uyy = sp.diff(uy, y)

        Lu = uxx+uyy
        self.u = sp.lambdify(('x', 'y'), u, 'numpy') 

        self.ux = sp.lambdify(('x', 'y'), ux, 'numpy') 
        self.uy = sp.lambdify(('x', 'y'), uy, 'numpy') 
        self.Lu = sp.lambdify(('x', 'y'), Lu, 'numpy')

    def init_mesh(self, n=1, meshtype='poly'):
        mesh = TriangleMesh.from_box([0, 1, 0, 1], 4, 4)
        return mesh

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        return -self.Lu(x, y)

    def solution(self, p):
        x = p[..., 0]
        y = p[..., 1]
        return self.u(x, y) 

    def gradient(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape, dtype=np.float_)
        val[..., 0] = self.ux(x, y)
        val[..., 1] = self.uy(x, y) 
        return val

    def dirichlet(self, p):
        return self.solution(p)

class DoubleLaplacePDE():
    def __init__(self, u):
        x = sp.symbols("x")
        y = sp.symbols("y")

        ux = sp.diff(u, x)
        uy = sp.diff(u, y)
        uxx = sp.diff(ux, x)
        uyy = sp.diff(uy, y)
        uxy = sp.diff(ux, y)

        Lu = uxx+uyy
        L2u = sp.diff(Lu, x, 2) + sp.diff(Lu, y, 2)

        self.u = sp.lambdify(('x', 'y'), u, 'numpy') 

        self.ux = sp.lambdify(('x', 'y'), ux, 'numpy') 
        self.uy = sp.lambdify(('x', 'y'), uy, 'numpy') 

        self.uxx = sp.lambdify(('x', 'y'), uxx, 'numpy') 
        self.uyy = sp.lambdify(('x', 'y'), uyy, 'numpy') 
        self.uxy = sp.lambdify(('x', 'y'), uxy, 'numpy') 

        self.L2u = sp.lambdify(('x', 'y'), L2u, 'numpy')

    def init_mesh(self, n=1, meshtype='poly'):
        mesh = TriangleMesh.from_box([0, 1, 0, 1], 4, 4)
        return mesh

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        return self.L2u(x, y)

    def solution(self, p):
        x = p[..., 0]
        y = p[..., 1]
        return self.u(x, y) 

    def gradient(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape, dtype=np.float_)
        val[..., 0] = self.ux(x, y)
        val[..., 1] = self.uy(x, y) 
        return val

    @cartesian
    def hessian(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape[:-1]+(3, ), dtype=np.float_)
        val[..., 0] = self.uxx(x, y) 
        val[..., 1] = self.uxy(x, y) 
        val[..., 2] = self.uyy(x, y) 
        return val

    def dirichlet(self, p):
        return self.solution(p)

class TripleLaplacePDE():
    def __init__(self, u):
        x = sp.symbols("x")
        y = sp.symbols("y")

        ux = sp.diff(u, x)
        uy = sp.diff(u, y)
        uxx = sp.diff(ux, x)
        uyy = sp.diff(uy, y)
        uxy = sp.diff(ux, y)

        uxxx = sp.diff(uxx, x)
        uxxy = sp.diff(uxx, y)
        uxyx = sp.diff(uxy, x)
        uxyy = sp.diff(uxy, y)
        uyxx = sp.diff(uxy, x)
        uyxy = sp.diff(uxy, y)
        uyyx = sp.diff(uyy, x)
        uyyy = sp.diff(uyy, y)

        Lu = uxx+uyy
        L2u = sp.diff(Lu, x, 2) + sp.diff(Lu, y, 2)
        L3u = sp.diff(L2u, x, 2) + sp.diff(L2u, y, 2)

        self.u = sp.lambdify(('x', 'y'), u, 'numpy') 

        self.ux = sp.lambdify(('x', 'y'), ux, 'numpy') 
        self.uy = sp.lambdify(('x', 'y'), uy, 'numpy') 

        self.uxx = sp.lambdify(('x', 'y'), uxx, 'numpy') 
        self.uyy = sp.lambdify(('x', 'y'), uyy, 'numpy') 
        self.uxy = sp.lambdify(('x', 'y'), uxy, 'numpy') 

        self.uxxx = sp.lambdify(('x', 'y'), uxxx, 'numpy') 
        self.uxxy = sp.lambdify(('x', 'y'), uxxy, 'numpy') 
        self.uxyx = sp.lambdify(('x', 'y'), uxyx, 'numpy') 
        self.uxyy = sp.lambdify(('x', 'y'), uxyy, 'numpy') 
        self.uyxx = sp.lambdify(('x', 'y'), uyxx, 'numpy') 
        self.uyxy = sp.lambdify(('x', 'y'), uyxy, 'numpy') 
        self.uyyx = sp.lambdify(('x', 'y'), uyyx, 'numpy') 
        self.uyyy = sp.lambdify(('x', 'y'), uyyy, 'numpy') 

        self.L3u = sp.lambdify(('x', 'y'), L3u, 'numpy')

    def init_mesh(self, n=1, meshtype='poly'):
        mesh = TriangleMesh.from_box([0, 1, 0, 1], 4, 4)
        return mesh

    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        return -self.L3u(x, y) 

    def solution(self, p):
        x = p[..., 0]
        y = p[..., 1]
        return self.u(x, y) 

    def gradient(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape, dtype=np.float_)
        val[..., 0] = self.ux(x, y)
        val[..., 1] = self.uy(x, y) 
        return val

    def hessian(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape[:-1]+(3, ), dtype=np.float_)
        val[..., 0] = self.uxx(x, y) 
        val[..., 1] = self.uxy(x, y) 
        val[..., 2] = self.uyy(x, y) 
        return val

    def grad_3(self, p):
        x = p[..., 0]
        y = p[..., 1]
        sin = np.sin
        pi = np.pi
        cos = np.cos
        val = np.zeros(p.shape[:-1]+(4, ), dtype=np.float_)
        val[..., 0] = self.uxxx(x, y) 
        val[..., 1] = self.uyxx(x, y)
        val[..., 2] = self.uyyx(x, y)
        val[..., 3] = self.uyyy(x, y)
        return val

    def dirichlet(self, p):
        return self.solution(p)
