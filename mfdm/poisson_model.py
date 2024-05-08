import numpy as np

from fealpy.decorator import cartesian

from fealpy.mesh.polygon_mesh import PolygonMesh

class CosCosData:
    """
    \\-Delta u = 5*pi**2 * uexact
        uexact = cos(2*pi*x) * cos(pi*y)
    """
    @cartesian
    def source(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = 5*np.pi*np.pi*np.cos(2*np.pi*x)*np.cos(np.pi*y)

        return val

    @cartesian
    def solution(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = np.cos(2*np.pi*x)*np.cos(np.pi*y)

        return val

    @cartesian
    def Dirichlet(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = np.cos(2*np.pi*x)*np.cos(np.pi*y)
        return val

class ExpSinData:
    """
    \\-Delta u = 0
        uexact = exp(x)*sin(y)
    """
    @cartesian
    def solution(self, p, index=None):
        x = p[..., 0]
        y = p[..., 1]
        val = np.exp(x) * np.sin(y)
        return val

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = np.zeros_like(x)
        return val

    @cartesian
    def Dirichlet(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = np.exp(x) * np.sin(y)
        return val

class SinSinData:
    """
    -\\Delta u = 2*pi**2 * sin(pi*x)*sin(pi*y)
        uexact = sin(pi*x)*sin(pi*y)
    """
    def __init__(self):
        pass

    @cartesian
    def solution(self, p, index=None):
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.sin(pi*x) * np.sin(pi*y)

        return val

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = 2*pi*pi * np.sin(pi*x) * np.sin(pi*y)
        return val

    @cartesian
    def Dirichlet(self, p):
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.sin(pi*x) * np.sin(pi*y)
        return val
    #@cartesian
    #def gradient_u(self, p):
    #    x = p[..., 0]
    #    y = p[..., 1]
    #    pi = np.pi
    #    val = np.zeros_like(p)
    #    val[..., 0] = pi*np.cos(pi*x)*np.sin(pi*y)
    #    val[..., 1] = pi*np.sin(pi*x)*np.cos(pi*y)
    #    return val # val.shape == p.shape

    #@cartesian
    #def div_u(self, p, index=None):
    #    x = p[..., 0]
    #    y = p[..., 1]
    #    value0 = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
    #    value1 = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
    #    value = value0+value1
    #    return value

