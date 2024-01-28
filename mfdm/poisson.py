import numpy as np

from fealpy.decorator import cartesian

from fealpy.mesh.polygon_mesh import PolygonMesh

class SinSinData:
    """
        -\\Delta u = f
        u = sin(pi*x)*sin(pi*y)
    """
    def __init__(self):
        pass

    def domain(self):
        return np.array([0, 1, 0, 1])

    def polygon_mesh(self):
        node = np.array([[0.0, 0.0], [0.0, 0.5], [0.0, 1.0],
                    [0.5, 0.0], [0.5, 0.5], [0.5, 1.0],
                    [1.0, 0.0], [1.0, 0.5], [1.0, 1.0]], dtype=np.float64)

        cell = np.array([0, 3, 4, 1, 3, 6, 7, 4, 1, 4, 5, 2, 4, 7, 8, 4, 8, 5], dtype=np.int_)
        cellLocation = np.array([0, 4, 8, 12, 15, 18], dtype=np.int_)
        mesh = PolygonMesh(node=node, cell=cell, cellLocation=cellLocation)

        return mesh

    def polygon_mesh_2(self, n=2):
        mesh = PolygonMesh.from_unit_square(nx=n, ny=n)

        return mesh

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
    def gradient_u(self, p):
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.zeros_like(p)
        val[..., 0] = pi*np.cos(pi*x)*np.sin(pi*y)
        val[..., 1] = pi*np.sin(pi*x)*np.cos(pi*y)
        return val # val.shape == p.shape

    @cartesian
    def div_u(self, p, index=None):
        x = p[..., 0]
        y = p[..., 1]
        value0 = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
        value1 = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
        value = value0+value1
        return value

    @cartesian
    def Dirichlet(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = 0
        return val

class SinSin5Data:
    """
        -\\Delta u = f
        u = sin(pi*x)*sin(pi*y)
    """
    def __init__(self):
        pass

    def domain(self):
        return np.array([0, 1, 0, 1])

    def polygon_mesh(self):
        node = np.array([[0.0, 0.0], [0.0, 0.5], [0.0, 1.0],
                    [0.5, 0.0], [0.5, 0.5], [0.5, 1.0],
                    [1.0, 0.0], [1.0, 0.5], [1.0, 1.0]], dtype=np.float64)

        cell = np.array([0, 3, 4, 1, 3, 6, 7, 4, 1, 4, 5, 2, 4, 7, 8, 4, 8, 5], dtype=np.int_)
        cellLocation = np.array([0, 4, 8, 12, 15, 18], dtype=np.int_)
        mesh = PolygonMesh(node=node, cell=cell, cellLocation=cellLocation)

        return mesh

    def polygon_mesh_2(self, n=2):
        mesh = PolygonMesh.from_unit_square(nx=n, ny=n)

        return mesh

    @cartesian
    def source(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = 5 * np.pi**2 *np.sin(2*np.pi*x) * np.sin(np.pi*y)
        return val

    @cartesian
    def Dirichlet(self, p):
        x = p[...,0]
        y = p[...,1]
        val = np.sin(2*np.pi*x) * np.sin(np.pi*y)
        #val = 0
        return val

    @cartesian
    def solution(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = np.sin(2*np.pi*x) * np.sin(np.pi*y)

        return val

    #@cartesian
    #def gradient_u(self, p, index=None):
    #    x = p[...,0]
    #    y = p[...,1]
    #    value = np.zeros_like(p)
    #    value[...,0] = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
    #    value[...,1] = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
    #    return value

    #@cartesian
    #def div_u(self, p, index=None):
    #    x = p[...,0]
    #    y = p[...,1]
    #    value0 = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
    #    value1 = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
    #    value = value0+value1
    #    return value

class CosCos5Data:
    """
        -\\Delta u = f
        u = sin(pi*x)*sin(pi*y)
    """
    def __init__(self):
        pass

    def domain(self):
        return np.array([0, 1, 0, 1])

    def polygon_mesh(self):
        node = np.array([[0.0, 0.0], [0.0, 0.5], [0.0, 1.0],
                    [0.5, 0.0], [0.5, 0.5], [0.5, 1.0],
                    [1.0, 0.0], [1.0, 0.5], [1.0, 1.0]], dtype=np.float64)

        cell = np.array([0, 3, 4, 1, 3, 6, 7, 4, 1, 4, 5, 2, 4, 7, 8, 4, 8, 5], dtype=np.int_)
        cellLocation = np.array([0, 4, 8, 12, 15, 18], dtype=np.int_)
        mesh = PolygonMesh(node=node, cell=cell, cellLocation=cellLocation)

        return mesh

    def polygon_mesh_2(self, n=2):
        mesh = PolygonMesh.from_unit_square(nx=n, ny=n)

        return mesh

    @cartesian
    def source(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = 5 * np.pi**2 *np.cos(2*np.pi*x) * np.cos(np.pi*y)
        return val

    @cartesian
    def Dirichlet(self, p):
        x = p[...,0]
        y = p[...,1]
        val = np.cos(2*np.pi*x) * np.cos(np.pi*y)
        return val

    @cartesian
    def solution(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = np.cos(2*np.pi*x) * np.cos(np.pi*y)
        return val

    #@cartesian
    #def gradient_u(self, p, index=None):
    #    x = p[...,0]
    #    y = p[...,1]
    #    value = np.zeros_like(p)
    #    value[...,0] = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
    #    value[...,1] = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
    #    return value

    #@cartesian
    #def div_u(self, p, index=None):
    #    x = p[...,0]
    #    y = p[...,1]
    #    value0 = np.pi*np.cos(np.pi*x)*np.sin(np.pi*y)
    #    value1 = np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)
    #    value = value0+value1
    #    return value

