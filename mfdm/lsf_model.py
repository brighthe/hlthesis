import numpy as np

from fealpy.decorator  import cartesian 

from fealpy.mesh.polygon_mesh import PolygonMesh

class ClassicalLsfData():
    def __init__(self):
        pass

    def domain(self):
        return [0, 1, 0, 1]

    def polygon_mesh(self):
        node = np.array([
                        [0.0, 0.0], [0.0, 0.5], [0.0, 1.0],
                        [0.5, 0.0], [0.5, 0.5], [0.5, 1.0],
                        [1.0, 0.0], [1.0, 0.5], [1.0, 1.0]
                        ], dtype=np.float64)

        cell = np.array([0, 3, 4, 1, 3, 6, 7, 4, 1, 4, 5, 2, 4, 7, 8, 4, 8, 5], dtype=np.int_)
        cellLocation = np.array([0, 4, 8, 12, 15, 18], dtype=np.int_)
        mesh = PolygonMesh(node=node, cell=cell, cellLocation=cellLocation)

        return mesh

    def polygon_mesh_2(self, n=2):
        mesh = PolygonMesh.from_unit_square(nx=n, ny=n)

        return mesh

# Define the velocity field $u$ for the evolution
    @cartesian
    def velocity_field(self, p):
        x = p[..., 0]
        y = p[..., 1]
        u = np.zeros(p.shape)
        u[..., 0] = np.sin((np.pi*x))**2 * np.sin(2*np.pi*y)
        u[..., 1] = -np.sin((np.pi*y))**2 * np.sin(2*np.pi*x)

        return u

# Initial level set function $\phi0$ representing the circle
    @cartesian
    def circle(self, p):
        x = p[...,0]
        y = p[...,1]
        val = np.sqrt((x-0.5)**2 + (y-0.75)**2) - 0.15

        return val
