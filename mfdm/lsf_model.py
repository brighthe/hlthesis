import numpy as np

from fealpy.decorator  import cartesian 

from fealpy.mesh.polygon_mesh import PolygonMesh

from fealpy.mesh import TriangleMesh

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

    def triangle_mesh(self, domain, nx, ny):
        mesh = TriangleMesh.from_box(box=domain, nx=nx, ny=ny)

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

    # Define the velocity field $u$ for the evolution
    @cartesian
    def velocity_field_t(self, p, t):
        x = p[..., 0]
        y = p[..., 1]
        u = np.zeros(p.shape)
        u[..., 0] = np.sin((np.pi*x))**2 * np.sin(2*np.pi*y) * np.cos(np.pi*t)
        u[..., 1] = -np.sin((np.pi*y))**2 * np.sin(2*np.pi*x) * np.cos(np.pi*t)
        return u

    # Initial level set function $\phi0$ representing the circle
    @cartesian
    def circle(self, p, index=None):
        x = p[...,0]
        y = p[...,1]
        val = np.sqrt((x-0.5)**2 + (y-0.75)**2) - 0.15

        return val

    @cartesian
    def grad_circle(self, p, index=None):
        x = p[..., 0]
        y = p[..., 1]
        val = np.zeros(p.shape, dtype=np.float64)
        denom = np.sqrt((x - 0.5)**2 + (y - 0.75)**2)
        denom[denom == 0] = 1e-12

        val[..., 0] = (x - 0.5) / denom
        val[..., 1] = (y - 0.75) / denom

        return val


    @cartesian
    def scalar_product(self, p, index=None):
        u = self.velocity_field(p)
        grad_phi = self.grad_circle(p)
        dot_product = np.sum(u * grad_phi, axis=-1)

        return dot_product

