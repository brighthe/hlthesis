import numpy as np

from fealpy.decorator  import cartesian 

from fealpy.mesh import TriangleMesh

class ClassicalLsfData():
    def __init__(self):
        pass

    def domain(self):
        return [0, 1, 0, 1]

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

# Initial level set function $\phi0$ representing the circle
    @cartesian
    def circle(self, p):
        x = p[...,0]
        y = p[...,1]
        val = np.sqrt((x-0.5)**2+(y-0.75)**2)-0.15
        return val
