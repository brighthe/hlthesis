import numpy as np

from fealpy.decorator import cartesian

class CosCosCosData:
    def __init__(self):
        pass
    def domain(self):
        return np.array([0, 1, 0, 1, 0, 1])

    @cartesian
    def solution(self, p):
        """ the exact solution
        """
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        u = np.cos(np.pi*x)*np.cos(np.pi*y)*np.cos(np.pi*z)
        return u

    @cartesian
    def gradient(self, p):
        """ The gradient of the exact solution
        """
        pi = np.pi
        sin = np.sin
        cos = np.cos
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        val = np.zeros(p.shape, dtype=p.dtype)
        val[..., 0] = -pi*sin(pi*x)*cos(pi*y)*cos(pi*z)
        val[..., 1] = -pi*cos(pi*x)*sin(pi*y)*cos(pi*z)
        val[..., 2] = -pi*cos(pi*x)*cos(pi*y)*sin(pi*z)
        return val

    @cartesian
    def flux(self, p):
        return -self.gradient(p)

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        val = 3*np.pi**2*np.cos(np.pi*x)*np.cos(np.pi*y)*np.cos(np.pi*z)
        return val

    @cartesian
    def dirichlet(self, p):
        """Dilichlet boundary condition
        """
        return self.solution(p)
    
    @cartesian
    def neumann(self, p, n):
        """ 
        Neuman  boundary condition

        Parameters
        ----------

        p: (NQ, NE, 3)
        n: (NE, 3)

        grad*n : (NQ, NE, 3)
        """
        grad = self.gradient(p) # (NQ, NE, 3)
        val = np.sum(grad*n, axis=-1) # (NQ, NE)
        return val

    @cartesian
    def robin(self, p, n):
        grad = self.gradient(p) # (NQ, NE, 3)
        val = np.sum(grad*n, axis=-1)
        shape = len(val.shape)*(1, )
        kappa = np.array([1.0], dtype=np.float64).reshape(shape)
        val += self.solution(p) 
        return val, kappa

