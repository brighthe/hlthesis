import numpy as np

from fealpy.decorator  import cartesian 
from fealpy.mesh import TetrahedronMesh

class BoxDomainData3d():
    def __init__(self, E = 1.0, nu = 0.3):
        """
        @brief 构造函数
        @param[in] E 弹性模量，默认值为 1.0
        @param[in] nu 泊松比，默认值为 0.3
        """
        self.E = E 
        self.nu = nu

        self.lam = self.nu*self.E/((1+self.nu)*(1-2*self.nu))
        self.mu = self.E/(2*(1+self.nu))

    def domain(self):
        return [0, 1, 0, 1, 0, 1]

    def triangle_mesh(self):
        mesh = TetrahedronMesh.from_box(box=self.domain, nx=2, ny=2, nz=2)
        return mesh


    @cartesian
    def source(self, p):
#        shape = len(p.shape[:-1])*(1,) + (-1, )
        val = self.d*self.g*self.rho
        return val 

    @cartesian
    def dirichlet(self, p):
        val = np.array([0.0, 0.0, 0.0])
        return val

    @cartesian
    def is_dirichlet_boundary(self, p):
        return np.abs(p[..., 0]) < 1e-12
    
    @cartesian
    def neumann(self, p, n):
        val = np.array([0.0, -50, 0.0], dtype=np.float64)
        return val

    @cartesian
    def is_neumann_boundary(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        flag = np.abs(y - 0.2) < 1e-13
        return flag
