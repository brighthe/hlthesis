import numpy as np

from fealpy.decorator  import cartesian 
from fealpy.mesh import TetrahedronMesh

class BoxDomainData3d():
    def __init__(self, E=1.0, nu=0.3):
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

    def init_mesh(self, nx, ny, nz):
        mesh = TetrahedronMesh.from_box(box=[0, 1, 0, 1, 0, 1], nx=nx, ny=ny, nz=nz)
        #mesh.uniform_refine(n)
        return mesh

    @cartesian
    def solution(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = x*(1-x)*y*(1-y)*z*(1-z)
        val[..., 1] = y*(1-y)*z*(1-z)*x*(1-x)
        val[..., 2] = z*(1-z)*x*(1-x)*y*(1-y)
        return val
    
    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        val = np.zeros(p.shape, dtype=np.float64)
        
        common_term = (1 - 2*x)*(1 - 2*y)*(1 - 2*z)
        factor = -(2*self.mu + self.lam)
        
        val[..., 0] = factor * common_term
        val[..., 1] = factor * common_term
        val[..., 2] = factor * common_term
        
        return val

    @cartesian
    def dirichlet(self, p):
        """
        @brief Dirichlet 边界条件
        """
        return self.solution(p)

    @cartesian
    def is_dirichlet_boundary(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        flag1 = np.abs(x) < 1e-12
        flag2 = np.abs(x - 1) < 1e-12
        flag3 = np.abs(y) < 1e-12
        flag4 = np.abs(y - 1) < 1e-12
        flag5 = np.abs(z) < 1e-12
        flag6 = np.abs(z - 1) < 1e-12
        return np.logical_or(np.logical_or(flag1, flag2),\
                             np.logical_or(np.logical_or(flag3, flag4), np.logical_or(flag5, flag6)))