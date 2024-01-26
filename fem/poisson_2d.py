import numpy as np

from fealpy.decorator import cartesian

class CosCosData:
    """
        -\\Delta u = f
        u = cos(pi*x)*cos(pi*y)
    """
    def __init__(self, kappa=1.0):
        self.kappa = kappa # Robin 条件中的系数

    def domain(self):
        """
        @brief 模型定义域
        """
        return np.array([0, 1, 0, 1])

    @cartesian
    def solution(self, p):
        """  
        @brief 模型真解
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.cos(pi*x)*np.cos(pi*y)

        return val


    @cartesian
    def source(self, p):
        """ 
        @brief 源项
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = 2*pi*pi*np.cos(pi*x)*np.cos(pi*y)

        return val

    @cartesian
    def gradient(self, p):
        """  
        @brief 真解梯度
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = -pi*np.sin(pi*x)*np.cos(pi*y)
        val[..., 1] = -pi*np.cos(pi*x)*np.sin(pi*y)
        return val # val.shape == p.shape

    @cartesian
    def flux(self, p):
        """
        @brief 真解通量
        """
        return -self.gradient(p)

    @cartesian
    def dirichlet(self, p):
        """
        @brief Dirichlet 边界条件 
        """
        return self.solution(p)

    @cartesian
    def is_dirichlet_boundary(self, p):
        """
        @brief Dirichlet 边界的判断函数
        """
        x = p[..., 0]
        y = p[..., 1]
        flag1 = (np.abs(x - 1.0) < 1e-12) | (np.abs( x -  0.0) < 1e-12)
        flag2 = (np.abs(y - 1.0) < 1e-12) | (np.abs( y -  0.0) < 1e-12)

        return flag1 | flag2

    @cartesian
    def neumann(self, p, n):
        """ 
        @brief Neumann 边界条件
        """
        grad = self.gradient(p) # (NQ, NE, 2)
        val = np.sum(grad*n, axis=-1) # (NQ, NE)
        return val

    @cartesian
    def is_neumann_boundary(self, p):
        """
        @brief Neumann 边界的判断函数
        """
        x = p[..., 0]
        return np.abs(x - 1.0) < 1e-12

    @cartesian
    def robin(self, p, n):
        """
        @brief Robin 边界条件
        """
        grad = self.gradient(p) # (NQ, NE, 2)
        val = np.sum(grad*n, axis=-1)
        val += self.kappa*self.solution(p) 
        return val

    @cartesian
    def is_robin_boundary(self, p):
        """
        @brief Robin 边界条件判断函数
        """
        x = p[..., 0]
        return np.abs(x - 0.0) < 1e-12
