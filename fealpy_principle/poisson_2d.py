import numpy as np

from fealpy.decorator import cartesian

class CosCosData:
    """
        -\\Delta u = f
        u = cos(pi*x)*cos(pi*y)
    """
    def __init__(self):
        pass

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


class CosCosData_test:
    """
        用于检测混合边界条件
        -\\Delta u = f
        u = cos(pi*x)*cos(pi*y)
    """
    def __init__(self):
        pass

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
        flag1 = (np.abs(x - 1.0) < 1e-12)
        flag2 = (np.abs( y -  0.0) < 1e-12)

        return flag1 | flag2