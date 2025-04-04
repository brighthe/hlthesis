import numpy as np

from fealpy.decorator import cartesian

class SinSin4PiExpData:
    def __init__(self, D=[0, 1, 0, 1], T=[0, 1]):
        """
        @brief 模型初始化函数
        
        @param[in] D 模型空间定义域
        @param[in] T 模型时间定义域
        """
        self._domain = D 
        self._duration = T 

    def domain(self):
        """
        @brief 空间区间
        """
        return self._domain

    def duration(self):
        """
        @brief 时间区间
        """
        return self._duration 

    @cartesian
    def solution(self, p, t):
        """
        @brief 真解函数

        @param[in] p numpy.ndarray, 空间点
        @param[in] t float, 时间点 

        @return 真解函数值
        """
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return np.sin(4*pi*x)*np.sin(4*pi*y)*np.exp(-20*t) 

    @cartesian
    def init_solution(self, p):
        """
        @brief 初始解

        @param[in] p numpy.ndarray, 空间点
        @param[in] t float, 时间点 

        @return 初始解函数值
        """
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]

        return np.sin(4*pi*x)*np.sin(4*pi*y)
        
    @cartesian
    def source(self, p, t):
        """
        @brief 方程右端项 

        @param[in] p numpy.ndarray, 空间点
        @param[in] t float, 时间点 

        @return 方程右端函数值
        """
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return -20*np.exp(-20*t)*np.sin(4*pi*x)*np.sin(4*pi*x) + 32*pi**2*np.exp(-20*t)*np.sin(4*pi*x)*np.sin(4*pi*y)
    
    @cartesian
    def gradient(self, p, t):
        """
        @brief 真解导数 

        @param[in] p numpy.ndarray, 空间点
        @param[in] t float, 时间点 

        @return 真解导函数值
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = 4*pi*np.cos(4*pi*x)*np.sin(4*pi*y)*np.exp(-20*t)
        val[..., 1] = 4*pi*np.sin(4*pi*x)*np.cos(4*pi*y)*np.exp(-20*t)
        return val
    
    @cartesian    
    def dirichlet(self, p, t):
        """
        @brief Dirichlet 边界条件

        @param[in] p numpy.ndarray, 空间点
        @param[in] t float, 时间点 
        """
        return self.solution(p, t)

    @cartesian
    def is_dirichlet_boundary(self, p):
        """
        @brief Dirichlet 边界的判断函数
        """
        x = p[..., 0]
        y = p[..., 1]
        flag1 = (np.abs(x - 1.0) < 1e-12) | (np.abs(x -  0.0) < 1e-12)
        flag2 = (np.abs(y - 1.0) < 1e-12) | (np.abs(y -  0.0) < 1e-12)

        return flag1 | flag2

class SinSinExpPDEData: 
    def __init__(self, D=[0, 1, 0, 1], T=[0, 1]):
        self._domain = D 
        self._duration = T 

    def domain(self):
        return self._domain
    
    def duration(self):
        return self._duration 
    
    @cartesian
    def solution(self, p, t):
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return np.sin(pi*x)*np.sin(pi*y)*np.exp(-2*pi*t) 
    
    @cartesian
    def init_solution(self, p):
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return np.sin(pi*x)*np.sin(pi*y)
    
    @cartesian
    def source(self, p, t):
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return np.zeros(x.shape)
    
    @cartesian
    def gradient(self, p, t):
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = pi*np.cos(pi*x)*np.sin(pi*y)*np.exp(-2*pi*t)
        val[..., 1] = pi*np.sin(pi*x)*np.cos(pi*y)*np.exp(-2*pi*t)
        return val
    
    @cartesian    
    def dirichlet(self, p, t):
        return self.solution(p, t)
