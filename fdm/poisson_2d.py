import numpy as np

from fealpy.decorator import cartesian

class SinSinData:
    def domain(self):
        """
        @brief 得到 PDE 模型的区域
        @return: 表示 PDE 模型的区域的列表
        """
        return np.array([0, 1, 0, 1])

    @cartesian
    def solution(self, p):
        """
        @brief 计算 PDE 模型的精确解
        @param p: 自标量 x,y 的数组
        @return: PDE 模型在给定点的精确解
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = np.sin(pi*x)*np.sin(pi*y)
        return val 
    
    @cartesian
    def source(self, p):
        """
        @brief: 计算 PDE 模型的原项 
        @param p: 自标量 x,y 的数组
        @return: PDE 模型在给定点处的源项
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        val = 2*pi*pi*np.sin(pi*x)*np.sin(pi*y)
        return val
    
    @cartesian    
    def dirichlet(self, p):
        return self.solution(p)
        
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

