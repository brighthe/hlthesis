import numpy as np

from fealpy.decorator import cartesian

class Sin4PiData:
    """
    -u''(x) = 16 pi^2 sin(4 pi x), 
       u(0) = 0, u(1) = 0.
    exact solution：
       u(x) = sin(4 pi x).
    """
    def domain(self):
        """
        @brief:    Get the domain of the PDE model

        @return:   A list representing the domain of the PDE model
        """
        return [0, 1]

    @cartesian    
    def solution(self, p):
        """
        @brief:    Calculate the exact solution of the PDE model

        @param p:  An array of the independent variable x
        @return:   The exact solution of the PDE model at the given points
        """
        val = np.sin(4*np.pi*p)

        return val
    
    @cartesian    
    def source(self, p):
        """
        @brief:    Calculate the source term of the PDE model

        @param p:  An array of the independent variable x
        @return:   The source term of the PDE model at the given points
        """
        val = 16*np.pi**2*np.sin(4*np.pi*p)

        return val
    
    @cartesian    
    def gradient(self, p):
        """
        @brief:    Calculate the gradient of the exact solution of the PDE model

        @param p:  An array of the independent variable x

        @return:   The gradient of the exact solution of the PDE model at the given points
        """
        val = 4*np.pi*np.cos(4*np.pi*p)

        return val

    @cartesian    
    def dirichlet(self, p):
        """
        @brief: Dirichlet BC
        """
        return self.solution(p)

    @cartesian
    def is_dirichlet_boundary(self, p):
        """
        @brief Dirichlet 边界的判断函数
        """
        flag = (np.abs(p - 1.0) < 1e-12) | (np.abs( p -  0.0) < 1e-12)

        return flag
