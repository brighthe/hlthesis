import numpy as np

from fealpy.decorator import cartesian 
from fealpy.mesh import EdgeMesh


class Truss_3d():
    def __init__(self):
        """
        初始化函数

        在此函数中设置桁架模型的基本参数
        """
        self.A: float = 2000 # Cross-sectional area - mm^2
        self.E: float = 1500 # Elastic Modulus newton - ton/mm^2

    def init_mesh(self):
        """
        初始化网格结构

        此函数用于初始化桁架的网格结构

        返回:
        mesh: EdgeMesh, 初始化的网格对象
        """
        mesh = EdgeMesh.from_tower()

        return mesh

class Truss_2d():
    def __init__(self):
        self.A = 100 # Cross-sectional - area mm^2
        self.E = 29.5e4 # Elastic Modulus - newton/mm^2

    def init_mesh(self):
        mesh = EdgeMesh.from_four_bar_mesh()

        return mesh


    #@cartesian
    #def source(self, p):
    #    shape = len(p.shape[:-1])*(1,) + (-1, )
    #    val = np.zeros(shape, dtype=np.float_)

    #    return val 

    #@cartesian
    #def force(self):
    #    '''
    #    施加 (0, 900, 0) 的力，即平行于 y 轴方向大小为 900N 的力
    #    '''
    #    val = np.array([0, 900, 0])
    #    return val

    #def is_force_boundary(self, p):
    #    '''
    #    对第 0，1 号节点施加力
    #    '''
    #    return np.abs(p[..., 2]) == 5080

    #@cartesian
    #def dirichlet(self, p):
    #    shape = len(p.shape)*(1, )
    #    val = np.array([0.0])
    #    return val.reshape(shape)

    #@cartesian
    #def is_dirichlet_boundary(self, p):
    #    return np.abs(p[..., 2]) < 1e-12
