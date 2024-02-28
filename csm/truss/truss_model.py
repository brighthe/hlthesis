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
