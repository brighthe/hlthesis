import numpy as np

from fealpy.decorator  import cartesian 
from fealpy.mesh import TriangleMesh
from fealpy.geometry.domain_2d import RectangleDomain

class BoxDomainData2d():
    """
    @brief 混合边界条件的线弹性问题模型
    @note 本模型假设在二维方形区域 [0,1] x [0,1] 内的线性弹性问题
    """
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
        return [0, 1, 0, 1]

    def init_mesh(self, n = 1):
        """
        @brief 初始化网格
        @param[in] n 网格加密的次数，默认值为 1
        @return 返回一个初始化后的网格对象
        """
        h = 0.5
        domain = RectangleDomain()
        mesh = TriangleMesh.from_domain_distmesh(domain, h, output=False)
        mesh.uniform_refine(n)

        return mesh 

    def triangle_mesh(self):
        mesh = TriangleMesh.from_box(box=[0, 1, 0, 1], nx=2, ny=2)

        return mesh

    def delaunay_mesh(self):
        import scipy.io as io
        data_1 = io.loadmat('p.mat')
        data_2 = io.loadmat('t.mat')

        node = np.transpose(data_1['p'])
        cell = np.transpose(data_2['t'][:3, ]) - 1
        cell = cell.astype(np.uint32)
        mesh = TriangleMesh(node, cell)

        return mesh

    @cartesian
    def source(self, p):
        """
        @brief 返回给定点的源项值 f
        @param[in] p 一个表示空间点坐标的数组
        @return 返回源项值
        """
        x = p[..., 0]
        y = p[..., 1]
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = 35/13*y - 35/13*y**2 + 10/13*x - 10/13*x**2
        val[..., 1] = -25/26*(-1+2*y) * (-1+2*x)

        return val

    @cartesian
    def solution(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = np.zeros(p.shape, dtype=np.float64)
        val[..., 0] = x*(1-x)*y*(1-y)
        val[..., 1] = 0

        return val


    @cartesian
    def dirichlet(self, p):
        """
        @brief 返回 Dirichlet 边界上的给定点的位移
        @param[in] p 一个表示空间点坐标的数组
        @return 返回位移值，这里返回常数向量 [0.0, 0.0]
        """
        val = np.zeros((p.shape[0], 2), dtype=np.float64)
        # val = np.array([0.0, 0.0], dtype=np.float64)

        return val

    @cartesian
    def is_dirichlet_boundary(self, p):
        """
        @brief 判断给定点是否在 Dirichlet 边界上
        @param[in] p 一个表示空间点坐标的数组
        @return 如果在 Dirichlet 边界上，返回 True，否则返回 False
        """
        x = p[..., 0]
        y = p[..., 1]
        flag1 = np.abs(x) < 1e-13
        flag2 = np.abs(x - 1) < 1e-13
        flagx = np.logical_or(flag1, flag2)
        flag3 = np.abs(y) < 1e-13
        flag4 = np.abs(y - 1) < 1e-13
        flagy = np.logical_or(flag3, flag4)
        flag = np.logical_or(flagx, flagy)

        return flag

