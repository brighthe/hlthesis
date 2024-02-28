import numpy as np

from fealpy.mesh import QuadrangleMesh

# Cantilever
class TopRBFPlsm:

    def __init__(self, nelx: int = 60, nely: int = 30, volfrac: float = 0.5):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): 沿设计区域水平方向的单元数.
        - nely (int): 沿设计区域垂直方向的单元数.
        - volfrac (float) : 规定的体积分数.
        '''

        self._nelx = nelx
        self._nely = nely
        self._volfrac = volfrac

        #nx = self._nelx
        #ny = self._nely
        #x = np.linspace(0, nx, nx + 1)
        #y = np.linspace(ny, 0, ny + 1)
        #xv, yv = np.meshgrid(x, y, indexing='ij')
        #nodes = np.vstack([xv.ravel(), yv.ravel()]).T
        #cells = []
        #for j in range(nx):
        #    for i in range(ny):
        #        top_left = i + ny * j + j
        #        top_right = top_left + 1
        #        bottom_left = top_left + ny + 1
        #        bottom_right = bottom_left + 1
        #        cells.append([top_left, bottom_left, bottom_right, top_right])
        #node = nodes
        #cell = np.array(cells)
        #self._mesh = QuadrangleMesh(node=node, cell=cell)

        self._mesh = QuadrangleMesh.from_box(box = [0, self._nelx, 0, self._nely], \
                                              nx = self._nelx, ny = self._nely)

    def lsf_init(self, nelx, nely, mesh):
        node = mesh.entity('node') # 按列增加
        # 网格中节点的 x 坐标
        X = node[:, 0].reshape(nelx+1, nely+1).T
        print("X:", X.shape, "\n", X)
        # 网格中节点的 y 坐标
        Y = node[:, 1].reshape(nelx+1, nely+1).T
        print("Y:", Y.shape, "\n", Y)
        # 初始孔洞的半径
        r = nely * 0.1
        # hX 和 hY 分别是初始孔洞的中心处的 x 和 y 坐标
        hX = nelx * np.concatenate([np.tile([1/6, 5/6], 3), np.tile([0, 1/3, 2/3, 1], 2), [1/2]])
        print("hX:", hX.shape, "\n", hX)
        hY = nely * np.concatenate([np.repeat([0, 1/2, 1], 2), np.repeat([1/4, 3/4], 4), [1/2]])
        print("hY:", hY.shape, "\n", hY)

        # dX 和 dY 分别是所有网格点在 x 方向上和 y 方向上与初始孔洞的中心之间的距离差
        dX = X[:, :, np.newaxis] - hX
        dY = Y[:, :, np.newaxis] - hY

        # 计算网格点到最近孔洞附近的欧氏距离，并限制在 -3 到 3 之间
        Phi = np.sqrt(dX**2 + dY**2) - r
        Phi = np.min(Phi, axis=2)
        Phi = np.clip(Phi, -3, 3)

        return Phi


