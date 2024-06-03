import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopSimp:
    def __init__(self, nelx, nely, volfrac, penal, rmin):
        '''
        初始化拓扑优化问题

        Parameters:
        - nelx (int): 水平方向的单元数.
        - nely (int): 垂直方向的单元数.
        - volfrac (float): Volume fraction, 表示材料将占据的设计空间的期望分数.
        - penal (float): Penalization power, 在 SIMP 中控制中间密度的 penalization.
        - rmin (float): Filter 半径.
        '''
        self._nelx = nelx
        self._nely = nely
        self._volfrac = volfrac
        self._penal = penal
        self._rmin = rmin
        node = np.array([[0, 2], [0, 1], [0, 0],
                         [1, 2], [1, 1], [1, 0],
                         [2, 2], [2, 1], [2, 0]], dtype=np.float64)
        cell = np.array([[0, 3, 4, 1],
                         [1, 4, 5, 2],
                         [3, 6, 7, 4],
                         [4, 7, 8, 5]], dtype=np.int_)
        self._mesh = QuadrangleMesh(node=node, cell=cell)

        nx = self._nelx
        ny = self._nely
        x = np.linspace(0, nx, nx + 1)
        y = np.linspace(ny, 0, ny + 1)
        xv, yv = np.meshgrid(x, y, indexing='ij')
        nodes = np.vstack([xv.ravel(), yv.ravel()]).T
        cells = []
        for j in range(nx):
            for i in range(ny):
                top_left = i + ny * j + j
                top_right = top_left + 1
                bottom_left = top_left + ny + 1
                bottom_right = bottom_left + 1
                cells.append([top_left, bottom_left, bottom_right, top_right])
        node = nodes
        cell = np.array(cells)
        self._mesh = QuadrangleMesh(node=node, cell=cell)

    def FE(self, mesh, x, penal, KE, F, fixeddofs):
        """
        有限元计算位移.

        Parameters:
        - mesh:
        - x (ndarray - (nely, nelx) ): 设计变量.
        - penal (float): penalization power.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        - F ( ndarray - (gdof*GD, nLoads) ): 节点荷载.
        - fixeddofs (ndarray): 位移约束(supports).

        Returns:
        - uh ( ndarray - (gdof, GD, nLoads) ): 总位移.
        - ue ( ndarray - (NC, ldof*GD, nLoads) ): 单元位移.
        """
        from simp_beam_operator_integrator import BeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import spdiags

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        nLoads = F.shape[-1]
        uh = np.repeat(uh[:, :, np.newaxis], nLoads, axis=2)
        vspace = GD*(space, )
        ldof = vspace[0].number_of_local_dofs()
        vldof = ldof * GD

        integrator = BeamOperatorIntegrator(x=x, penal=penal, KE=KE)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        bform.assembly()
        K = bform.get_matrix()

        dflag = fixeddofs
        F = F - K@uh.reshape(-1, nLoads)
        bdIdx = np.zeros(K.shape[0], dtype=np.int_)
        bdIdx[dflag.flat] = 1
        D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
        D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
        K = D0@K@D0 + D1
        F[dflag.flat] = uh.reshape(-1, nLoads)[dflag.flat]

        # 线性方程组求解
        uh.flat[:] = spsolve(K, F)

        cell2dof = vspace[0].cell_to_dof()
        NC = mesh.number_of_cells()
        ue = np.zeros((NC, vldof, nLoads))

        for i in range(nLoads):
            reshaped_uh = uh[:, : ,i].reshape(-1)
            # 每个单元的自由度（每个节点两个自由度）
            updated_cell2dof = np.repeat(cell2dof*GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
            #print("updated_cell2dof:", updated_cell2dof.shape, "\n", updated_cell2dof)
            idx = np.array([0, 1, 4, 5, 6, 7, 2, 3], dtype=np.int_)
            # 用 Top 中的自由度替换 FEALPy 中的自由度
            updated_cell2dof = updated_cell2dof[:, idx]
            ue[:, :, i] = reshaped_uh[updated_cell2dof]

        return uh, ue

    def check(self, rmin, x, dc):
        """
        应用 mesh-independency filter 进行每个单元的灵敏度过滤.

        Parameters:
        - rmin (float): Filter 半径.
        - x (ndarray - (nely, nelx) ): 密度分布矩阵.
        - dc (ndarray - (nely, nelx) ): 原始的灵敏度矩阵.

        Returns:
        - dcn(ndarray - (nely, nelx) ): 过滤后的灵敏度矩阵.
        """

        nely, nelx = x.shape
        dcn = np.zeros((nely, nelx))

        # 计算过滤器半径
        r = int(rmin)

        for i in range(nelx):
            for j in range(nely):
                sum_val = 0.0
                # 确定邻域的范围
                min_x = max(i - r, 0)
                max_x = min(i + r + 1, nelx)
                min_y = max(j - r, 0)
                max_y = min(j + r + 1, nely)

                for k in range(min_x, max_x):
                    for l in range(min_y, max_y):

                        # Calculate convolution operator value for the element (k,l) with respect to (i,j)
                        fac = rmin - np.sqrt((i - k)**2 + (j - l)**2)

                        # Accumulate the convolution sum
                        sum_val += max(0, fac)

                        # 基于 convolution 算子的值修改单元的灵敏度
                        dcn[j, i] += max(0, fac) * x[l, k] * dc[l, k]

                # Normalize the modified sensitivity for element (i,j)
                dcn[j, i] /= (x[j, i] * sum_val)

        return dcn

    def OC(self, volfrac, x, dc, passive=None):
        """
        使用 Optimality Criteria Method 执行设计更新.

        Parameters:
        - volfrac (float): Volume fraction, 表示材料将占据的设计空间的期望分数.
        - x (ndarray - (nely, nelx) ): 密度分布矩阵.
        - dc (ndarray - (nely, nelx) ): 原始的灵敏度矩阵.
        - passive (ndarray - (nely, nelx) ): 一个区域，在自由变化的单元处为 0，在固定的单元处为 1.

        Returns:
        - xnew (ndarray - (nely, nelx) ): 优化步骤后更新的设计变量.
        """
        nely, nelx = x.shape

        # 初始化 Lagrange 乘子
        l1, l2 = 0, 1e5

        # 定义设计变量中允许的最大变化
        move = 0.2

        xnew = np.copy(x)

        # 二分法以找到满足体积约束的 Lagrange 乘子
        while (l2 - l1) > 1e-4:
            # 计算当前 Lagrange 乘子区间的中点
            lmid = 0.5 * (l2 + l1)
            # Lower limit move restriction
            tmp0 = x - move
            # Upper limit move restriction
            tmp1 = x + move

            # Design variable update (intermediate step) using OC update scheme
            be = -dc / lmid
            tmp2 = x * np.sqrt(be)
            tmp3 = np.minimum(tmp1, tmp2)
            tmp4 = np.minimum(1, tmp3) 
            tmp5 = np.maximum(tmp0, tmp4)

            xnew = np.maximum(0.001, tmp5)

            # 寻找 Passive 单元，并将其密度设置为最小密度
            if passive is not None:
                xnew[passive == 1] = 0.001

            # 检查当前设计是否满足体积约束
            if np.sum(xnew) - volfrac * nelx * nely > 0:
                l1 = lmid
            else:
                l2 = lmid

        return xnew
