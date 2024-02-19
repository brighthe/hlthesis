import numpy as np

from fealpy.mesh import QuadrangleMesh

# Short Cantilever
class TopSimp:
    def __init__(self, nelx: int = 32, nely: int = 20, volfrac: float = 0.4, 
                penal: float = 3.0, rmin: float = 1.2):
        '''
        初始化拓扑优化问题

        Parameters:
        - nelx (int): 水平方向的单元数. Defaults to 32.
        - nely (int): 垂直方向的单元数. Defaults to 20.
        - volfrac (float): Volume fraction, 表示材料将占据的设计空间的期望分数. Defaults to 0.4.
        - penal (float): Penalization power, 在 SIMP 中控制中间密度的 penalization. Defaults to 3.0.
        - rmin (float): Filter 半径. Defaults to 1.2.
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
        self._mesh_top2 = QuadrangleMesh(node=node, cell=cell)

        #self._mesh = QuadrangleMesh.from_box(box = [0, nelx+1, 0, nely+1], nx = nelx, ny = nely)

    def FE(self, mesh, x, penal):
        from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import spdiags

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        vspace = GD*(space, )
        gdof = vspace[0].number_of_global_dofs()
        vgdof = gdof * GD
        ldof = vspace[0].number_of_local_dofs()
        vldof = ldof * GD

        E0 = 1.0
        nu = 0.3
        nely, nelx = x.shape
        integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, penal=penal, x=x)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        bform.assembly()
        K = bform.get_matrix()

        # 定义荷载 - Short Cantilever
        F = np.zeros(vgdof)
        F[vgdof-1] = 1
        #print("F:", F.shape, "\n", F.round(4))

        # 定义支撑(边界处理) - Short Cantilever
        fixeddofs = np.arange(0, 2*(nely+1), 1)
        dflag = fixeddofs
        #print("dflag:", dflag)
        F = F - K@uh.flat
        bdIdx = np.zeros(K.shape[0], dtype=np.int_)
        bdIdx[dflag.flat] = 1
        D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
        D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
        K = D0@K@D0 + D1
        F[dflag.flat] = uh.ravel()[dflag.flat]

        # 线性方程组求解
        uh.flat[:] = spsolve(K, F)

        reshaped_uh = uh.reshape(-1)
        cell2dof = vspace[0].cell_to_dof()
        NC = mesh.number_of_cells()
        updated_cell2dof = np.repeat(cell2dof*GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
        idx = np.array([0, 1, 4, 5, 6, 7, 2, 3], dtype=np.int_)
        # 用 Top 中的自由度替换 FEALPy 中的自由度
        updated_cell2dof = updated_cell2dof[:, idx]
        ue = reshaped_uh[updated_cell2dof] # (NC, ldof*GD)

        return uh, ue

    def check(self, rmin, x, dc):
        """
        应用 mesh-independency filter 进行每个单元的灵敏度过滤.

        Parameters:
            - rmin (float): Filter 半径.
            - x (ndarray): 密度分布矩阵 (nely, nelx).
            - dc (ndarray): 原始的灵敏度矩阵 (nely, nelx).

        Returns:
            - dcn(ndarray): 过滤后的灵敏度矩阵 (nely, nelx).
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

    def OC(self, volfrac, x, dc):
        """
        使用 Optimality Criteria Method 执行设计更新.

        Parameters:
            - volfrac (float): Volume fraction, 表示材料将占据的设计空间的期望分数.
            - x (ndarray): 密度分布矩阵 (nely, nelx).
            - dc (ndarray): 原始的灵敏度矩阵 (nely, nelx).

        Returns:
            - xnew (ndarray): 优化步骤后更新的设计变量 (nely, nelx).
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

            # 检查当前设计是否满足体积约束
            if np.sum(xnew) - volfrac * nelx * nely > 0:
                l1 = lmid
            else:
                l2 = lmid

        return xnew
