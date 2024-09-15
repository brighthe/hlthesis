import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsfShapeGrad:

    def __init__(self, domain_width, domain_hight, nelx, nely, lagV, lagCur):
        '''
        初始化拓扑优化问题 - 形状梯度的水平集方法

        Parameters:
        - domain_width: 设计区域的宽度;
        - domain_hight: 设计区域的高度;
        - nelx (int): x 方向的单元数;
        - nely (int): y 方向的单元数;
        - lagV: Lagrange multiplier for volume constraint;
        - lagCur: Lagrange multiplier for perimeter constraint whose shape sensitivity is curvature;
        '''

        self.domain_width = domain_width
        self.domain_hight = domain_hight
        self.nelx = nelx
        self.nely = nely
        self.lagV = lagV
        self.lagCur = lagCur


    def generate_mesh(self, domain, nelx, nely):

        mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)

        return mesh


    def init_lsf(self, mesh):
        '''
        初始化设计区域内的水平集函数，以分布初始孔洞的符号距离函数表示.

        该函数通过定义一系列圆形孔洞，并计算网格点到这些孔洞边界的最短距离，来初始化水平集函数
        水平集函数在孔洞内部为负，在孔洞外部为正.

        Parameters:
        - mesh (object): 初始的水平集网格.

        Returns:
        - ls_Phi ( ndarray - (ls_NC, ) ): 初始的水平集函数.

        '''
        domain_width = self.domain_width
        domain_hight = self.domain_hight
        ls_node = mesh.entity('node')
        ls_x = ls_node[:, 0]
        ls_y = ls_node[:, 1]

        # 定义初始孔洞的圆心
        cx = domain_width/200 * np.array([33.33,  100,  166.67,   0,    66.67,  133.33,
                                           200,  33.33,  100,   166.67,   0,    66.67,
                                         133.33,  200,  33.33,   100,   166.67], dtype=np.float64)
        cy = domain_hight/100 * np.array([0,  0,  0,   25,  25, 25,
                                          25, 50, 50,  50,  75, 75,
                                          75, 75, 100, 100, 100], dtype=np.float64)
        # 定义初始孔洞的半径
        cr = domain_hight / 10

        # 计算每个网格点与所有圆心的欧式距离
        tmpPhi = -np.sqrt((ls_x[:, None] - cx)**2 + (ls_y[:, None] - cy)**2) + cr

        # 对于每个节点，取其到所有圆心距离的最大负值，以确保节点若在任一孔洞内部，则其水平集函数值为负
        ls_Phi = -np.max(tmpPhi, axis=1)

        return ls_Phi


    def fe_analysis(self, mesh, E0, E1, nu, ew, eh, phi, F, fixeddofs):
        """
        有限元计算位移.

        Parameters:
        - mesh (object): 初始的水平集网格;
        - E0 (float): Young's modulus of void material;
        - E1 (float): Young's modulus of elastic material;
        - nu (float): 泊松比;
        - ew (float): 有限元单元的几何宽度;
        - eh (float): 有限元单元的几何高度;
        - phi ( ndarray - (nely+2, nelx+2) ): 水平集函数在有限元网格节点上的值;
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵;
        - F ( ndarray - (gdof*GD, nLoads) ): 节点荷载;
        - fixeddofs (ndarray): 位移约束(supports).

        Returns:
        - uh ( ndarray - (gdof, GD, nLoads) ): 总位移.
        """
        from shape_gradient_operator_integrator import BeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse import spdiags
        from scipy.sparse.linalg import spsolve

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        nLoads = F.shape[-1]
        uh = np.repeat(uh[:, :, np.newaxis], nLoads, axis=2)
        vspace = GD*(space, )
        e0 = E0
        e1 = E1
        integrator = BeamOperatorIntegrator(E0=e0, E1=e1, nu=nu, ew=ew, eh=eh, phi=phi)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        K = bform.assembly()

        dflag = fixeddofs
        F = F - K@uh.reshape(-1, nLoads)
        bdIdx = np.zeros(K.shape[0], dtype=np.int_)
        bdIdx[dflag.flat] = 1
        D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
        D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
        K = D0@K@D0 + D1
        F[dflag.flat] = uh.reshape(-1, nLoads)[dflag.flat]

        #print("K:", K.shape, "\n", K.toarray())
        #print("F:", F.shape, "\n", F)
        # 线性方程组求解
        uh.flat[:] = spsolve(K, F)

        return uh


    def upwind_diff(self, phi, d, direction):
        """
        使用迎风格式计算向前和向后有限差分.

        Parameters:
        - phi ( ndarray - (nely+2, nlex+2) ): 水平集函数在水平集网格点上的值;
        - d (float): x 或 y 方向上相邻网格的间距, 取决于 'direction' 参数;
        - direction (str): 'x': 沿 x 方向计算差分, 'y': 表示沿 y 方向计算差分.

        Returns:
        - back_diff ( ndarray - (nely+2, nlex+2) ): 向后差分矩阵: ( phi(i, j) - phi(i-1, j) ) / dx
                               或 ( phi(i, j) - phi(i, j-1) ) / dx, 取决于 'direction'.
        - fawd_diff ( ndarray - (nely+2, nlex+2) ): 向前差分矩阵: ( phi(i+1, j) - phi(i, j) ) / dx
                               或 ( phi(i, j+1) - phi(i, j) ) / dx, 取决于 'direction'.
        """

        # 根据指定的方向计算后向和前向差分
        if direction == 'x':
            x_minus_1 = np.roll(phi, 1, axis=1) # x 方向向右位移
            x_plus_1 = np.roll(phi, -1, axis=1) # x 方向向左位移
            dx = d
            back_diff = (phi - x_minus_1) / dx
            fawd_diff = (x_plus_1 - phi) / dx
        elif direction == 'y':
            y_minus_1 = np.roll(phi, 1, axis=0) # y 方向向下位移
            y_plus_1 = np.roll(phi, -1, axis=0) # y 方向向上位移
            dy = d
            back_diff = (phi - y_minus_1) / dy
            fawd_diff = (y_plus_1 - phi) / dy
        else:
            raise ValueError("direction 必须是 'x' 或 'y'")
        
        return back_diff, fawd_diff


    def calc_curvature(self, phi, dx, dy):
        """
        计算水平集函数的曲率.

        Parameters:
        - phi ( ndarray - (nely+2, nlex+2) ): 水平集函数在水平集网格点上的值;
        - dx (float): x 轴方向上相邻网格的间距;
        - dy (float): y 轴方向上相邻网格的间距.

        Returns:
        - ls_curv ( ndarray - (ls_NC, ) ): 水平集网格点处的平均曲率值.
        """
        # 计算 phi 关于 x 和 y 的一阶偏导数
        x_minus_1 = np.roll(phi, 1, axis=1) # x 方向向右位移
        x_plus_1 = np.roll(phi, -1, axis=1) # x 方向向左位移
        phix = (x_plus_1 - x_minus_1) / (2 * dx)
        #print("phix:\n", phix)
        y_minus_1 = np.roll(phi, 1, axis=0) # y 方向向下位移
        y_plus_1 = np.roll(phi, -1, axis=0) # y 方向向上位移
        phiy = (y_plus_1 - y_minus_1) / (2 * dy)

        # 使用迎风差分格式计算 phiy 关于 x 的偏导数
        phiyx_bk, phiyx_fw = self.upwind_diff(phiy, dx, 'x')
        phixy = (phiyx_bk + phiyx_fw) / 2
        
        # 计算 phi 的二阶偏导数
        phixx = (x_plus_1 - 2*phi + x_minus_1) / dx**2
        phiyy = (y_plus_1 - 2*phi + y_minus_1) / dy**2
        
        # 根据曲率公式计算每个点的曲率
        curvature = (phixx * phiy**2 - 2*phix*phiy*phixy + phiyy * phix**2) / \
                    (np.power(phix**2 + phiy**2, 1.5) + 100*np.finfo(float).eps)
        #print("curvature:\n", curvature)
        ls_curv = curvature.flatten('F')
        
        return ls_curv


    def sensi_analysis(self, fe_mesh, ls_mesh, E1, E0, u, v, ew, eh, nu, lag4Vol, lag4Curv,
                             ele_lsgrid_id, phi, curvature,):
        """
        计算应变能和特定水平集函数网格上的速度场.

        Parameters:
        - E1 (float): Young's modulus of elastic material;
        - E0 (float): Young's modulus of void material;
        - u: 网格节点的位移向量的 x 分量;
        - v: 网格节点的位移向量的 y 分量;
        - ew (float): 有限元单元的几何宽度;
        - eh (float): 有限元单元的几何高度;
        - nu (float): 泊松比;
        - lag4Vol (float): Lagrange multiplier for volume constraint;
        - lag4Curv (float): Lagrange multiplier for perimeter constraint whose shape sensitivity is curvature;
        - ele_lsgrid_id ( ndarray - (ls_NC, ) ): 有限元网格单元中心对应的水平集网格节点的索引;
        - phi ( ndarray - (nely+2, nlex+2) ): 水平集函数在水平集网格点上的值;
        - curvature ( ndarray - (ls_NC, ) ): 水平集函数的曲率.

        Returns:
        - beta ( ndarray - (ls_NC, ) ): 计算得到的速度场.
        """
        from fealpy.functionspace import LagrangeFESpace as Space

        p = 1
        space = Space(fe_mesh, p=p, doforder='vdims')
        GD = 2
        vspace = GD*(space, )
        ldof = vspace[0].number_of_local_dofs()
        vldof = ldof * GD
        fe_cell2node = fe_mesh.ds.cell_to_node()
        fe_NC = fe_mesh.number_of_cells()

        # 计算单元内每个自由度处的位移
        ae = np.zeros((fe_NC, vldof))
        ae[:, ::2] = u[fe_cell2node].reshape(fe_NC, ldof) # 偶数索引处填充 u
        ae[:, 1::2] = v[fe_cell2node].reshape(fe_NC, ldof) # 奇数索引处填充 v
        #print("ae:", ae.shape, "\n", ae)

        # 构建应变矩阵 B
        ew = ew
        eh = eh
        B = 1/2 * np.array([[-1/ew,   0,   1/ew,   0,   1/ew,  0,   -1/ew,  0],
                            [  0,   -1/eh,  0,   -1/eh,  0,   1/eh,   0,   1/eh],
                            [-1/eh, -1/ew, -1/eh, 1/ew, 1/eh, 1/ew, 1/eh, -1/ew]], dtype=np.float64)
        #print("B:", B.shape, "\n", B)
        
        # 计算应变
        strain = np.einsum('ij, kj -> ki', B, ae)
        #print("strain:", strain.shape, "\n", strain)
        
        ls_NC = ls_mesh.number_of_nodes()
        beta = np.zeros(ls_NC)
        phi = phi[ele_lsgrid_id]
        curvature = curvature[ele_lsgrid_id]

        # 根据 phi 值确定材料的弹性模量 E
        E = np.zeros_like(phi)

        delta = 0.75 * ew
        E[phi > delta] = E1
        E[phi < -delta] = E0

        # 对于界面区域内的材料属性进行插值
        density_min = 1e-3
        mask = (phi <= delta) & (phi >= -delta)
        xd = phi[mask] / delta
        E[mask] = E1 * 0.75 * (1.0 - density_min) * (xd - xd**3/3.0) + 0.5 * (1 + density_min)
        #print("E:", E.shape, "\n", E)

        # 构建本构矩阵 D
        D = E[:, np.newaxis, np.newaxis] / (1 - nu**2) * \
                              np.array([[1,  nu,    0],
                                        [nu, 1,     0],
                                        [0,  0, (1-nu)/2]], dtype=np.float64)
        #print("D:", D.shape, "\n", D)

        # 计算变形能量
        strain_energy = np.einsum('ki, kij, kj -> k', strain, D, strain)
        #print("strain_energy:", strain_energy.shape, "\n", strain_energy)

        # 计算速度场 beta
        beta[ele_lsgrid_id] = lag4Vol - strain_energy - lag4Curv * curvature

        return beta


    def level_set_evolve(self, phi0, vn, dx, dy, loop_num):
        """
        使用一阶 space convex 更新水平集界面.

        Parameters:
        - phi0: 演化前的水平集界面;
        - vn: 法向速度场;
        - dx: 有限元单元的宽度;
        - dy: 有限元单元的高度;
        - loop_num: 演化的循环次数.

        Returns:
        - phi1 (ndarray - (ls_NN, )): 演化后的水平集界面.
        """
        detT = 0.5 * min(dx, dy) / np.max(np.abs(vn))
        for _ in range(loop_num):
            dx_L, dx_R = self.upwind_diff(phi0, dx, 'x')
            dy_L, dy_R = self.upwind_diff(phi0, dy, 'y')

            grad_Plus = np.sqrt(np.maximum(dx_L, 0)**2 + np.minimum(dx_R, 0)**2 +
                                np.maximum(dy_L, 0)**2 + np.minimum(dy_R, 0)**2)
            grad_Minus = np.sqrt(np.minimum(dx_L, 0)**2 + np.maximum(dx_R, 0)**2 +
                                 np.minimum(dy_L, 0)**2 + np.maximum(dy_R, 0)**2)

            phi0 = phi0 - detT * (np.maximum(vn, 0) * grad_Plus + np.minimum(vn, 0) * grad_Minus)

        phi1 = phi0.flatten('F')

        return phi1


    def reinitialize(self, phi0, dx, dy, loop_num):
        """
        将水平集函数重初始化为符号距离函数.

        Parameters:
        - phi0: 演化前的水平集界面;
        - dx: 有限元单元的宽度;
        - dy: 有限元单元的高度;
        - loop_num: 重初始化的时间步数.

        Returns:
        - sign_dist_phi (ndarray - (ls_NN, )): 重置化后的水平集界面.
        """
        for _ in range(loop_num + 1):
            # 水平集函数沿 x 和 y 方向的向前和向后差分算子
            dx_L, dx_R = self.upwind_diff(phi0, dx, 'x')
            dy_L, dy_R = self.upwind_diff(phi0, dy, 'y')
            
            # 水平集函数沿 x 和 y 方向的中心差分算子
            dx_C = (dx_L + dx_R) / 2
            dy_C = (dy_L + dy_R) / 2
            
            # sign(Phi) 的数值计算
            signPhi = phi0 / (np.sqrt(phi0**2 + (dx_C**2 + dy_C**2) * dx**2) + np.finfo(float).eps)

            # CFL 时间步长
            detT = 0.5 * min(dx, dy) / np.max(np.abs(signPhi))
            
            grad_plus = np.sqrt(np.maximum(dx_L, 0)**2 + np.minimum(dx_R, 0)**2 +
                                np.maximum(dy_L, 0)**2 + np.minimum(dy_R, 0)**2)
            grad_minus = np.sqrt(np.minimum(dx_L, 0)**2 + np.maximum(dx_R, 0)**2 +
                                 np.minimum(dy_L, 0)**2 + np.maximum(dy_R, 0)**2)

            phi0 = phi0 - detT * \
                (np.maximum(signPhi, 0)*grad_plus + np.minimum(signPhi, 0)*grad_minus - signPhi)

        sign_dist_phi = phi0.flatten('F')

        return sign_dist_phi










