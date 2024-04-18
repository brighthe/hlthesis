import numpy as np
from fealpy.mesh import QuadrangleMesh


class TopLSM:
    def __init__(self, nelx, nely, xlength, yheight):

        self.nelx = nelx
        self.nely = nely
        self.xlength = xlength
        self.yheight = yheight
        domain = [0, self.xlength, 0, self.yheight]
        self.mesh = QuadrangleMesh.from_box(box = domain, nx = self.nelx, ny = self.nely)

    def lk(self):
        """
        刚度矩阵
        这个矩阵通过符号操作得到，是使用有限元方法分析我们的弹性结构的关键。
        """
        E = 1.0  # 弹性模量
        nu = 0.3  # 泊松比
        k = [1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3*nu/8,
            -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3*nu/8]

        KE = E / (1 - nu**2) * np.array([
            [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
            [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
            [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
            [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
            [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
            [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
            [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
            [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
        ])

        return KE

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
        xlength, yheight, mesh = self.xlength, self.yheight, self.mesh
        node = mesh.entity('node')
        nodex = node[:, 0]
        nodey = node[:, 1]
        ls_node = mesh.entity('node')

        # 定义初始孔洞的圆心
        cx = xlength/200 * np.array([33.33,  100,  166.67,   0,    66.67,  133.33,
                                           200,  33.33,  100,   166.67,   0,    66.67,
                                         133.33,  200,  33.33,   100,   166.67], dtype=np.float64)
        cy = yheight/100 * np.array([0,  0,  0,   25,  25, 25,
                                          25, 50, 50,  50,  75, 75,
                                          75, 75, 100, 100, 100], dtype=np.float64)
        # 定义初始孔洞的半径
        cr = yheight / 10

        # 计算每个网格点与所有圆心的欧式距离
        tmpPhi = np.sqrt((nodex[:, None] - cx)**2 + (nodey[:, None] - cy)**2) - cr # (NN, 17)

        # 对于每个节点，取其到所有圆心距离的最小值
        ls_Phi = np.min(tmpPhi, axis=1)
        #print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi)

        # 调整阈值，将材料密度标准化到 -0.1（孔洞） 和 0.1（实体材料）
        ls_Phi = 0.2 * np.ceil(np.maximum(ls_Phi, 0)) - 0.1
        #print("ls_Phi:", ls_Phi.shape, "\n", ls_Phi)


        return ls_Phi


    def mesh0(self, hx, hy, r):
        """
        网格中均匀分布孔洞, 使用三角函数来在结构的设计区域内创建周期性的孔洞模式

        Parameters:
        - hx: 水平方向的孔洞数量
        - hy: 垂直方向的孔洞数量
        - r: 孔洞的尺寸调整参数，范围 (0, 1)，其中较小的值产生小孔洞，较大的值产生大孔洞

        Returns:
        - phi - ( ndarray - (nely+1, nelx+1) ): 设计区域中每点的材料密度，其中 0.1 表示孔洞，-0.1 表示实体材料

        """
        nelx, nely, xlength, yheight, mesh = self.nelx, self.nely, self.xlength, self.yheight, self.mesh
        node = mesh.entity('node')
        nodex = node[:, 0]
        nodey = node[:, 1]

        # 创建一个初始化的设计域，其中所有值为零（即无材料）
        phi0 = np.zeros((nely + 1, nelx + 1))

        # 使用一个 3d 函数将 x 和 y 坐标映射到 z = cos(x)*cos(y)
        # 计算每个节点的材料密度值
        phi0 = -np.cos((hy + 1) * (nodey * np.pi) / yheight) * \
                np.cos((hx + 1) * (nodex * np.pi) / xlength) + r - 1

        # 调整阈值，将材料密度标准化到 -0.1（孔洞） 和 0.1（实体材料）
        phi0 = 0.2 * np.ceil(np.maximum(phi0, 0)) - 0.1
        phi0 = -phi0

        return phi0


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


    def reinitialize(self, phi0, dx, dy, loop_num):
        """ 
        将水平集函数重初始化为符号距离函数.

        parameters:
        - phi0: 演化前的水平集界面;
        - dx: 有限元单元的宽度;
        - dy: 有限元单元的高度;
        - loop_num: 重初始化的时间步数.

        returns:
        - sign_dist_phi (ndarray - (ls_nn, )): 重置化后的水平集界面.
        """
        for _ in range(loop_num + 1):
            # 水平集函数沿 x 和 y 方向的向前和向后差分算子
            dx_l, dx_r = self.upwind_diff(phi0, dx, 'x')
            dy_l, dy_r = self.upwind_diff(phi0, dy, 'y')
            
            # 水平集函数沿 x 和 y 方向的中心差分算子
            dx_c = (dx_l + dx_r) / 2
            dy_c = (dy_l + dy_r) / 2
            
            # sign(phi) 的数值计算
            signphi = phi0 / (np.sqrt(phi0**2 + (dx_c**2 + dy_c**2) * dx**2) + np.finfo(float).eps)

            # cfl 时间步长
            dett = 0.5 * min(dx, dy) / np.max(np.abs(signphi))
            
            grad_plus = np.sqrt(np.maximum(dx_l, 0)**2 + np.minimum(dx_r, 0)**2 +
                                np.maximum(dy_l, 0)**2 + np.minimum(dy_r, 0)**2)
            grad_minus = np.sqrt(np.minimum(dx_l, 0)**2 + np.maximum(dx_r, 0)**2 +
                                 np.minimum(dy_l, 0)**2 + np.maximum(dy_r, 0)**2)

            phi0 = phi0 - dett * \
                (np.maximum(signphi, 0)*grad_plus + np.minimum(signphi, 0)*grad_minus - signphi)

        sign_dist_phi = phi0.flatten('f')

        return sign_dist_phi


    def g(self, u1, u2, v1, v2):
        """
        数值通量函数
        这是 Hamilton-Jacobi 方程数值方案中的通量函数。
        """
        return np.sqrt(np.maximum(u1, 0)**2 + np.minimum(u2, 0)**2 + 
                    np.maximum(v1, 0)**2 + np.minimum(v2, 0)**2)

    def minmod(self, phi1, phi2):
        """
        The minmod limiter function, which is used to prevent numerical
        oscillations in high-resolution shock-capturing schemes.
        
        Args:
        phi1, phi2 (numpy arrays): Input values to compare.
        
        Returns:
        numpy array: Minmod value of phi1 and phi2.
        """
        sphi1 = np.sign(phi1)
        sphi2 = np.sign(phi2)
        mout = np.maximum(0, sphi1 * sphi2) * sphi1 * np.minimum(np.abs(phi1), np.abs(phi2))
        return mout

    def mesh00(self, phi, RIiter):
        """
        重初始化水平集函数
        在初始化之后以及在求解运输方程的过程中周期性地重初始化水平集函数。
        这意味着我们求解方程 d(phi)/dt + sign(phi0)(|grad(phi)|-1) = 0，
        其中 phi(t=0, x) = phi_0(x)。
        我们这样做是为了确保水平集函数在零水平集附近（即我们的结构边界）不会变得过平或过陡。
        """
        xlength, nelx, nely = self.xlength, self.nelx, self.nely
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        cfl = 0.5  # CFL 条件
        dt0 = min(dx, dy) * cfl  # 根据 CFL 条件定义时间步长

        phi = phi.reshape((nelx+1, nely+1)).T
        for _ in range(RIiter):
            phin = np.roll(phi, -1, axis=0)
            phis = np.roll(phi, 1, axis=0)
            phie = np.roll(phi, -1, axis=1)
            phiw = np.roll(phi, 1, axis=1)

            phinn = np.roll(phin, -1, axis=0)
            phiss = np.roll(phis, 1, axis=0)
            phiee = np.roll(phie, -1, axis=1)
            phiww = np.roll(phiw, 1, axis=1)

            dxm = (phi - phie) / dx
            dxp = (phiw - phi) / dx
            dym = (phi - phis) / dy
            dyp = (phin - phi) / dy

            dxmxm = (phi - 2 * phie + phiee) / (dx**2)
            dxpxp = (phiww - 2 * phiw + phi) / (dx**2)
            dxpxm = (phie - 2 * phi + phiw) / (dx**2)

            dymym = (phi - 2 * phis + phiss) / (dy**2)
            dypyp = (phinn - 2 * phin + phi) / (dy**2)
            dypym = (phin - 2 * phi + phis) / (dy**2)

            partA = dxm + .5 * dx * np.minimum(np.abs(dxmxm), np.abs(dxpxm))
            partB = dxp - .5 * dx * np.minimum(np.abs(dxpxp), np.abs(dxpxm))
            partC = dym + .5 * dy * np.minimum(np.abs(dymym), np.abs(dypym))
            partD = dyp - .5 * dy * np.minimum(np.abs(dypyp), np.abs(dypym))

            delp2 = self.g(partA, partB, partC, partD)
            delm2 = self.g(partB, partA, partD, partC)

            nabla = 0.5 * (dxm**2 + dxp**2 + dym**2 + dyp**2)
            sphi = phi / np.sqrt(phi**2 + np.sqrt(dx**2 + dy**2) * nabla / 10)
            sphip = np.maximum(sphi, 0)
            sphim = np.minimum(sphi, 0)

            phi = phi - dt0 * (sphip * delp2 + sphim * delm2 - sphi)

        return phi

    def fe_density(self, phi, eps):
        """
        根据水平集函数在节点的值，为每个矩形有限元分配密度值。
        参数:
        - phi: numpy array, 水平集函数值，维度为 (nely+1, nelx+1)
        - eps: float, 微小材料密度用于近似空洞区域

        返回:
        - FEtheta: numpy array, 每个单元的材料密度，维度为 (nely, nelx)
        """
        nelx, nely = self.nelx, self.nely
        fe_theta = np.zeros((nely, nelx))

        # 处理 phi 中值为 0 的节点
        e = 1e-6
        indices = np.where(phi == 0)
        import random
        for index in zip(*indices):
            phi[index] += (2 * random.randint(0, 1) - 1) * e

        # 遍历每个单元上四个节点处的 phi 值
        for i in range(nely):
            for j in range(nelx):

                phis = [phi[i, j], phi[i, j + 1], phi[i + 1, j + 1], phi[i + 1, j]]
                phis_sign_sum = np.sum(np.sign(phis))

                if phis_sign_sum == 4:
                    fe_theta[i, j] = 1  # 全部节点均为负值，单元完全由材料填充
                elif phis_sign_sum == -4:
                    fe_theta[i, j] = eps  # 全部节点均为正值，单元几乎为空洞
                else:
                    # 通过切割单元的两条主对角线将单元分割成四个三角形
                    # 然后对周围的节点求和，得到第 5 个值，代表单元的中心值
                    # We don't want it to be 0, so we bump it positive or negative 
                    # according to equal probabilities.
                    center_phi = np.sum(phis) / 4
                    if center_phi == 0:
                        center_phi = (2 * random.randint(0, 1) - 1) * e
                    phis.append(center_phi)

                    # 创建三角形矩阵
                    triMat = np.array([
                        [center_phi, phis[0], phis[1]],
                        [center_phi, phis[1], phis[2]],
                        [center_phi, phis[2], phis[3]],
                        [center_phi, phis[3], phis[0]]
                    ])
                    
                    # 处理每个三角形
                    for tri in triMat:
                        tri_sign_sum = np.sum(np.sign(tri))
                        if tri_sign_sum == 3:
                            fe_theta[i, j] += 0.25
                        elif tri_sign_sum == -3:
                            fe_theta[i, j] += 0.25 * eps
                        else:
                            # 使用线性插值估算密度贡献
                            fe_theta[i, j] += self.interpolate_density(tri, eps)

        return fe_theta

    def interpolate_density(self, tri, eps):
        """
        使用线性插值计算三角形的密度贡献。
        """
        n = np.where(np.sign(tri) != np.sign(np.sum(tri)))[0][0]
        phimid = tri[n]
        phicc = tri[(n + 1) % 3]
        phic = tri[(n - 1) % 3]
        f1 = phimid / (phimid - phicc)
        f2 = phimid / (phimid - phic)
        if np.sum(np.sign(tri)) == 1:
            return 0.25 * ((1 - f1 * f2) * eps + f1 * f2)
        else:
            return 0.25 * ((1 - f1 * f2) + f1 * f2 * eps)

    def fe_analysis(self, mesh, fe_theta, KE, F, fixeddofs):
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
        from allaire_operator_integrator import BeamOperatorIntegrator
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
        ldof = vspace[0].number_of_local_dofs()
        vldof = ldof * GD
        integrator = BeamOperatorIntegrator(struc=fe_theta, KE=KE)
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

    def volume(self, fe_theta):
        """
        计算结构的体积
        """
        xlength, nelx = self.xlength, self.nelx
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        # 单元的体积
        dV = dx * dy
        totvol = np.sum(fe_theta) * dV

        return totvol

    def perimeter(self, phi):
        """
        计算结构的周长。
        参数:
        - phi: numpy array, 水平集函数值，维度为 (nely+1, nelx+1)
        - dx, dy: float, 网格在 x 和 y 方向上的步长

        返回:
        - totperim: float, 结构的估计周长
        """
        xlength, nelx = self.xlength, self.nelx
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        epsperim = min(dx, dy) / 20
        sx = phi / np.sqrt(phi**2 + epsperim**2)

        # 计算 sx 在四个方向上的位移
        sxn = np.roll(sx, -1, axis=0)  # 向北移动
        sxs = np.roll(sx, 1, axis=0)   # 向南移动
        sxe = np.roll(sx, -1, axis=1)  # 向东移动
        sxw = np.roll(sx, 1, axis=1)   # 向西移动

        # 处理边界
        sxn[-1, :] = sx[-1, :]  # 北边界
        sxs[0, :] = sx[0, :]    # 南边界
        sxe[:, -1] = sx[:, -1]  # 东边界
        sxw[:, 0] = sx[:, 0]    # 西边界

        # 计算 d(phi)/dx 和 d(phi)/dy
        dsxx = (sxw - sxe) / (2 * dx)
        dsxy = (sxn - sxs) / (2 * dy)

        # 积分并计算周长
        totperim = 0.5 * np.sum(np.sqrt(dsxx**2 + dsxy**2)) * dx * dy

        return totperim

    def regularize(self, phi, V, e2):
        nelx, nely, xlength = self.nelx, self.nely, self.xlength
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        
        k1 = e2 / (dx**2)
        k2 = e2 / (dy**2)

        # Calculate the surface Dirac function
        epsperim = min(dx, dy) / 20
        sx = phi / np.sqrt(phi**2 + epsperim**2)
        
        dsxx = (np.roll(sx, -1, axis=1) - np.roll(sx, 1, axis=1)) / (2 * dx)
        dsxy = (np.roll(sx, -1, axis=0) - np.roll(sx, 1, axis=0)) / (2 * dy)
        
        delta = 0.5 * np.sqrt(dsxx**2 + dsxy**2)
        
        b = -V * delta
        b = b.ravel()  # Flatten the matrix to a vector

        from scipy.sparse import lil_matrix
        from scipy.sparse.linalg import spsolve
        from scipy.linalg import LinAlgError
# Create the regularization matrix M with Neumann conditions
        M = lil_matrix(((nelx+1) * (nely+1), (nelx+1) * (nely+1)))
        for j in range(nelx+1):
            for i in range(nely+1):
                k = j * (nely+1) + i
                M[k, k] = 1 + 2 * k1 + 2 * k2
                
                if i > 0:
                    M[k, k-1] = -k2
                if i < nely:
                    M[k, k+1] = -k2
                if j > 0:
                    M[k, k-(nely+1)] = -k1
                if j < nelx:
                    M[k, k+(nely+1)] = -k1
    # Solve the linear system
        M = M.tocsr()  # Convert to CSR format for solving
        J = spsolve(M, b)

        # Reshape back to matrix form and invert the direction
        v = -J.reshape((nely+1, nelx+1))

        return v


    def curv(self, phi):
        """
        Calculate the mean curvature of a level set function.
        
        Args:
        phi (numpy.ndarray): Level set function.
        dx (float): Grid spacing in the x-direction.
        dy (float): Grid spacing in the y-direction.
        
        Returns:
        numpy.ndarray: Mean curvature of the level set function.
        """
        nelx, nely, xlength = self.nelx, self.nely, self.xlength
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        # Ensure the gradient of phi never goes to 0
        epscurv = min(dx, dy) / 20
        
        # First derivatives to find the gradient of phi
        phin = np.roll(phi, -1, axis=0)  # Shift north
        phis = np.roll(phi, 1, axis=0)   # Shift south
        phie = np.roll(phi, -1, axis=1)  # Shift east
        phiw = np.roll(phi, 1, axis=1)   # Shift west
        
        dphix = (phiw - phie) / (2 * dx)
        dphiy = (phin - phis) / (2 * dy)
        
        # Magnitude of the gradient of phi
        mag = np.sqrt(dphix**2 + dphiy**2 + epscurv**2)
        nx = dphix / mag
        ny = dphiy / mag
        
        # Divergence of the normal vector field
        nxe = np.roll(nx, -1, axis=1)  # East shift of nx
        nxw = np.roll(nx, 1, axis=1)   # West shift of nx
        nyn = np.roll(ny, -1, axis=0)  # North shift of ny
        nys = np.roll(ny, 1, axis=0)   # South shift of ny
        
        divnx = (nxw - nxe) / (2 * dx)
        divny = (nyn - nys) / (2 * dy)
        
        # Mean curvature
        H = divnx + divny
        
        return H

    def solvelevelset(self, phi, V, dt, HJiter, lagP):
        #print("phi:", phi.shape)
        nelx, nely, xlength = self.nelx, self.nely, self.xlength
        dx = xlength / nelx  # x 方向的步长
        dy = dx  # y 方向的步长，与 x 方向相同
        Vp = np.maximum(V, 0)
        Vm = np.minimum(V, 0)
        
        for _ in range(HJiter):
            # 计算一阶导数
            phin = np.roll(phi, -1, axis=0)
            phis = np.roll(phi, 1, axis=0)
            phie = np.roll(phi, -1, axis=1)
            phiw = np.roll(phi, 1, axis=1)
            
            # 计算二阶导数
            dxm = (phi - phie) / dx
            dxp = (phiw - phi) / dx
            dym = (phi - phis) / dy
            dyp = (phin - phi) / dy
            
            dxmxm = (phi - 2 * phie + np.roll(phie, -1, axis=1)) / (dx**2)
            dxpxp = (np.roll(phiw, 1, axis=1) - 2 * phiw + phi) / (dx**2)
            dymym = (phi - 2 * phis + np.roll(phis, 1, axis=0)) / (dy**2)
            dypyp = (np.roll(phin, -1, axis=0) - 2 * phin + phi) / (dy**2)
            
            # 应用通量函数
            partA = dxm + 0.5 * dx * self.minmod(dxmxm, dxm - dxp)
            partB = dxp - 0.5 * dx * self.minmod(dxpxp, dxp - dxm)
            partC = dym + 0.5 * dy * self.minmod(dymym, dym - dyp)
            partD = dyp - 0.5 * dy * self.minmod(dypyp, dyp - dym)

            delp2 = self.g(partA, partB, partC, partD)
            delm2 = self.g(partB, partA, partD, partC)

            # 计算梯度的模 |grad(phi)|
            epscurv = min(dx, dy) / 20
            dphix = (phiw - phie) / (2 * dx)
            dphiy = (phin - phis) / (2 * dy)
            mag = np.sqrt(dphix**2 + dphiy**2 + epscurv**2)

            # 更新水平集函数
            phi = phi - dt * (delp2 * Vp + delm2 * Vm) + dt * lagP * self.curv(phi) * mag

        solvephi = phi
        return solvephi

