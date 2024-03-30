import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsfShapeGrad:

    def __init__(self, domain_width, domain_hight, nelx, nely, volReq, stepLength, numReinit):
        '''
        初始化拓扑优化问题 - 形状梯度的水平集方法

        Parameters: 
        - domain: 设计区域.
        - nelx (int): 沿设计区域水平方向的单元数.
        - nely (int): 沿设计区域垂直方向的单元数.
        - volReq (float) : 最终设计所需的体积分数. 0.2 < volReq < 0.7.
        - stepLength (int): 每次迭代中求解演化方程的 CFL 时间步长.
                            min(nelx,nely)/10 < stepLength < max(nelx,nely)/5.
        - numReinit (int): 水平集函数重置化为符号距离函数的频率. 2 < numReinit < 6.

        Note:
            numReinit 和 topWeight 会影响设计中形成孔洞的能力.
        '''

        self.domain_width = domain_width
        self.domain_hight = domain_hight
        self.nelx = nelx
        self.nely = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit

    def generate_mesh(self, domain, nelx, nely):

        mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)

        return mesh

    def plot(self, x, y, z, label):

        import matplotlib.pyplot as plt

        plt.scatter(x, y, c=z, cmap='viridis', label=label, edgecolors='k', s=50)

        plt.colorbar(label='Phi values')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot_mesh(self, x0, y0, label0, x, y, z, label):

        import matplotlib.pyplot as plt

        plt.scatter(x, y, c=z, cmap='cool', label=label, edgecolors='r', marker='s', s=50)

        plt.scatter(x0, y0, color='green', marker='o', label=label0)

        plt.colorbar(label='Phi values')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.legend()
        plt.grid(True)
        plt.show()


    def init_lsf(self, mesh):
        '''
        初始化设计区域内的水平集函数，以分布初始孔洞的符号距离函数表示.

        该函数通过定义一系列圆形孔洞，并计算网格点到这些孔洞边界的最短距离，来初始化水平集函数
        水平集函数在孔洞内部为负，在孔洞外部为正

        Parameters:
        - mesh (object): 初始的水平集网格.

        Returns:
        - ls_Phi ( ndarray - ( (domain_width+1)*(domain_hight+1), ) ): 初始的水平集函数.

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


        #def basic_ke(self, E, nu, a, b):
        #    """
        #    Returns the element stiffness matrix of a full/empty element.

        #    采用平面应力假设，根据二维线弹性理论计算局部刚度矩阵.
        #    
        #    Parameters:
        #    - E : Young's modulus, 代表材料的刚性
        #    - nu: Poisson ratio;
        #    - a: 有限元单元的几何宽度.
        #    - b: 有限元单元的几何高度.
        #    
        #    Returns:
        #    - KE : a 8-by-8 stiffness matrix.
        #    """
        #    # 刚度矩阵的系数
        #    k = np.array([-1/(6*a*b) * ( nu*a**2 - 2*b**2 - a**2 ), 1/8*nu + 1/8,
        #                -1/(12*a*b) * ( nu*a**2 + 4*b**2 - a**2 ), 3/8*nu - 1/8,
        #                1/(12*a*b) * ( nu*a**2 - 2*b**2 - a**2 ), -1/8*nu - 1/8,
        #                1/(6*a*b) * ( nu*a**2 + b**2 - a**2 ), -3/8*nu + 1/8])
        #    
        #    KE = E / (1 - nu**2) * np.array([
        #        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
        #        [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
        #        [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        #        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        #        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        #        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        #        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        #        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
        #    ])

        #    return KE

        #def ele_stiff_matrix(self, ew, eh, E1, E0, nu, Phi, ele_node_id, i):
        #    """
        #    返回单元刚度矩阵，根据它们与边界的相对位置，包含 3 种情况：
        #    inside the boundary, outside the boundary, or on the boundary.

        #    Parameters:
        #    - ew: 有限元单元的几何宽度.
        #    - eh: 有限元单元的几何高度.
        #    - E1: Young's modulus of elastic material.
        #    - E0: Young's modulus of void material.
        #    - nu: Poisson ratio, assuming the two materials have the same nu.
        #    - Phi: 水平集函数在有限元节点处的值.
        #    - ele_node_id: The corresponding serial number of the 4 nodes in a finite element.
        #    - i: 当前处理的单元索引.

        #    Returns:
        #    - ke: Element stiffness matrix.
        #    """
        #    if np.min( Phi[ ele_node_id[i] ] ) > 0:  # Element is inside the boundary
        #        ke = self.basic_ke(E1, nu, ew, eh)

        #    elif np.max( Phi[ ele_node_id[i] ] ) < 0:  # Element is outside the boundary
        #        ke = self.basic_ke(E0, nu, ew, eh)

        #    else:  # Element is cut by the boundary
        #        s, t = np.meshgrid(np.linspace(-1, 1, 21), np.linspace(-1, 1, 21))
        #        print("ele_node_id:", ele_node_id.shape, "\n", ele_node_id)

        #        # 插值水平集函数, 将每个节点的水平集函数值与其相应的形函数值相乘，然后对所有节点求和
        #        c0 = ( (1 - s.flatten('F')) * (1 - t.flatten('F')) )[:, np.newaxis]
        #        c1 = ( (1 + s.flatten('F')) * (1 - t.flatten('F')) )[:, np.newaxis]
        #        c2 = ( (1 + s.flatten('F')) * (1 + t.flatten('F')) )[:, np.newaxis]
        #        c3 = ( (1 - s.flatten('F')) * (1 + t.flatten('F')) )[:, np.newaxis]
        #        tmpPhi = c0 / 4 * Phi.flatten('F')[ele_node_id[:, 0]] + \
        #                 c1 / 4 * Phi.flatten('F')[ele_node_id[:, 1]] + \
        #                 c2 / 4 * Phi.flatten('F')[ele_node_id[:, 2]] + \
        #                 c3 / 4 * Phi.flatten('F')[ele_node_id[:, 3]]

        #        #tmpPhi = ((1 - s.flatten('F')) * (1 - t.flatten('F')) / 4 * Phi[ele_node_id[i, 0]] + 
        #        #          (1 + s.flatten('F')) * (1 - t.flatten('F')) / 4 * Phi[ele_node_id[i, 1]] +
        #        #          (1 + s.flatten('F')) * (1 + t.flatten('F')) / 4 * Phi[ele_node_id[i, 2]] +
        #        #          (1 - s.flatten('F')) * (1 + t.flatten('F')) / 4 * Phi[ele_node_id[i, 3]])
        #        print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))

        #        # 计算覆盖材料区域的面积比
        #        area_ratio = np.sum(tmpPhi >= 0) / len(s.flatten('F'))
        #        ke = area_ratio * self.basic_ke(E1, nu, ew, eh)
        #    
        #    return ke




    def fe_analysis(self, mesh, E0, E1, nu, ew, eh, Phi, F, fixeddofs):
        """
        有限元计算位移.

        Parameters:
        - mesh:
        - struc ( ndarray - (nely, nelx) ): 表示结构的 solid(1) 和 void(0) 单元.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        - F ( ndarray - (gdof*GD, nLoads) ): 节点荷载.
        - fixeddofs (ndarray): 位移约束(supports).

        Returns:
        - uh ( ndarray - (gdof, GD, nLoads) ): 总位移.
        - ue ( ndarray - (NC, ldof*GD, nLoads) ): 单元位移.
        """
        from shape_gradient_operator_integrator import BeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse import spdiags

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        nLoads = F.shape[-1]
        uh = np.repeat(uh[:, :, np.newaxis], nLoads, axis=2)
        vspace = GD*(space, )
        e0 = E0
        e1 = E1
        phi = Phi
        integrator = BeamOperatorIntegrator(E0=e0, E1=e1, nu=nu, ew=ew, eh=eh, Phi=phi)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        #print("KK:", KK.shape, "\n", KK)
        K = bform.assembly()
        print("K0:", K.shape, "\n", K.toarray())

        dflag = fixeddofs
        F = F - K@uh.reshape(-1, nLoads)
        bdIdx = np.zeros(K.shape[0], dtype=np.int_)
        bdIdx[dflag.flat] = 1
        D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
        D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
        K = D0@K@D0 + D1
        print("K1:", K.shape, "\n", K.toarray())
        F[dflag.flat] = uh.reshape(-1, nLoads)[dflag.flat]

        ## 线性方程组求解
        #uh.flat[:] = spsolve(K, F)


    import numpy as np

    def Assemble(self, KE, ke, elementsNodeID, EleID):
        """
        Assembles the element stiffness matrix ke of the quadrilateral element
        into the global stiffness matrix K.

        Parameters:
        - K : 全局刚度矩阵.
        - KE : 单元刚度矩阵.
        - elementsNodeID: Array storing the node IDs for a specified element, shape (TotalEleNum, 4).
        - EleID: The serial number of a specified finite element.

        Returns:
        - K: Updated global stiffness matrix after assembling the element stiffness matrix.
        """
        m = elementsNodeID[EleID, :]
        for i in range(len(m)):
            for j in range(len(m)):
                K[2*m[i]-2 : 2*m[i], 2*m[j]-2: 2*m[j]] += ke[2*i-2: 2*i, 2*j-2:2*j]

        return K










    def lsf_init(self, mesh, r):
        '''
        水平集函数初始化为设计区域内分布初始孔洞的符号距离函数.

        Parameters:
        - mesh:
        - r (float): 初始孔洞的半径.

        Returns:
        - Phi ( ndarray - (nely+1, nelx+1) ): 初始的水平集函数，值在 (-3, 3) 之间.

        '''
        nelx = self._nelx
        nely = self._nely

        node = mesh.entity('node') # 按列增加
        # 网格中节点的 x 坐标 - (nely+1, nelx+1)
        X = node[:, 0].reshape(nelx+1, nely+1).T
        #print("X:", X.shape, "\n", X)
        # 网格中节点的 y 坐标 - (nely+1, nelx+1)
        Y = node[:, 1].reshape(nelx+1, nely+1).T
        #print("Y:", Y.shape, "\n", Y)

        # hX 是初始孔洞的中心处的 x 坐标 - (15, )
        hX = nelx * np.concatenate([np.tile([1/6, 5/6], 3), np.tile([0, 1/3, 2/3, 1], 2), [1/2]])
        #print("hX:", hX.shape, "\n", hX)
        # hY 是初始孔洞的中心处的 y 坐标 - (15, )
        hY = nely * np.concatenate([np.repeat([0, 1/2, 1], 2), np.repeat([1/4, 3/4], 4), [1/2]])
        #print("hY:", hY.shape, "\n", hY)

        # dX 是所有网格点在 x 方向上与初始孔洞的中心之间的距离差, 形状为 (nely+1, nelx+1, 15)
        dX = X[:, :, np.newaxis] - hX
        #print("dX:", dX.shape, "\n", dX)
        # dY 是所有网格点在 Y 方向上与初始孔洞的中心之间的距离差, 形状为 (nely+1, nelx+1, 15)
        dY = Y[:, :, np.newaxis] - hY
        #print("dY:", dY.shape, "\n", dY)

        # 计算网格点到最近孔洞附近的欧氏距离，并限制在 -3 到 3 之间
        Phi = np.sqrt(dX**2 + dY**2) - r
        Phi = np.min(Phi, axis=2)
        Phi = np.clip(Phi, -3, 3)

        return Phi

    def reinit(self, struc):
        """
        根据给定的结构重置化水平集函数.

        该函数通过添加 void 单元的边界来扩展输入结构，计算到最近的 solid 和 void 单元
        的欧几里得距离，并计算水平集函数，该函数在 solid phase 内为负，在 void phase 中为正.

        Parameters:
        - struc ( ndarray - (nely, nelx) ): 表示结构的 solid(1) 和 void(0) 单元.

        Returns:
        - lsf ( ndarray - (nely+2, nelx+2) ): 表示重置化后的水平集函数.
        """
        from scipy import ndimage

        nely, nelx = struc.shape
        strucFull = np.zeros((nely + 2, nelx + 2))
        strucFull[1:-1, 1:-1] = struc

        # Compute the distance to the nearest void (0-valued) cells.
        dist_to_0 = ndimage.distance_transform_edt(strucFull)

        # Compute the distance to the nearest solid (1-valued) cells.
        dist_to_1 = ndimage.distance_transform_edt(strucFull - 1)

        # Offset the distances by 0.5 to center the level set function on the boundaries.
        temp_1 = dist_to_1 - 0.5
        temp_2 = dist_to_0 - 0.5

        # Calculate the level set function, ensuring the correct sign inside and outside the structure.
        lsf = (~strucFull.astype(bool)).astype(int) * temp_1 - strucFull * temp_2

        return lsf

    def FE(self, mesh, struc, KE, F, fixeddofs):
        """
        有限元计算位移.

        Parameters:
        - mesh:
        - struc ( ndarray - (nely, nelx) ): 表示结构的 solid(1) 和 void(0) 单元.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        - F ( ndarray - (gdof*GD, nLoads) ): 节点荷载.
        - fixeddofs (ndarray): 位移约束(supports).

        Returns:
        - uh ( ndarray - (gdof, GD, nLoads) ): 总位移.
        - ue ( ndarray - (NC, ldof*GD, nLoads) ): 单元位移.
        """
        from lsf_beam_operator_integrator import BeamOperatorIntegrator
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

        integrator = BeamOperatorIntegrator(struc=struc, KE=KE)
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


    def updateStep(self, lsf, shapeSens, stepLength, loadBearingIndices):
        """
        使用形状灵敏度和拓扑灵敏度执行设计更新
        
        Parameters:
        - lsf ( ndarray - (nely+2, nelx+2) ): 水平集函数，描述当前结构的界面.
        - shapeSens ( ndarray - (nely, nelx) ): 当前设计的形状灵敏度.
        - topSens ( ndarray - (nely, nelx) ): 当前设计的拓扑灵敏度.
        - stepLength (float): The step length parameter controlling the extent of the evolution.
        - topWeight (float): The weighting factor for the topological sensitivity in the evolution.
        - loadBearingIndices (tuples): 荷载支撑单元的索引.

        Returns:
        - struc ( ndarray - (nely, nelx) ): 设计更新后的新结构，只能为 0 或 1.
        - lsf ( ndarray - (nely+2, nelx+2) ): 设计更新后的新水平集函数.
        """

        # Load bearing pixels must remain solid - short cantilever
        shapeSens[loadBearingIndices] = 0

        # 求解水平集函数的演化方程以更新结构
        # 形状灵敏度的负数作为速度场
        struc, lsf = self.evolve(-shapeSens, lsf, stepLength)

        return struc, lsf

    def evolve(self, v, lsf, stepLength):
        """
        求解水平集函数演化方程，以模拟网格上界面的移动
        
        Parameters:
        - v (ndarray - (nely, nelx) ): 表示每个网格点的速度场.
        - lsf (ndarray - (nely+2, nelx+2) ): 表示界面的水平集函数.
        - stepLength (float): The total time for which the level set function should be evolved.
        - w (float): A weighting parameter for the influence of the force term on the evolution.

        Returns:
        - struc ( ndarray - (nely, nelx) ): 演化后的更新结构，只能为 0 或 1.
        - lsf ( ndarray - (nely+2, nelx+2) ): 演化的水平集函数.
        """
        # 用零边界填充速度场和 forcing 项
        vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
        print("vFull", vFull.shape)
        print("lsf", lsf.shape)

        # 基于 CFL 值选择演化的时间步
        frac_time_step = 0.1 # Fraction of the CFL time step to use as a time step 
                        # for solving the evolution equation
        dt = frac_time_step / np.max(np.abs(v))

        # 基于演化方程迭代更新水平集函数
        num_time_step = 1 / frac_time_step # Number of time steps required to evolve
                        # the level-set function for a time equal to the CFL time step
        for _ in range(int(num_time_step * stepLength)):
            # 计算 x 方向和 y 方向的向前和向后差分
            dpx = np.roll(lsf, shift=(0, -1), axis=(0, 1)) - lsf # forward differences in x directions
            dmx = lsf - np.roll(lsf, shift=(0, 1), axis=(0, 1)) # backward differences in x directions
            dpy = np.roll(lsf, shift=(-1, 0), axis=(0, 1)) - lsf # forward differences in y directions
            dmy = lsf - np.roll(lsf, shift=(1, 0), axis=(0, 1)) # backward differences in y directions

            # 使用迎风格式更新水平集函数
            lsf = lsf \
                - dt * np.minimum(vFull, 0) * \
                np.sqrt( np.minimum(dmx, 0)**2 + np.maximum(dpx, 0)**2 + \
                         np.minimum(dmy, 0)**2 + np.maximum(dpy, 0)**2 ) \
                - dt * np.maximum(vFull, 0) * \
                np.sqrt( np.maximum(dmx, 0)**2 + np.minimum(dpx, 0)**2 + \
                        np.maximum(dmy, 0)**2 + np.minimum(dpy, 0)**2 )

        # 基于零水平集导出新结构
        strucFULL = (lsf < 0).astype(int)
        struc = strucFULL[1:-1, 1:-1]
        
        return struc, lsf






