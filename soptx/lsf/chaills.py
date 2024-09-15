import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsf:

    def __init__(self, nelx, nely, volReq, stepLength, numReinit, topWeight):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): 沿设计区域水平方向的单元数. nelx > 20, nelx*nely < 5000.
        - nely (int): 沿设计区域垂直方向的单元数. nelx > 20, nelx*nely < 5000.
        - volReq (float) : 最终设计所需的体积分数. 0.2 < volReq < 0.7.
        - stepLength (int): 每次迭代中求解演化方程的 CFL 时间步长.
                            min(nelx,nely)/10 < stepLength < max(nelx,nely)/5.
        - numReinit (int): 水平集函数重置化为符号距离函数的频率. 2 < numReinit < 6.
        - topWeight (int): 演化方程中 forcing 项的权重. 1 < topWeight < 4.

        Note:
            numReinit 和 topWeight 会影响设计中形成孔洞的能力.
        '''

        self._nelx = nelx
        self._nely = nely
        self._domain_x = nelx
        self._domain_y = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit
        self._topWeight = topWeight

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

    def generate_mesh(self, domain, nelx, nely):

        mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)

        return mesh


    def reinit(self, struc):
        """
        根据给定的结构重置化水平集函数.

        该函数通过添加 void 单元的边界来扩展输入结构，计算到最近的 solid 和 void 单元
        的欧几里得距离，并计算水平集函数，该函数在 solid phase 内为正，在 void phase 中为负

        Parameters:
        - struc ( ndarray - (nely, nelx) ): 表示结构的 solid(1) 和 void(0) 单元

        Returns:
        - lsf ( ndarray - (nely+2, nelx+2) ): 表示重置化后的水平集函数
        """
        from scipy import ndimage

        nely, nelx = self._nely, self._nelx
        strucFull = np.zeros((nely + 2, nelx + 2))
        strucFull[1:-1, 1:-1] = struc
        print("strucFull:", strucFull.shape, "\n", strucFull)

        # Compute the distance to the nearest void (0-valued) cells.
        dist_to_0 = ndimage.distance_transform_edt(strucFull)
        print("dist_to_0:", dist_to_0.shape, "\n", dist_to_0)

        # Compute the distance to the nearest solid (1-valued) cells.
        dist_to_1 = ndimage.distance_transform_edt(strucFull - 1)
        print("dist_to_1:", dist_to_1.shape, "\n", dist_to_1.round(2))

        # Offset the distances by 0.5 to center the level set function on the boundaries.
        element_length = self._domain_x / (2*self._nelx)
        temp_0 = dist_to_0 - element_length
        print("temp0:", temp_0.shape, "\n", temp_0)
        temp_1 = dist_to_1 - element_length
        print("temp1:", temp_1.shape, "\n", temp_1.round(2))

        # Calculate the level set function, ensuring the correct sign inside and outside the structure.
        #lsf = -(~strucFull.astype(bool)).astype(int) * temp_1 + strucFull * temp_0
        lsf = -(1 - strucFull) * temp_1 + strucFull * temp_0
        print("lsf", lsf.shape, "\n", lsf.round(2))

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

    def smooth_sens(self, sens, kernel_size=3, padding_mode='edge'):
        """
        Smooth the sensitivity using convolution with a predefined kernel.

        Parameters:
        - sens ( ndarray - (nely, nelx) ): Sensitivity to be smoothed.
        - kernel_size : The size of the convolution kernel. Default is 3.
        - padding_mode : The mode used for padding. Default is 'edge' 
        which pads with the edge values of the array.

        Returns:
        - smoothed_sens ( ndarray - (nely, nelx) ) : Smoothed sensitivity.
        """
        from scipy.signal import convolve2d
        # Convolution filter to smooth the sensitivities
        kernel_value = 1 / (2*kernel_size)
        kernel = kernel_value * np.array([[0, 1, 0], 
                                          [1, 2, 1], 
                                          [0, 1, 0]], dtype=np.int_)

        # Apply padding to the sensitivity array
        padded_sens = np.pad(sens, ((1, 1), (1, 1)), mode=padding_mode)

        # Perform the convolution using the padded array and the kernel
        smoothed_sens = convolve2d(padded_sens, kernel, mode='valid')

        return smoothed_sens

    def updateStep(self, lsf, shapeSens, topSens, stepLength, topWeight, loadBearingIndices):
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
        topSens[loadBearingIndices] = 0

        # 求解水平集函数的演化方程以更新结构
        # 形状灵敏度的负数作为速度场
        # 拓扑灵敏度按 topWeight 因子缩放，仅应用于结构的 solid 部分作为 forcing 项
        struc, lsf = self.evolve(-shapeSens, topSens*(lsf[1:-1, 1:-1] < 0), lsf, stepLength, topWeight)

        return struc, lsf

    def evolve(self, v, g, lsf, stepLength, w):
        """
        求解水平集函数演化方程，以模拟网格上界面的移动
        
        Parameters:
        - v (ndarray - (nely, nelx) ): 表示每个网格点的速度场.
        - g (ndarray - (nely, nelx) ): 表示每个网格点的 forcing 项.
        - lsf (ndarray - (nely+2, nelx+2) ): 表示界面的水平集函数.
        - stepLength (float): The total time for which the level set function should be evolved.
        - w (float): A weighting parameter for the influence of the force term on the evolution.

        Returns:
        - struc ( ndarray - (nely, nelx) ): 演化后的更新结构，只能为 0 或 1.
        - lsf ( ndarray - (nely+2, nelx+2) ): 演化的水平集函数.
        """
        # 用零边界填充速度场和 forcing 项
        vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
        gFull = np.pad(g, ((1,1),(1,1)), mode='constant', constant_values=0)

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
                np.sqrt( np.minimum(dmx, 0)**2 + np.maximum(dpx, 0)**2 + np.minimum(dmy, 0)**2 + np.maximum(dpy, 0)**2 ) \
            - dt * np.maximum(vFull, 0) * \
                np.sqrt( np.maximum(dmx, 0)**2 + np.minimum(dpx, 0)**2 + np.maximum(dmy, 0)**2 + np.minimum(dpy, 0)**2 ) \
            - dt * w * gFull

        # 基于零水平集导出新结构
        strucFULL = (lsf < 0).astype(int)
        struc = strucFULL[1:-1, 1:-1]
        
        return struc, lsf






