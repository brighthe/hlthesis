import numpy as np

from fealpy.mesh import QuadrangleMesh
from numpy.ma import shape

class TopLsf:

    def __init__(self, nelx: int = 32, nely: int = 20, volReq: float = 0.4, 
                stepLength: int = 2, numReinit: int = 3, topWeight: int = 2):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): 水平方向的单元数. Defaults to 32.
        - nely (int): 垂直方向的单元数. Defaults to 20.
        - volReq (float) : 最终设计所需的体积分数. Defaults to 0.4.
        - stepLength (int): 演化方程每次迭代中使用的 CFL 时间步长. Defaults to 2.
        - numReinit (int): The number of algorithm iterations before performing level set reinitialization. Defaults to 3.
        - topWeight (int): 演化方程中 forcing 项的权重. Defaults to 2.
        '''

        self._nelx = nelx
        self._nely = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit
        self._topWeight = topWeight
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

        self._mesh2 = QuadrangleMesh.from_box(box = [0, self._nelx, 0, self._nely], \
                                              nx = self._nelx, ny = self._nely)

    def reinit(self, struc):
        """
        根据给定的结构重置化水平集函数

        该函数通过添加 void 单元的边界来扩展输入结构，计算到最近的 solid 和 void 单元
        的欧几里得距离，并计算水平集函数，该函数在 solid phase 内为负，在 void phase 中为正.

        Parameters:
        - struc (numpy.ndarray): 表示结构的 solid(1) 和 void(0) 单元.

        Returns:
        - lsf (ndarray): A 2D array of the same shape as 'struc', 表示重置化后的水平集函数
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

    def FE(self, mesh, struc):
        from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse.linalg import spsolve
        from scipy.sparse import spdiags

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        #print("uh:", uh.shape)
        vspace = GD*(space, )
        gdof = vspace[0].number_of_global_dofs()
        vgdof = gdof * GD
        ldof = vspace[0].number_of_local_dofs()
        vldof = ldof * GD
        #print("vgdof", vgdof)
        #print("vldof", vldof)

        E0 = 1.0
        nu = 0.3
        nely, nelx = struc.shape
        integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        #print("KK:", KK.shape, "\n", KK.round(4))
        bform.assembly()
        K = bform.get_matrix()
        #print("K:", K.shape, "\n", K.toarray().round(4))

        #import scipy.io as io
        #K_test = io.loadmat('KFull.mat')
        #K2_test = K_test['kFull']
        #error = np.sum(np.abs(K2-K2_test))
        #print("error:", error)


        # 定义荷载 - Cantilever MBB
        F = np.zeros(vgdof)
        F[vgdof-1] = 1
        #print("F:", F.shape, "\n", F.round(4))

        #from scipy.sparse import csr_matrix
        #K2_test = csr_matrix(K2_test)
        #K = K2_test
        # 定义支撑(边界处理) - Cantilever MBB
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
        #print("uh:", uh.shape, "\n", uh)

        reshaped_uh = uh.reshape(-1)
        cell2dof = vspace[0].cell_to_dof()
        NC = mesh.number_of_cells()
        updated_cell2dof = np.repeat(cell2dof*GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
        idx = np.array([0, 1, 4, 5, 6, 7, 2, 3], dtype=np.int_)
        # 用 Top 中的自由度替换 FEALPy 中的自由度
        updated_cell2dof = updated_cell2dof[:, idx]
        ue = reshaped_uh[updated_cell2dof] # (NC, ldof*GD)

        return uh, ue

    def smooth_sens(self, sens, kernel_size=3, padding_mode='edge'):
        """
        Smooth the sensitivity using convolution with a predefined kernel.

        Parameters:
        - sens : Sensitivity to be smoothed.
        - kernel_size : The size of the convolution kernel. Default is 3.
        - padding_mode : The mode used for padding. Default is 'edge' 
        which pads with the edge values of the array.

        Returns:
        - Smoothed sensitivity.
        """
        from scipy.signal import convolve2d
        # 定义 convolution kernel
        kernel_value = 1 / (2*kernel_size)
        kernel = kernel_value * np.array([[0, 1, 0], 
                                          [1, 2, 1], 
                                          [0, 1, 0]], dtype=np.int_)

        # Apply padding to the sensitivity array
        padded_sens = np.pad(sens, ((1, 1), (1, 1)), mode=padding_mode)

        # Perform the convolution using the padded array and the kernel
        smoothed_sens = convolve2d(padded_sens, kernel, mode='valid')

        return smoothed_sens

    def updateStep(self, lsf, shapeSens, topSens, stepLength, topWeight):
        """
        使用形状灵敏度和拓扑灵敏度执行设计更新
        
        Parameters:
        - lsf (ndarray): The level set function, which describes the interface of the current structure.
        - shapeSens : 当前设计的形状灵敏度.
        - topSens : 当前设计的拓扑灵敏度.
        - stepLength (float): The step length parameter controlling the extent of the evolution.
        - topWeight (float): The weighting factor for the topological sensitivity in the evolution.

        Returns:
        - struc (numpy.ndarray): The updated structure after the design step, represented as a 2D array.
        - lsf (numpy.ndarray): The updated level set function after the design step.
        """

        # Load bearing pixels must remain solid - Cantilever MBB
        shapeSens[-1, -1] = 0
        #print("shapeSens4:", shapeSens.shape, "\n", shapeSens)
        topSens[-1, -1] = 0
        #print("topSens4:", topSens.shape, "\n", topSens)

        # 求解水平集函数的演化方程以更新结构
        # 形状灵敏度的负数作为速度场
        # 拓扑灵敏度按 topWeight 因子缩放，仅应用于结构的 solid 部分作为 forcing 项
        #print("lsf_origin:\n", lsf.round(4))
        struc, lsf = self.evolve(-shapeSens, topSens*(lsf[1:-1, 1:-1] < 0), lsf, stepLength, topWeight)

        return struc, lsf

    def evolve(self, v, g, lsf, stepLength, w):
        """
        求解水平集函数演化方程，以模拟网格上界面的移动
        
        Parameters:
        - v (ndarray): 表示每个网格点的速度场.
        - g (ndarray): 表示每个网格点的 forcing 项.
        - lsf (ndarray): 表示界面的水平集函数.
        - stepLength (float): The total time for which the level set function should be evolved.
        - w (float): A weighting parameter for the influence of the force term on the evolution.

        Returns:
        - struc (ndarray): 演化后的更新结构，只能为 0 或 1.
        - lsf (ndarray): 演化的水平集函数.
        """
        # 用零边界填充速度场和 forcing 项
        #print("v:", v)
        #print("g:", g)
        vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
        gFull = np.pad(g, ((1,1),(1,1)), mode='constant', constant_values=0)
        #print("vFull:", vFull)
        #print("gFull:", gFull)

        # 基于 CFL 值选择演化的时间步
        dt = 0.1 / np.max(np.abs(v))

        # 基于演化方程迭代更新水平集函数
        for _ in range(int(10 * stepLength)):
            # Compute forward and backward differences in the x and y directions
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
        #print("lsf_updated:\n", lsf.round(4))
        
        return struc, lsf






