import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsf:

    def __init__(self, nelx: int = 60, nely: int = 30, volReq: float = 0.3, 
                stepLength: int = 3, numReinit: int = 2, topWeight: int = 2):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 30.
        - volReq : 最终设计所需的体积分数. Defaults to 0.3.
        - stepLength (int): The CFL time step size used in each iteration of the evolution equation. Defaults to 3.
        - numReinit (int): The number of algorithm iterations before performing level set reinitialization. Defaults to 2.
        - topWeight (int): Weight of the topology derivative term in the evolution equation. Defaults to 2.
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
        #self._mesh = QuadrangleMesh.from_box(box = [0, self._nelx+2, 0, self._nely+2], 
        #nx = self._nelx + 2, ny = self._nely + 2)

    def reinit(self, strucFull):
        """
        根据给定的结构重置化水平集函数

        该函数通过添加 void 单元的边界来扩展输入结构，计算到最近的 solid 和 void 单元
        的欧几里得距离，并计算水平集函数，该函数在 solid phase 内为负，在 void phase 中为正.

        Parameters:
        - strucFull (numpy.ndarray): A 2D array representing the solid (1) and void (0) cells of the structure.

        Returns:
        numpy.ndarray: A 2D array of the same shape as 'struc', 表示重置化后的水平集函数
        """
        from scipy import ndimage

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
        nelx, nely = struc.shape
        integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        bform.assembly()
        K = bform.get_matrix()
        #print("K:", K.shape, "\n", K.toarray().round(4))

        # 定义荷载
        F = np.zeros(vgdof)
        F[vgdof-1] = 1
        #print("F:", F.shape, "\n", F.round(4))

        # 定义支撑(边界处理)
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
        updated_cell2dof = np.repeat(cell2dof * GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
        idx = np.array([0, 1, 4, 5, 6, 7, 2, 3], dtype=np.int_)
        # 用 Top 中的自由度替换 FEALPy 中的自由度
        updated_cell2dof = updated_cell2dof[:, idx]
        ue = reshaped_uh[updated_cell2dof] # (NC, ldof*GD)


        return uh, ue




