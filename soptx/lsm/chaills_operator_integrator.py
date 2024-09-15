import numpy as np

class OperatorIntegrator:
    def __init__(self, struc, KE):
        """
        Parameters:
        - struc ( ndarray - (nely, nelx) ): 表示结构的 solid(1) 和 void(0) 单元.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        """
        self.struc = struc
        self.KE = KE

    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """
        Parameters:
        - space (Tuple): 有限元空间.
        - index (Union[np.s_, ndarray]): 选定的单元索引，默认为全部单元.
        - cellmeasure (Optional[ndarray]): 对应单元的度量，默认为 None.
        - out (Optional[ndarray]): 输出矩阵，默认为 None.

        Returns:
        - Optional[ndarray - (NC, ldof*GD, ldof*GD) ]: 如果 out 参数为 None，则返回梁的单元刚度矩阵，否则不返回.
        """

        mesh = space[0].mesh
        ldof = space[0].number_of_local_dofs()
        GD = mesh.geo_dimension()
        NC = mesh.number_of_cells()

        if cellmeasure is None:
            cellmeasure = mesh.entity_measure('cell', index=index)

        if out is None:
            K = np.zeros((NC, GD*ldof, GD*ldof), dtype=np.float64)
        else:
            assert out.shape == (NC, GD*ldof, GD*ldof)
            K = out

        if space[0].doforder == 'sdofs':
            pass
        elif space[0].doforder == 'vdims':
            KE = self.KE
            struc = self.struc

            # 确保矩阵非奇异
            struc = np.maximum(struc, 0.0001)
            # 在将结构乘上去时，注意单元是按列排序的，所以 struc 要按列展开
            K[:] = np.einsum('i, jk -> ijk', struc.flatten(order='F'), KE)

        if out is None:
            return K
