import numpy as np

class BeamOperatorIntegrator:
    def __init__(self, E, KE):
        """
        初始化 BeamOperatorIntegrator 类

        Parameters:
        - E ( ndarray - (NC, ) ): 单元刚度矩阵中的等效杨氏模量.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 具有 unit 杨氏模量的单元刚度矩阵.
        """
        self.E = E
        self.KE = KE

    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """
        构建 PLSF 中梁的单元刚度矩阵.

        参数:
        - space (tuple): 有限元空间.
        - index (Union[np.s_, ndarray]): 选定的单元索引，默认为全部单元.
        - cellmeasure (Optional[ndarray]): 对应单元的度量，默认为 None.
        - out (Optional[ndarray]): 输出矩阵，默认为 None.

        返回:
        - Optional[ndarray]: 如果 out 参数为 None，则返回梁的单元刚度矩阵，否则不返回.
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
            E = self.E

            # 用 FEALPy 中的自由度替换 Top 中的自由度
            idx = np.array([0, 1, 6, 7, 2, 3, 4, 5], dtype=np.int_)
            KE = KE[idx, :][:, idx]

            K[:] = np.einsum('i, jk -> ijk', E, KE)

        if out is None:
            return K

