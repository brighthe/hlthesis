import numpy as np


class BeamOperatorIntegrator:
    def __init__(self, x, penal, KE):
        """
        初始化 BeamOperatorIntegrator 类

        Parameters:
        - x (ndarray - (nely, nelx) ): 设计变量.
        - penal (float): penalization power.
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        """
        self.x = x
        self.penal = penal
        self.KE = KE


    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """

        Parameters:
        - space (Tuple): 有限元空间.
        - index (Union[np.s_, ndarray]): 选定的单元索引，默认为全部单元.
        - cellmeasure (Optional[ndarray]): 对应单元的度量，默认为 None.
        - out (Optional[ndarray]): 输出矩阵，默认为 None.

        Returns:
        - Optional[ndarray]: 如果 out 参数为 None，则返回 SIMP 中梁有限元矩阵，否则不返回.
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
            x = self.x
            penal = self.penal
            KE = self.KE

            # 用 FEALPy 中的自由度替换 Top 中的自由度
            idx = np.array([0, 1, 6, 7, 2, 3, 4, 5], dtype=np.int_)
            KE = KE[idx, :][:, idx]

            coef = x ** penal
            # 在将结构乘上去时，注意单元是按列排序的，所以 coef 要按列展开
            K[:] = np.einsum('i, jk -> ijk', coef.flatten(order='F'), KE)

        if out is None:
            return K
