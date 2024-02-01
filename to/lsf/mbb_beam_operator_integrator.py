import numpy as np


class MbbBeamOperatorIntegrator:
    def __init__(self, nu, E0, nelx, nely, struc):
        """
        初始化 MbbBeamOperatorIntegrator 类

        参数:
        nu: 泊松比
        """
        self.nu = nu # Poisson's ration
        self.E0 = E0 # Young's module
        self.nelx = nelx
        self.nely = nely
        self.struc = struc # 设计变量

    def matrix(self):
        """
        Note:
        考虑四个单元的四边形网格中的 0 号单元：
        0,1 - 6,7
        2,3 - 8,9
        拓扑 : cell2dof : 0,1,6,7,8,9,2,3
        FEALPy : cell2dof : 0,1,2,3,6,7,8,9
        """
        nu = self.nu
        E0 = self.E0
        k = np.array([1/2 - nu/6,   1/8 + nu/8,   -1/4 - nu/12, -1/8 + 3 * nu/8,
                    -1/4 + nu/12,  -1/8 - nu/8,    nu/6,         1/8 - 3 * nu/8])
        KE = E0 / (1 - nu**2) * np.array([
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


    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """
        构建 SIMP 中 MBB 梁有限元矩阵

        参数:
        space (Tuple): 有限元空间
        index (Union[np.s_, np.ndarray]): 选定的单元索引，默认为全部单元
        cellmeasure (Optional[np.ndarray]): 对应单元的度量，默认为 None
        out (Optional[np.ndarray]): 输出矩阵，默认为 None

        返回:
        Optional[np.ndarray]: 如果 out 参数为 None，则返回 SIMP 中 MBB 梁有限元矩阵，否则不返回
        """
        struc = self.struc

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
            KE = self.matrix()

            # 用 FEALPy 中的自由度替换 Top 中的自由度
            idx = np.array([0, 1, 6, 7, 2, 3, 4, 5], dtype=np.int_)
            KE = KE[idx, :][:, idx]

            # 确保矩阵非奇异
            struc = np.maximum(struc, 0.0001).ravel()

            print("K:", K.shape)
            K[:] = np.einsum('i, jk -> ijk', struc, KE)

        if out is None:
            return K
