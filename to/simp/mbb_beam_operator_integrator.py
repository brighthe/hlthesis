import numpy as np
from typing import Optional, Tuple

from fealpy.fem.precomp_data import data

class MbbBeamOperatorIntegrator:
    def __init__(self, nu, nelx, nely, xPhys, penal, E0, Emin,):
        """
        初始化 MbbBeamOperatorIntegrator 类

        参数:
        nu: 泊松比
        """
        self.nu = nu
        self.nelx = nelx
        self.nely = nely
        self.xPhys = xPhys
        self.penal = penal
        self.E0 = E0
        self.Emin = Emin

    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """
        构建 MBB 梁有限元矩阵

        参数:
        space (Tuple): 有限元空间
        index (Union[np.s_, np.ndarray]): 选定的单元索引，默认为全部单元
        cellmeasure (Optional[np.ndarray]): 对应单元的度量，默认为 None
        out (Optional[np.ndarray]): 输出矩阵，默认为 None

        返回:
        Optional[np.ndarray]: 如果 out 参数为 None，则返回 MBB 梁有限元矩阵，否则不返回
        """
        nu = self.nu
        nelx = self.nelx
        nely = self.nely
        xPhys = self.xPhys
        penal = self.penal
        E0 = self.E0
        Emin = self.Emin
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
            A11 = np.array([[12, 3, -6, -3], [3, 12, 3, 0], [-6, 3, 12, -3], [-3, 0, -3, 12]])
            A12 = np.array([[-6, -3, 0, 3], [-3, -6, -3, -6], [0, -3, -6, 3], [3, -6, 3, -6]])
            B11 = np.array([[-4, 3, -2, 9], [3, -4, -9, 4], [-2, -9, -4, -3], [9, 4, -3, -4]])
            B12 = np.array([[2, -3, 4, -9], [-3, 2, 9, -2], [4, 9, 2, 3], [-9, -2, 3, 2]])

            KE = 1/(1-nu**2)/24 * (np.block([[A11, A12], [A12.T, A11]]) +
                nu * np.block([[B11, B12], [B12.T, B11]]))

            xPhys = np.full((nely, nelx), 0.5)

            E = Emin + xPhys.ravel()**penal * (E0 - Emin)

            K[:] = np.einsum('i, jk -> ijk', E, KE)

        if out is None:
            return K








