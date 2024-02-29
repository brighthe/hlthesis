import numpy as np

class BeamOperatorIntegrator:
    def __init__(self, nu, E0):
        """
        初始化 BeamOperatorIntegrator 类

        参数:
        nu: 泊松比
        """
        self.nu = nu # Poisson's ration
        self.E0 = E0 # Young's module

    def stiff_matrix(self):
        """
        Note:
            考虑四个单元的四边形网格中的 0 号单元：
            0,1 - 6,7
            2,3 - 8,9
            拓扑 : cell2dof : 0,1,6,7,8,9,2,3
            FEALPy : cell2dof : 0,1,2,3,6,7,8,9
        """
        nu = self.nu

        A11 = np.array([[12,  3, -6, -3],
                        [ 3, 12,  3,  0],
                        [-6,  3, 12, -3],
                        [-3,  0, -3, 12]])
        A12 = np.array([[-6, -3,  0,  3],
                        [-3, -6, -3, -6],
                        [ 0, -3, -6,  3],
                        [ 3, -6,  3, -6]])
        B11 = np.array([[-4,  3, -2,  9],
                        [ 3, -4, -9,  4],
                        [-2, -9, -4, -3],
                        [ 9,  4, -3, -4]])
        B12 = np.array([[ 2, -3,  4, -9],
                        [-3,  2,  9, -2],
                        [ 4,  9,  2,  3],
                        [-9, -2,  3,  2]])

        KE = 1 / (1-nu**2) / 24 * (np.block([[A11, A12], [A12.T, A11]]) +\
                              nu * np.block([[B11, B12], [B12.T, B11]]))

        return  KE

    def stiff_matrix_2(self):
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

