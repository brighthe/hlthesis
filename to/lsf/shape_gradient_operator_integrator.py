import numpy as np

class BeamOperatorIntegrator:
    def __init__(self, E0, E1, nu, ew, eh, Phi):
        """
        初始化 BeamOperatorIntegrator 类

        Parameters:
        - KE ( ndarray - (ldof*GD, ldof*GD) ): 单元刚度矩阵.
        - Phi
        """
        self.E0 = E0
        self.E1 = E1
        self.nu = nu
        self.ew = ew
        self.eh = eh
        self.Phi = Phi


    def basic_KE(self, E, nu, a, b):
        """
        Returns the element stiffness matrix of a full/empty element.

        采用平面应力假设，根据二维线弹性理论计算局部刚度矩阵.
        
        Parameters:
        - E : Young's modulus, 代表材料的刚性
        - nu: Poisson ratio;
        - a: 有限元单元的几何宽度.
        - b: 有限元单元的几何高度.
        
        Returns:
        - KE : a 8-by-8 stiffness matrix.
        """
        # 刚度矩阵的系数
        k = np.array([-1/(6*a*b) * ( nu*a**2 - 2*b**2 - a**2 ), 1/8*nu + 1/8,
                    -1/(12*a*b) * ( nu*a**2 + 4*b**2 - a**2 ), 3/8*nu - 1/8,
                    1/(12*a*b) * ( nu*a**2 - 2*b**2 - a**2 ), -1/8*nu - 1/8,
                    1/(6*a*b) * ( nu*a**2 + b**2 - a**2 ), -3/8*nu + 1/8])
        
        KE = E / (1 - nu**2) * np.array([
            [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
            [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
            [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
            [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
            [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
            [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
            [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
            [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
        ])

        # 用 FEALPy 中的自由度替换 Top 中的自由度
        idx = np.array([0, 1, 6, 7, 2, 3, 4, 5], dtype=np.int_)
        KE = KE[idx, :][:, idx]


        return KE


    def assembly_cell_matrix(self, space, index=np.s_[:], cellmeasure=None, out=None):
        """
        构建 LSM 中梁的单元刚度矩阵.

        Parameters:
        - space (Tuple): 有限元空间.
        - index (Union[np.s_, ndarray]): 选定的单元索引，默认为全部单元.
        - cellmeasure (Optional[ndarray]): 对应单元的度量，默认为 None.
        - out (Optional[ndarray]): 输出矩阵，默认为 None.

        Returns:
        - Optional[ndarray - (NC, ldof*GD, ldof*GD) ]: 如果 out 参数为 None，则返回梁的单元刚度矩阵，否则不返回.
        """

        Phi = self.Phi
        E1 = self.E1
        E0 = self.E0
        nu = self.nu
        ew = self.ew
        eh = self.eh

        mesh = space[0].mesh
        cell = mesh.entity('cell')
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

            # Element is inside the boundary
            inside = np.min(Phi[cell[:]], axis=1) > 0
            K[inside, :, :] = self.basic_KE(E=E1, nu=nu, a=ew, b=eh)

            # Element is outside the boundary
            outside = np.max(Phi[cell], axis=1) < 0
            K[outside, :, :] = self.basic_KE(E=E0, nu=nu, a=ew, b=eh)

            # Element is cut by the boundary
            s, t = np.meshgrid(np.linspace(-1, 1, 21), np.linspace(-1, 1, 21))
            # 插值水平集函数, 将每个节点的水平集函数值与其相应的形函数值相乘，然后对所有节点求和
            cut = ~(inside | outside)
            c0 = ( (1 - s.flatten('F')) * (1 - t.flatten('F')) )[:, np.newaxis]
            c1 = ( (1 + s.flatten('F')) * (1 - t.flatten('F')) )[:, np.newaxis]
            c2 = ( (1 + s.flatten('F')) * (1 + t.flatten('F')) )[:, np.newaxis]
            c3 = ( (1 - s.flatten('F')) * (1 + t.flatten('F')) )[:, np.newaxis]
            tmpPhi = c0 / 4 * Phi.flatten('F')[cell[cut, 0]] + \
                     c1 / 4 * Phi.flatten('F')[cell[cut, 1]] + \
                     c2 / 4 * Phi.flatten('F')[cell[cut, 2]] + \
                     c3 / 4 * Phi.flatten('F')[cell[cut, 3]]
            print("tmpPhi:", tmpPhi.shape, "\n", tmpPhi.round(4))
            # 计算覆盖材料区域的面积比
            area_ratio = np.sum(tmpPhi >= 0, axis=0) / len(s.flatten('F'))
            print("area_ratio:", area_ratio)
            K[cut, :, :] = area_ratio.reshape(-1, 1, 1) * self.basic_KE(E=E1, nu=nu, a=ew, b=eh)

        if out is None:
            return K
