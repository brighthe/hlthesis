import numpy as np

from fealpy.mesh import QuadrangleMesh

class TopLsfShapeGrad:

    def __init__(self, domain_width, domain_hight, nelx, nely, volReq, stepLength, numReinit):
        '''
        初始化拓扑优化问题 - 形状梯度的水平集方法

        Parameters: 
        - domain: 设计区域.
        - nelx (int): 沿设计区域水平方向的单元数.
        - nely (int): 沿设计区域垂直方向的单元数.
        - volReq (float) : 最终设计所需的体积分数. 0.2 < volReq < 0.7.
        - stepLength (int): 每次迭代中求解演化方程的 CFL 时间步长.
                            min(nelx,nely)/10 < stepLength < max(nelx,nely)/5.
        - numReinit (int): 水平集函数重置化为符号距离函数的频率. 2 < numReinit < 6.

        Note:
            numReinit 和 topWeight 会影响设计中形成孔洞的能力.
        '''

        self.domain_width = domain_width
        self.domain_hight = domain_hight
        self.nelx = nelx
        self.nely = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit

    def generate_mesh(self, domain, nelx, nely):

        mesh = QuadrangleMesh.from_box(box = domain, nx = nelx, ny = nely)

        return mesh

    def plot(self, x, y, z, label):

        import matplotlib.pyplot as plt

        plt.scatter(x, y, c=z, cmap='viridis', label=label, edgecolors='k', s=50)

        plt.colorbar(label='Phi values')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot_mesh(self, x0, y0, label0, x, y, z, label):

        import matplotlib.pyplot as plt

        plt.scatter(x, y, c=z, cmap='cool', label=label, edgecolors='r', marker='s', s=50)

        plt.scatter(x0, y0, color='green', marker='o', label=label0)

        plt.colorbar(label='Phi values')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.legend()
        plt.grid(True)
        plt.show()


    def init_lsf(self, mesh):
        '''
        初始化设计区域内的水平集函数，以分布初始孔洞的符号距离函数表示.

        该函数通过定义一系列圆形孔洞，并计算网格点到这些孔洞边界的最短距离，来初始化水平集函数
        水平集函数在孔洞内部为负，在孔洞外部为正

        Parameters:
        - mesh (object): 初始的水平集网格.

        Returns:
        - ls_Phi ( ndarray - ( (domain_width+1)*(domain_hight+1), ) ): 初始的水平集函数.

        '''
        domain_width = self.domain_width
        domain_hight = self.domain_hight
        ls_node = mesh.entity('node')
        ls_x = ls_node[:, 0]
        ls_y = ls_node[:, 1]
        # 定义初始孔洞的圆心
        cx = domain_width/200 * np.array([33.33,  100,  166.67,   0,    66.67,  133.33,
                                           200,  33.33,  100,   166.67,   0,    66.67,
                                         133.33,  200,  33.33,   100,   166.67], dtype=np.float64)
        cy = domain_hight/100 * np.array([0,  0,  0,   25,  25, 25,
                                          25, 50, 50,  50,  75, 75,
                                          75, 75, 100, 100, 100], dtype=np.float64)
        # 定义初始孔洞的半径
        cr = domain_hight / 10

        # 计算每个网格点与所有圆心的欧式距离
        tmpPhi = -np.sqrt((ls_x[:, None] - cx)**2 + (ls_y[:, None] - cy)**2) + cr

        # 对于每个节点，取其到所有圆心距离的最大负值，以确保节点若在任一孔洞内部，则其水平集函数值为负
        ls_Phi = -np.max(tmpPhi, axis=1)

        return ls_Phi


    def fe_analysis(self, mesh, E0, E1, nu, ew, eh, Phi, F, fixeddofs):
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
        from shape_gradient_operator_integrator import BeamOperatorIntegrator
        from fealpy.fem import BilinearForm
        from fealpy.functionspace import LagrangeFESpace as Space
        from scipy.sparse import spdiags

        p = 1
        space = Space(mesh, p=p, doforder='vdims')
        GD = 2
        uh = space.function(dim=GD)
        nLoads = F.shape[-1]
        uh = np.repeat(uh[:, :, np.newaxis], nLoads, axis=2)
        vspace = GD*(space, )
        e0 = E0
        e1 = E1
        phi = Phi
        integrator = BeamOperatorIntegrator(E0=e0, E1=e1, nu=nu, ew=ew, eh=eh, Phi=phi)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator)
        KK = integrator.assembly_cell_matrix(space=vspace)
        #print("KK:", KK.shape, "\n", KK)
        K = bform.assembly()
        print("K0:", K.shape, "\n", K.toarray())

        dflag = fixeddofs
        F = F - K@uh.reshape(-1, nLoads)
        bdIdx = np.zeros(K.shape[0], dtype=np.int_)
        bdIdx[dflag.flat] = 1
        D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
        D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
        K = D0@K@D0 + D1
        print("K1:", K.shape, "\n", K.toarray())
        F[dflag.flat] = uh.reshape(-1, nLoads)[dflag.flat]

        ## 线性方程组求解
        #uh.flat[:] = spsolve(K, F)








