import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

from scipy.sparse.linalg import spsolve

from fealpy.pde.linear_elasticity_model import BoxDomainData
from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import VectorSourceIntegrator
from fealpy.fem import VectorMassIntegrator
from fealpy.fem import BilinearForm
from fealpy.fem import LinearForm
from fealpy.fem import DirichletBC


## 参数解析
parser = argparse.ArgumentParser(description=
        """
        单纯形网格（三角形、四面体）网格上任意次有限元方法
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Lagrange 有限元空间的次数, 默认为 1 次.')

parser.add_argument('--GD',
        default=2, type=int,
        help='模型问题的维数, 默认求解 2 维问题.')

parser.add_argument('--nrefine',
        default=2, type=int,
        help='初始网格加密的次数, 默认初始加密 2 次.')

parser.add_argument('--scale',
        default=1, type=float,
        help='网格变形系数，默认为 1')

parser.add_argument('--doforder',
        default='vdims', type=str,
        help='自由度排序的约定，默认为 vdims')

args = parser.parse_args()
p = args.degree
GD = args.GD
n = args.nrefine
scale = args.scale
doforder = args.doforder

pde = BoxDomainData()

mu = pde.mu
lambda_ = pde.lam
domain = pde.domain()
mesh = pde.triangle_mesh()
#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_cell(axes, showindex=True, color='k', marker='s', markersize=2, fontsize=8, fontcolor='k')
#mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=2, fontsize=8, fontcolor='r')
#plt.show()
NN = mesh.number_of_nodes()
print("NN:", NN)
NC = mesh.number_of_cells()
print("NC:", NC)
node = mesh.entity('node')
cell = mesh.entity('cell')
print("cell:\n", cell.shape, "\n", cell)

#output = './mesh/'
#if not os.path.exists(output):
#    os.makedirs(output)
#fname = os.path.join(output, 'TriangleMesh.vtu')

space = Space(mesh, p=p, doforder=doforder)
uh = space.function(dim=GD)
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
ldof = vspace[0].number_of_local_dofs()
vldof = ldof * GD
print("vgdof", vgdof)
print("vldof", vldof)

from fealpy.fem import ProvidesSymmetricTangentOperatorIntegrator
def linear_tangent_matrix(lam, mu, GD):
    #lam = self.model.lam # 拉梅第一参数
    #mu = self.model.mu # 拉梅第二参数
    #mesh = self.mesh
    #GD = mesh.geo_dimension()
    n = GD*(GD+1)//2
    D = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        D[i, i] += mu
    for i in range(GD):
        for j in range(i, GD):
            if i == j:
                D[i, i] += mu+lam
            else:
                D[i, j] += lam
                D[j, i] += lam
    return D

def dsigma_depsilon(c0, D0):
    """
    @brief 计算应力关于应变的导数矩阵
    @param phi 单元重心处的相场函数值, (NC, )
    @param uh 位移
    @return D 单元刚度系数矩阵
    """
    #eps = 1e-10
    #lam = self.model.lam # 拉梅第一参数
    #mu = self.model.mu # 拉梅第二参数

    #qf =  self.mesh.integrator(1, 'cell')
    #bc, ws = qf.get_quadrature_points_and_weights()

#   #     bc = np.array([1 / 3, 1 / 3, 1 / 3], dtype=np.float64)
    #c0 = (1 - phi(bc[-1])) ** 2 + eps
    D = np.einsum('i, jk -> ijk', c0, D0)
#        c0 = (1 - phi(bc)) ** 2 + eps
#        D0 = np.array([[2*mu+lam, lam, 0], [lam, 2*mu+lam, 0], [0, 0, mu]],
#                dtype=np.float_)
#        D = np.einsum('i, jk -> ijk', c0, D0)
    return D

for i in range(2):
    integrator1 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+1)
    bform1 = BilinearForm(vspace)
    bform1.add_domain_integrator(integrator1)
    KK1 = integrator1.assembly_cell_matrix(space=vspace)
    print("KK1", KK1.shape, "\n", KK1[0])
    bform1.assembly()
    K1 = bform1.get_matrix()
    print("K1:", K1.shape, "\n", K1.toarray().round(4))

    coef = np.full(NC, 0.5)
    integrator2 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+1, c=coef)
    bform2 = BilinearForm(vspace)
    bform2.add_domain_integrator(integrator2)
    KK2 = integrator2.assembly_cell_matrix(space=vspace)
    print("KK2", KK2.shape, "\n", KK2[0])
    bform2.assembly()
    K2 = bform2.get_matrix()
    print("K2:", K2.shape, "\n", K2.toarray().round(4))

    D0 = linear_tangent_matrix(lam=lambda_, mu=mu, GD=GD)
    D = dsigma_depsilon(c0=coef, D0=D0)
    bform3 = BilinearForm(vspace)
    integrator3 = ProvidesSymmetricTangentOperatorIntegrator(D)
    bform3.add_domain_integrator(integrator3)
    KK3 = integrator3.assembly_cell_matrix(space=vspace)
    print("KK3", KK3.shape, "\n", KK3[0])
    bform3.assembly()
    K3 = bform3.get_matrix()
    print("K3:", K3.shape, "\n", K3.toarray().round(4))

    KK13_test = KK1 - KK3
    print("KK13_test:\n", KK13_test)
    KK23_test = KK2 - KK3
    print("KK23_test:\n", KK23_test)

    K12_test = K1 - K2
    print("K12_test:\n", K12_test.toarray().round(4))
    K13_test = K1 - K3
    print("K13_test:\n", K13_test.toarray().round(4))
    K23_test = K2 - K3
    print("K23_test:\n", K23_test.toarray().round(4))





