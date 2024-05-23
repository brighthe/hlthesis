import argparse
import numpy as np

from linear_elasticity_model2d import BoxDomainData2d
from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.fem import LinearElasticityOperatorIntegrator, VectorDiffusionIntegrator
from fealpy.fem import VectorMassIntegrator
from fealpy.fem import BilinearForm


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
        default=1, type=int,
        help='初始网格加密的次数, 默认初始加密 1 次.')

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

pde = BoxDomainData2d()

mu = pde.mu
lambda_ = pde.lam
domain = pde.domain()
mesh = pde.triangle_mesh()
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
print("NC:", NC)
node = mesh.entity('node')
cell = mesh.entity('cell')

space = Space(mesh, p=p, doforder=doforder)
uh = space.function(dim=GD)
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
ldof = vspace[0].number_of_local_dofs()
vldof = ldof * GD
print("vgdof", vgdof)
print("vldof", vldof)

integrator1 = VectorMassIntegrator(c=1, q=p+2)

bform1 = BilinearForm(vspace)
bform1.add_domain_integrator(integrator1)
MK = integrator1.assembly_cell_matrix(space=vspace)
#print("MK:", MK.shape, "\n", MK)
bform1.assembly()
M = bform1.get_matrix()
#print("M:", M.shape, "\n", M.toarray().round(4))

bform2 = BilinearForm(vspace)
bform2.add_domain_integrator(integrator1)
MK_fast = integrator1.assembly_cell_matrix_fast(space=vspace)
#print("MK_fast:", MK_fast.shape, "\n", MK_fast)
bform2.fast_assembly()
M_fast = bform1.get_matrix()
#print("M_fast:", M_fast.shape, "\n", M_fast.toarray().round(4))

MK_equal = MK - MK_fast
#print("MK_equal:\n", MK_equal.round(4))
M_equal = M - M_fast
#print("M_equal:\n", M_equal.toarray().round(4))
MK_equal_test = np.allclose(MK, MK_fast, rtol=1e-05, atol=1e-08)
print("MK_equal_test:\n", MK_equal_test)
M_equal_test = np.allclose(M_fast.toarray(), M_fast.toarray(), rtol=1e-05, atol=1e-08)
print("M_equal_test:\n", M_equal_test)


integrator2 = VectorDiffusionIntegrator(q=p+2)

bform3 = BilinearForm(vspace)
bform3.add_domain_integrator(integrator2)
DK = integrator2.assembly_cell_matrix(space=vspace)
#print("DK:", DK.shape, "\n", DK)
bform3.assembly()
D = bform3.get_matrix()
#print("M:", M.shape, "\n", M.toarray().round(4))

bform4 = BilinearForm(vspace)
bform4.add_domain_integrator(integrator1)
DK_fast = integrator2.assembly_cell_matrix_fast(space=vspace)
#print("DK_fast:", DK_fast.shape, "\n", DK_fast[0])
bform4.fast_assembly()
D_fast = bform1.get_matrix()
#print("D_fast:", D_fast.shape, "\n", D_fast.toarray().round(4))

DK_equal = DK - DK_fast
#print("MK_equal:\n", MK_equal.round(4))
D_equal = D - D_fast
#print("M_equal:\n", M_equal.toarray().round(4))
DK_equal_test = np.allclose(DK, DK_fast, rtol=1e-05, atol=1e-08)
print("DK_equal_test:\n", DK_equal_test)
D_equal_test = np.allclose(D_fast.toarray(), D_fast.toarray(), rtol=1e-05, atol=1e-08)
print("D_equal_test:\n", D_equal_test)


from fealpy.decorator import cartesian, barycentric
@cartesian
def func_coef(p):
    x = p[..., 0]
    y = p[..., 1]
    return x + y

vector_coef = np.full(NC, 2)

integrator3 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+2, c=func_coef)

bform5 = BilinearForm(vspace)
bform5.add_domain_integrator(integrator3)
LK = integrator3.assembly_cell_matrix(space=vspace)
#print("LK", LK.shape, "\n", LK[0])
bform5.assembly()
L = bform5.get_matrix()
#print("L:", L.shape, "\n", L.toarray().round(4))

bform6 = BilinearForm(vspace)
bform6.add_domain_integrator(integrator3)
LK_fast = integrator3.assembly_cell_matrix_fast(space=vspace)
#print("LK", LK.shape)
bform6.fast_assembly()
L_fast = bform5.get_matrix()
#print("L_fast:", L_fast.shape, "\n", L_fast.toarray().round(4))

LK_equal = LK - LK_fast
#print("MK_equal:\n", MK_equal.round(4))
L_equal = L - L_fast
#print("M_equal:\n", M_equal.toarray().round(4))
LK_equal_test = np.allclose(LK, LK_fast, rtol=1e-05, atol=1e-08)
print("LK_equal_test:\n", LK_equal_test)
L_equal_test = np.allclose(L_fast.toarray(), L_fast.toarray(), rtol=1e-05, atol=1e-08)
print("L_equal_test:\n", L_equal_test)
