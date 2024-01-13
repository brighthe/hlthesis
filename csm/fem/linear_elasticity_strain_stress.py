import argparse
import os
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
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
node = mesh.entity('node')
cell = mesh.entity('cell')

output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'TraingleMesh.vtu')

import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_cell(axes, showindex=True, color='b', marker='s', markersize=5, fontsize=15, fontcolor='b')
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=5, fontsize=15, fontcolor='r')
plt.show()

space = Space(mesh, p=p, doforder=doforder)
uh = space.function(dim=GD)
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
ldof = vspace[0].number_of_local_dofs()
vldof = ldof * GD
print("vgdof", vgdof)
print("vldof", vldof)

integrator1 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+1)

bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator1)
KK = integrator1.assembly_cell_matrix(space=vspace)
bform.assembly()
K = bform.get_matrix()

integrator2 = VectorMassIntegrator(c=1, q=5)

bform2 = BilinearForm(vspace)
bform2.add_domain_integrator(integrator2)
MK = integrator2.assembly_cell_matrix(space=vspace)
bform2.assembly()
M = bform2.get_matrix()

integrator3 = VectorSourceIntegrator(f = pde.source, q=5)

lform = LinearForm(vspace)
lform.add_domain_integrator(integrator3)
FK = integrator3.assembly_cell_vector(space = vspace)
lform.assembly()
F = lform.get_vector()

ipoints = space.interpolation_points()
fh = pde.source(p=ipoints)
fh_1 = np.zeros(M.shape[0])
fh_1[::GD] = fh[:,0]
fh_1[1::GD] = fh[:,1]
Fh = M @ fh_1

if hasattr(pde, 'dirichlet'):
    bc = DirichletBC(space=vspace, gD=pde.dirichlet, threshold=pde.is_dirichlet_boundary)
    K, Fh = bc.apply(K, Fh, uh)

uh.flat[:] = spsolve(K, Fh)
print("uh:", uh.shape, uh)

reshaped_uh = uh.reshape(-1)
cell2dof = space.cell_to_dof()
print("cell2dof:", cell2dof.shape, "\n", cell2dof)
updated_cell2dof = np.repeat(cell2dof * GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
cell_displacements = reshaped_uh[updated_cell2dof] # (NC, ldof*GD)
print("cell_displacements:", cell_displacements.shape, "\n", cell_displacements)

ldof = (p+1)*(p+2)//2
idx = np.arange(0, ldof)
idx0 = np.floor((-1 + np.sqrt(1 + 8*idx))/2)
multiIndex = np.zeros((ldof, 3), dtype=np.int_)
multiIndex[:,2] = idx - idx0*(idx0 + 1)/2
multiIndex[:,1] = idx0 - multiIndex[:,2]
multiIndex[:,0] = p - multiIndex[:, 1] - multiIndex[:, 2]
bcs = multiIndex/p
grad = space.grad_basis(bcs) # (ldof, NC, ldof, GD)
grad_1 = np.swapaxes(grad, 0, 1)  # 将形状变为 (328, 3, 3, 2)
# 提取对角线元素
grad_1 = np.diagonal(grad_1, axis1=1, axis2=2)  # 形状变为 (328, 3, 2)
grad_1 = np.swapaxes(grad_1, 1, 2)
print("grad_1:", grad_1.shape, "\n", grad_1)

grad_2 = mesh.grad_lambda()
print("grad_2:", grad_2.shape, "\n", grad_2)

def compute_strain(grad, cell_displacement):
    """
    参数:
    grad (np.ndarray): 形状函数梯度，形状为 (NC, ldof, GD)
    cell_displacement (np.ndarray): 单元位移，形状为 (NC, ldof*GD)
    返回:
    np.ndarray: 应变矩阵
    """
    NC, ldof, GD = grad.shape
    strain_dim = GD * (GD + 1) // 2
    strain = np.zeros((NC, strain_dim))

    # (NC, ldof*GD) -> (NC, ldof, GD)
    cell_displacement_reshaped = cell_displacement.reshape(NC, ldof, GD)

    # 对于二维和三维问题，使用不同的索引策略
    if GD == 2:
        # 二维应变计算 (ε_xx, ε_yy, γ_xy)
        strain[:, 0] = np.einsum('ni, ni -> n', grad[:, :, 0], cell_displacement_reshaped[:, :, 0])  # ε_xx
        strain[:, 1] = np.einsum('ni, ni -> n', grad[:, :, 1], cell_displacement_reshaped[:, :, 1])  # ε_yy
        strain[:, 2] = np.einsum('ni, ni -> n', grad[:, :, 0], cell_displacement_reshaped[:, :, 1]) + \
                       np.einsum('ni, ni -> n', grad[:, :, 1], cell_displacement_reshaped[:, :, 0])  # γ_xy
    elif GD == 3:
        # 三维应变计算 (ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz)
        for i in range(GD):
            strain[:, i] = np.einsum('ni, ni -> n', grad[:, :, i], cell_displacement_reshaped[:, :, i])  # ε_xx, ε_yy, ε_zz
        strain[:, 3] = np.einsum('ni, ni -> n', grad[:, :, 0], cell_displacement_reshaped[:, :, 1]) + \
                       np.einsum('ni, ni -> n', grad[:, :, 1], cell_displacement_reshaped[:, :, 0])  # γ_xy
        strain[:, 4] = np.einsum('ni, ni -> n', grad[:, :, 1], cell_displacement_reshaped[:, :, 2]) + \
                       np.einsum('ni, ni -> n', grad[:, :, 2], cell_displacement_reshaped[:, :, 1])  # γ_yz
        strain[:, 5] = np.einsum('ni, ni -> n', grad[:, :, 0], cell_displacement_reshaped[:, :, 2]) + \
                       np.einsum('ni, ni -> n', grad[:, :, 2], cell_displacement_reshaped[:, :, 0])  # γ_xz

    return strain

strain = compute_strain(grad=grad_1, cell_displacement=cell_displacements)
print("strain:", strain.shape, "\n", strain)

def strain_2(uh):
    """
    @brief 给定一个位移，计算相应的应变，这里假设是线性元
    """
    cell = mesh.entity('cell')
    NC = mesh.number_of_cells()
    gphi = mesh.grad_lambda()  # NC x 3 x 2

    s = np.zeros((NC, GD, GD), dtype=np.float64)
    if uh.space.doforder == 'sdofs':
        uh = uh.T
    for i in range(GD):
        for j in range(i, GD):
            if i ==j:
                s[:, i, i] = np.sum(uh[:, i][cell] * gphi[:, :, i], axis=-1)
            else:
                val = np.sum(uh[:, i][cell] * gphi[:, :, j], axis=-1)
                val += np.sum(uh[:, j][cell] * gphi[:, :, i], axis=-1)
                val /= 2.0
                s[:, i, j] = val
                s[:, j, i] = val
    return s

strain_1 = strain_2(uh=uh)
print("strain_1:", strain_1.shape, "\n", strain_1)

def compute_stress(GD, mu, lam, strain):
    """
    参数:
    strain (np.ndarray): 应变矩阵，形状为 (NC, strain_dim)
    返回:
    np.ndarray: 应力矩阵
    """
    if GD == 2:
        # 二维应力计算
        D = np.array([
            [2*mu+lam, lam, 0],
            [lam, 2*mu+lam, 0],
            [0, 0, mu]
        ])
    elif GD == 3:
        # 三维应力计算
        D = np.array([
            [2*mu+lam, lam, lam, 0, 0, 0],
            [lam, 2*mu+lam, lam, 0, 0, 0],
            [lam, lam, 2*mu+lam, 0, 0, 0],
            [0, 0, 0, mu, 0, 0],
            [0, 0, 0, 0, mu, 0],
            [0, 0, 0, 0, 0, mu]
        ])
    stress = np.einsum('nc, cd -> nd', strain, D)

    return stress

stress = compute_stress(GD=GD, mu=mu, lam=lambda_, strain=strain)
print("stress:", stress.shape, "\n", stress)

mesh.nodedata['u'] = uh[:, 0]
mesh.nodedata['v'] = uh[:, 1]

mesh.to_vtk(fname=fname)
