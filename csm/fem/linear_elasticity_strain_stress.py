import argparse
import os
import numpy as np

from scipy.sparse.linalg import spsolve

from linear_elasticity_model2d import BoxDomainData
from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import VectorSourceIntegrator
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
print("NN:", NN)

output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'TraingleMesh.vtu')


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

integrator3 = VectorSourceIntegrator(f = pde.source, q=5)

lform = LinearForm(vspace)
lform.add_domain_integrator(integrator3)
FK = integrator3.assembly_cell_vector(space = vspace)
lform.assembly()
F = lform.get_vector()

if hasattr(pde, 'dirichlet'):
    bc = DirichletBC(space=vspace, gD=pde.dirichlet, threshold=pde.is_dirichlet_boundary)
    K, F = bc.apply(K, F, uh)

uh.flat[:] = spsolve(K, F)
print("uh:", uh.shape, "\n", uh)

reshaped_uh = uh.reshape(-1)
cell2dof = space.cell_to_dof()
updated_cell2dof = np.repeat(cell2dof * GD, GD, axis=1) + np.tile(np.array([0, 1]), (NC, ldof))
cell_displacements = reshaped_uh[updated_cell2dof] # (NC, ldof*GD)
print("cell_displacements:", cell_displacements.shape, "\n", cell_displacements)

def compute_strain(uh):
    """
    @brief 给定一个位移，计算相应的应变，这里假设是线性元

    Parameters:
    uh: ndarray
        位移场，形状为 (NN, GD)，其中 NN 是节点数，GD 是空间维度

    Returns:
    s: ndarray
        计算得到的应变张量，形状为 (NC, GD, GD)，其中 NC 是单元数，GD 是空间维度
    """
    cell2dof = space.cell_to_dof()
    NC = mesh.number_of_cells()
    #gphi = mesh.grad_lambda()  # NC x 3 x 2
    #print("gphi:", gphi.shape, "\n", gphi)

    ldof = (p+1)*(p+2)//2
    idx = np.arange(0, ldof)
    idx0 = np.floor((-1 + np.sqrt(1 + 8*idx))/2)
    multiIndex = np.zeros((ldof, 3), dtype=np.int_)
    multiIndex[:,2] = idx - idx0*(idx0 + 1)/2
    multiIndex[:,1] = idx0 - multiIndex[:,2]
    multiIndex[:,0] = p - multiIndex[:, 1] - multiIndex[:, 2]
    bcs = multiIndex/p
    grad = space.grad_basis(bcs) # (ldof, NC, ldof, GD)
    grad = np.diagonal(grad, axis1=0, axis2=2)  # 形状变为 (NC, GD, ldof)
    gphi = np.swapaxes(grad, 1, 2)

    s = np.zeros((NC, GD, GD), dtype=np.float64)
    if uh.space.doforder == 'sdofs':
        uh = uh.T
    for i in range(GD):
        for j in range(i, GD):
            if i ==j:
                s[:, i, i] = np.sum(uh[:, i][cell2dof] * gphi[:, :, i], axis=-1)
            else:
                val = np.sum(uh[:, i][cell2dof] * gphi[:, :, j], axis=-1)
                val += np.sum(uh[:, j][cell2dof] * gphi[:, :, i], axis=-1)
                val /= 2.0
                s[:, i, j] = val
                s[:, j, i] = val
    return s

strain = compute_strain(uh=uh)
print("strain:", strain.shape, "\n", strain)

def compute_stress(mu, lam, strain):
    """
    参数:
    strain (np.ndarray): 应变矩阵，形状为 (NC, GD, GD)

    返回:
    np.ndarray: 应力矩阵
    """
    if GD == 2:
        D = np.array([
            [2 * mu + lam, 0, 0, lam],
            [0, mu, mu, 0],
            [0, mu, mu, 0],
            [lam, 0, 0, 2 * mu + lam]
        ]).reshape(GD, GD, GD, GD)

        stress = np.einsum('ijkl, cjl -> cik', D, strain)

    elif GD == 3:
        D = np.array([
            [2 * mu + lam, 0, 0, 0, lam, 0, 0, 0, lam],
            [0, mu, 0, mu, 0, 0, 0, 0, 0],
            [0, 0, mu, 0, 0, 0, mu, 0, 0],
            [0, mu, 0, mu, 0, 0, 0, 0, 0],
            [lam, 0, 0, 0, 2 * mu + lam, 0, 0, 0, lam],
            [0, 0, 0, 0, 0, mu, 0, mu, 0],
            [0, 0, mu, 0, 0, 0, mu, 0, 0],
            [0, 0, 0, 0, 0, mu, 0, mu, 0],
            [lam, 0, 0, 0, lam, 0, 0, 0, 2 * mu + lam],
        ]).reshape(GD, GD, GD, GD)

        stress = np.einsum('ijkl, cjl -> cik', D, strain)

    return stress

stress = compute_stress(mu=mu, lam=lambda_, strain=strain)
print("stress:", stress.shape, "\n", stress)

mesh.nodedata['u'] = uh[:, 0]
mesh.nodedata['v'] = uh[:, 1]

import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_cell(axes, showindex=True, color='b', marker='s', markersize=5, fontsize=15, fontcolor='b')
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=5, fontsize=15, fontcolor='r')
plt.show()

mesh.to_vtk(fname=fname)
