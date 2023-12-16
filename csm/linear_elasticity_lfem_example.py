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
mesh = pde.delaunay_mesh()
#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_cell(axes, showindex=True, color='k', marker='s', markersize=2, fontsize=8, fontcolor='k')
#mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=2, fontsize=8, fontcolor='r')
#plt.show()
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
node = mesh.entity('node')
cell = mesh.entity('cell')

output = './csm_results/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'DelaunayMesh.vtu')

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
print("KK", KK.shape)
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))

integrator2 = VectorMassIntegrator(c=1, q=5)

bform2 = BilinearForm(vspace)
bform2.add_domain_integrator(integrator2)
MK = integrator2.assembly_cell_matrix(space=vspace)
print("MK:", MK.shape)
bform2.assembly()
M = bform2.get_matrix()
print("M:", M.shape)

integrator3 = VectorSourceIntegrator(f = pde.source, q=5)

lform = LinearForm(vspace)
lform.add_domain_integrator(integrator3)
FK = integrator3.assembly_cell_vector(space = vspace)
print("FK[0]:", FK.shape)
lform.assembly()
F = lform.get_vector()
print("F:", F.shape, "\n", F.round(4))

ipoints = space.interpolation_points()
fh = pde.source(p=ipoints)
fh_1 = np.zeros(M.shape[0])
fh_1[::GD] = fh[:,0]
fh_1[1::GD] = fh[:,1]
Fh = M @ fh_1
print("Fh:", Fh.shape, "\n", Fh.round(4))

print("error:", np.sum(np.abs(F - Fh)))

if hasattr(pde, 'dirichlet'):
    bc = DirichletBC(space=vspace, gD=pde.dirichlet, threshold=pde.is_dirichlet_boundary)
    K, Fh = bc.apply(K, Fh, uh)

print("K:", K.shape, "\n", K.toarray().round(4))
print("Fh:", Fh.shape, "\n", Fh.round(4))

uh.flat[:] = spsolve(K, Fh)
print("uh:", uh.shape, uh)
mesh.nodedata['u'] = uh[:, 0]
mesh.nodedata['v'] = uh[:, 1]

mesh.to_vtk(fname=fname)
