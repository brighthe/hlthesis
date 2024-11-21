import argparse
import os

import sympy as sp
import numpy as np

from fealpy.mesh import TriangleMesh

from fealpy.fem import BilinearForm
from fealpy.fem import LinearForm
from fealpy.fem import DirichletBC
from fealpy.fem import VectorSourceIntegrator

from HuZhangFiniteElementSpace2d import HuZhangFiniteElementSpace2d
from linear_elasticity_model2d import BoxDomainData2d
from hu_zhang_mass_operator_integrator import HuZhangMassOperatorIntegrator

from scipy.sparse import bmat, coo_matrix, spdiags
from scipy.sparse.linalg import spsolve


from fealpy.functionspace.lagrange_fe_space import LagrangeFESpace

from fealpy.pde.linear_elasticity_model2D import GenLinearElasticitymodel2D

parser = argparse.ArgumentParser(description=
        """
        三角形网格上用胡张元求解线弹性力学问题
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='间断 Lagrange 有限元空间的次数, 默认为 1 次.')

parser.add_argument('--GD',
        default=2, type=int,
        help='模型问题的维数, 默认求解 2 维问题.')

parser.add_argument('--nrefine',
        default=1, type=int,
        help='初始网格加密的次数, 默认初始加密 1 次.')

parser.add_argument('--maxit',
        default=1, type=int,
        help='默认网格加密求解的次数, 默认加密求解 1 次')

parser.add_argument('--bdtype',
        default='displacement', type=str,
        help='边界条件, 默认为位移边界')

args = parser.parse_args()

p = args.degree
GD = args.GD
nrefine = args.nrefine
maxit = args.maxit
bdtype = args.bdtype

pde = BoxDomainData2d()
mesh = pde.triangle_mesh(nx=2, ny=2)
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
node = mesh.entity('node')
cell = mesh.entity('cell')
print("NN:", NN)
print("NC:", NC)

output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'TraingleMesh.vtu')

spacetype = 'D'
doforder = 'sdofs'
space = LagrangeFESpace(mesh, p=p, spacetype=spacetype, doforder=doforder)
uh = space.function(dim=GD)
print("uh:", uh.shape)
vspace = GD*(space, )
vldof = vspace[0].number_of_local_dofs()
print("vldof:", vldof)
vgdof = vspace[0].number_of_global_dofs()
print("vgdof:", vgdof)

tspace = HuZhangFiniteElementSpace2d(mesh, p=p+1)
sigmah = tspace.function(dim=GD)
print("sigmah:", sigmah.shape)
tldof = tspace.number_of_local_dofs()
print("tldof:", tldof)
tgdof = tspace.number_of_global_dofs()
print("tgdof:", tgdof)

integrator1 = HuZhangMassOperatorIntegrator(mesh, p=p+1)
bform1 = BilinearForm(tspace)
bform1.add_domain_integrator(integrator1)
a0, A0 = integrator1.mass_matrix(space=tspace)
print("a0:", a0.shape, "\n", a0[0].round(4))
print("A0:", A0.shape, "\n", A0.toarray().round(4))

a1, A1 = integrator1.tr_mass_matrix(space=tspace)
print("a1:", a1.shape, "\n", a1[0].round(4))
print("A1:", A1.shape, "\n", A1.toarray().round(4))

A = A0 + A1
print("A:", A.shape, "\n", A.toarray().round(4))

b, B = integrator1.div_matrix(space=tspace, spaces=vspace)
print("b:", b.shape, "\n", b[0].round(4))
print("B:", B.shape, "\n", B.toarray().round(4))

integrator2 = VectorSourceIntegrator(f = pde.source, q=p+2)

lform = LinearForm(vspace)
lform.add_domain_integrator(integrator2)
f = integrator2.assembly_cell_vector(space = vspace)
print("f:", f.shape, "\n", f)
lform.assembly()
F = lform.get_vector()
print("F:", F.shape, "\n", F.round(4))

Z = coo_matrix((vgdof*GD, vgdof*GD))
BT = B.copy().T
Stiffness_matrix = bmat([[A, B], [BT, Z]], format='csr')
print("Stiffness_matrix:", Stiffness_matrix.shape, "\n", Stiffness_matrix.toarray())
g = np.zeros(tgdof)
F_rhs = np.hstack((g, F))
print("F_rhs:", F_rhs.shape)

Uh = np.hstack((sigmah, uh.flat))

isBdDof_1 = tspace.set_dirichlete_bc(uh=sigmah, gD=0)
print("isBdDof_1:", isBdDof_1.shape, "\n", isBdDof_1)
isBdDof_2 = np.zeros(vgdof*GD, dtype=bool)
print("isBdDof_2:", isBdDof_2.shape, "\n", isBdDof_2)
isBdDof = np.hstack((isBdDof_1, isBdDof_2))
print("isBdDof:", isBdDof.shape, "\n", isBdDof)

F_rhs -= Stiffness_matrix @ Uh
F_rhs[isBdDof] = Uh[isBdDof]
print("F_rhs:", F_rhs.shape)

bdIdx = np.zeros(Stiffness_matrix.shape[0], dtype=np.int_)
bdIdx[isBdDof] = 1
print("bdIdx:", bdIdx.shape, "\n", bdIdx)
D0 = spdiags(1-bdIdx, 0, Stiffness_matrix.shape[0], Stiffness_matrix.shape[0])
D1 = spdiags(bdIdx, 0, Stiffness_matrix.shape[0], Stiffness_matrix.shape[0])
Stiffness_matrix = D0@Stiffness_matrix@D0 + D1
print("Stiffness_matrix:", Stiffness_matrix.shape, "\n", Stiffness_matrix.toarray())

Uh.flat[:] = spsolve(Stiffness_matrix, F_rhs)
print("Uh:\n", Uh.shape, "\n", Uh)

sigmah = Uh.flat[:tgdof]
print("sigmah:", sigmah.shape, "\n", sigmah)
# 因为 uh 属于的间断 Lagrange 有限元空间的自由度排序是 sdofs，所以 uh 的形状应该为 (GD, vgdof)
uh[:] = Uh.flat[:-tgdof].reshape(2, -1)
print("uh:", uh.shape, "\n", uh)

print("error", mesh.error(pde.solution, uh))



import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_cell(axes, showindex=True, color='k', marker='s', markersize=8, fontsize=16, fontcolor='k')
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=8, fontsize=16, fontcolor='r')
mesh.find_edge(axes, showindex=True, color='g', marker='*', markersize=8, fontsize=16, fontcolor='g')
plt.show()

mesh.to_vtk(fname=fname)

