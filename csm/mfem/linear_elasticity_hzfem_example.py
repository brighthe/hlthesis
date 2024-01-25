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

from scipy.sparse import bmat, coo_matrix
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
Stiffness_matrix = bmat([[A, B], [BT, Z]], format='coo')
print("Stiffness_matrix:", Stiffness_matrix.shape, "\n", Stiffness_matrix.toarray())
g = np.zeros(tgdof)
F_rhs = np.hstack((g, F))
print("F_rhs:", F_rhs.shape)

Uh = np.hstack((sigmah, uh.flat))

ipoints = vspace[0].interpolation_points()
ipoints_1 = node[cell].reshape(-1, GD)
print("ipoints:", ipoints.shape, "\n", ipoints)
print("ipoints_1-ipoints:", np.max(np.abs(ipoints_1-ipoints)))
#isDDof = vspace[0].is_boundary_dof()
#print("isDDof:", isDDof.shape, "\n", isDDof)
cell2dof = vspace[0].cell_to_dof()
cell2dof_1 = np.arange(NC*vldof).reshape(NC, vldof)
print("cell2dof:", cell2dof.shape, "\n", cell2dof)
print("cell2dof-cell2dof_1:", np.max(np.abs(cell2dof-cell2dof_1)))
#idx = mesh.ds.boundary_edge_index()
#print("idx:", idx)
#isBdDof = np.zeros(vgdof, dtype=np.bool_)
#edge2dof = vspace[0].edge_to_dof()
#isBdDof[edge2dof[idx]] = True
#print("isBdDof:", isBdDof)



import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_cell(axes, showindex=True, color='k', marker='s', markersize=8, fontsize=16, fontcolor='k')
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=8, fontsize=16, fontcolor='r')
mesh.find_edge(axes, showindex=True, color='g', marker='*', markersize=8, fontsize=16, fontcolor='g')
plt.show()

Uh.flat[:] = spsolve(Stiffness_matrix, F_rhs)
print("Uh:\n", Uh.shape, "\n", Uh)



#lam = 1
#mu = 1
#
#ln = sp.ln
#x = sp.symbols('x0:2')
#u = [-(1-x[0])*ln(1.5-x[0]),-(1-x[0])*ln(1.5-x[1])]
#
#if bdtype == 'displacement': 
#        pde = GenLinearElasticitymodel2D(u, x, lam=lam, mu=mu,
#                Dirichletbd_n='(x0==1)|(x0==0)|(x1==0)|(x1==1)',
#                Dirichletbd_t='(x0==1)|(x0==0)|(x1==0)|(x1==1)')
#
#elif bdtype =='stress_and_displacement':
#        pde = GenLinearElasticitymodel2D(u,x,lam=lam,mu=mu,
#                Dirichletbd_n='(x1==0)|(x1==1)',Dirichletbd_t='(x1==0)|(x1==1)',
#                Neumannbd_nn='(x0==1)|(x0==0)',Neumannbd_nt='(x0==1)|(x0==0)')
#
#elif bdtype =='stress_and_displacement_corner_point':
#        pde = GenLinearElasticitymodel2D(u,x,lam=lam,mu=mu,
#                Dirichletbd_n='(x0==1)|(x1==1)',Dirichletbd_t='(x0==1)|(x1==1)',
#                Neumannbd_nn='(x0==0)|(x1==0)',Neumannbd_nt='(x0==0)|(x1==0)')
#
#mesh = TriangleMesh.from_box(box=pde.domain(), nx=2, ny=2)
#mesh.uniform_refine(nrefine)
#
#output = './mesh/'
#if not os.path.exists(output):
#    os.makedirs(output)
#fname = os.path.join(output, 'TriangleMesh.vtu')
#
#tspace = HuZhangFiniteElementSpace(mesh, degree)
#vspace = LagrangeFESpace(mesh, degree-1, spacetype='D') 
#
#tgdof = tspace.number_of_global_dofs()
#vgdof = vspace.number_of_global_dofs()
#
#sh = tspace.function()
#uh = vspace.function(dim=GD)
#
#mu = pde.mu
#lambda_ = pde.lam
#M = tspace.parallel_compliance_tensor_matrix(mu=mu, lam=lam)
#print("M:", M.shape, "\n", M)

mesh.to_vtk(fname=fname)

