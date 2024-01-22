import argparse
import os

import sympy as sp
import numpy as np

from fealpy.mesh import TriangleMesh

from fealpy.fem import BilinearForm

from HuZhangFiniteElementSpace2d import HuZhangFiniteElementSpace2d
from linear_elasticity_model2d import BoxDomainData2d
from hu_zhang_mass_operator_integrator import HuZhangMassOperatorIntegrator



from fealpy.functionspace.lagrange_fe_space import LagrangeFESpace

from fealpy.pde.linear_elasticity_model2D import GenLinearElasticitymodel2D

parser = argparse.ArgumentParser(description=
        """
        三角形网格上用胡张元求解线弹性力学问题
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

parser.add_argument('--maxit',
        default=1, type=int,
        help='默认网格加密求解的次数, 默认加密求解 1 次')

parser.add_argument('--bdtype',
        default='displacement', type=str,
        help='边界条件, 默认为位移边界')

args = parser.parse_args()

degree = args.degree
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

p = 1
spacetype = 'D'
doforder = 'sdofs'
space = LagrangeFESpace(mesh, p=p, spacetype=spacetype, doforder=doforder)
vspace = GD*(space, )

tspace = HuZhangFiniteElementSpace2d(mesh, p=p)
tldof = tspace.number_of_local_dofs()
print("tldof:", tldof)
tgdof = tspace.number_of_global_dofs()
print("tgdof:", tgdof)
tgdof = tspace.number_of_global_dofs()
#bcs, ws = tspace.integrator.get_quadrature_points_and_weights()
#phi = tspace.basis(bcs)
integrator1 = HuZhangMassOperatorIntegrator(mesh ,p=p)
bform = BilinearForm(tspace)
bform.add_domain_integrator(integrator1)
a0, A0 = integrator1.mass_matrix(space=tspace)
print("a0:", a0.shape, "\n", a0.round(4))
print("A0:", A0.shape, "\n", A0.toarray().round(4))
#bform.assembly()
#K = bform.get_matrix()





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

