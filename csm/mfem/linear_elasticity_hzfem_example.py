import argparse
import os

import sympy as sp
import numpy as np

from fealpy.mesh import TriangleMesh

from fealpy.functionspace.HuZhangFiniteElementSpace2D import HuZhangFiniteElementSpace
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

lam = 1
mu = 1

ln = sp.ln
x = sp.symbols('x0:2')
u = [-(1-x[0])*ln(1.5-x[0]),-(1-x[0])*ln(1.5-x[1])]

if bdtype == 'displacement': 
        pde = GenLinearElasticitymodel2D(u, x, lam=lam, mu=mu,
                Dirichletbd_n='(x0==1)|(x0==0)|(x1==0)|(x1==1)',
                Dirichletbd_t='(x0==1)|(x0==0)|(x1==0)|(x1==1)')

elif bdtype =='stress_and_displacement':
        pde = GenLinearElasticitymodel2D(u,x,lam=lam,mu=mu,
                Dirichletbd_n='(x1==0)|(x1==1)',Dirichletbd_t='(x1==0)|(x1==1)',
                Neumannbd_nn='(x0==1)|(x0==0)',Neumannbd_nt='(x0==1)|(x0==0)')

elif bdtype =='stress_and_displacement_corner_point':
        pde = GenLinearElasticitymodel2D(u,x,lam=lam,mu=mu,
                Dirichletbd_n='(x0==1)|(x1==1)',Dirichletbd_t='(x0==1)|(x1==1)',
                Neumannbd_nn='(x0==0)|(x1==0)',Neumannbd_nt='(x0==0)|(x1==0)')

mesh = TriangleMesh.from_box(box=pde.domain(), nx=2, ny=2)
mesh.uniform_refine(nrefine)

output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'TriangleMesh.vtu')

tspace = HuZhangFiniteElementSpace(mesh, degree)
vspace = LagrangeFESpace(mesh, degree-1, spacetype='D') 

tgdof = tspace.number_of_global_dofs()
vgdof = vspace.number_of_global_dofs()

sh = tspace.function()
uh = vspace.function(dim=GD)

mu = pde.mu
lambda_ = pde.lam
M = tspace.parallel_compliance_tensor_matrix(mu=mu, lam=lam)
print("M:", M.shape, "\n", M)

mesh.to_vtk(fname=fname)

