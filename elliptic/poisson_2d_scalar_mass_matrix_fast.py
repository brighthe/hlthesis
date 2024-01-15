import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

from scipy.sparse.linalg import spsolve

from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.mesh import TriangleMesh

from fealpy.pde.poisson_2d import CosCosData

from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import VectorSourceIntegrator
from fealpy.fem import VectorMassIntegrator
from fealpy.fem import ScalarMassIntegrator
from fealpy.fem import BilinearForm
from fealpy.fem import LinearForm
from fealpy.fem import DirichletBC


# Argument Parsing
parser = argparse.ArgumentParser(description=
        """
        Finite element method on a TriangleMesh of arbitrary order.
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Degree of the Lagrange finite element space. Default is 1.')

parser.add_argument('--GD',
        default=2, type=int,
        help='模型问题的维数, 默认求解 2 维问题.')

parser.add_argument('--nx',
        default=2, type=int,
        help='Number of initial mesh divisions along x.')

parser.add_argument('--ny',
        default=2, type=int,
        help='Number of initial mesh divisions along y.')

args = parser.parse_args()

p = args.degree
GD = args.GD
nx = args.nx
ny = args.ny

# Initialize the problem with given true solution
pde = CosCosData()
domain = pde.domain()

# Create the initial triangle mesh
mesh = TriangleMesh.from_box(box = domain, nx = nx, ny = ny)
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
print("NC:", NC)
node = mesh.entity('node')
cell = mesh.entity('cell')

output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'TriangleMesh.vtu')

space = Space(mesh, p=p)
gdof = space.number_of_global_dofs()
ldof = space.number_of_local_dofs()
print("ldof", ldof)

from fealpy.decorator import cartesian
@cartesian
def func_coef(p):
    x = p[..., 0]
    y = p[..., 1]
    return x + y
    #return np.cos(x) + np.sin(y)
    #return x**2 + y**2
u = space.interpolate(func_coef)
cell2dof = space.cell_to_dof()
uc = u[cell2dof]
print("uc:", uc.shape)

integrator_scalar_mass = ScalarMassIntegrator(c=func_coef, q=p+3)

bform_scalar_mass_1 = BilinearForm(space)
bform_scalar_mass_1.add_domain_integrator(integrator_scalar_mass)
MK_scalar_mass_1 = integrator_scalar_mass.assembly_cell_matrix(space=space)
print("MK_scalar_mass_1:\n", MK_scalar_mass_1.shape, "\n", MK_scalar_mass_1[0])
bform_scalar_mass_1.assembly()
M_scalar_mass_1 = bform_scalar_mass_1.get_matrix()
#print("M_scalar_mass_1:\n", M_scalar_mass_1.shape, "\n", M_scalar_mass_1.toarray())

bform_scalar_mass_2 = BilinearForm(space)
bform_scalar_mass_2.add_domain_integrator(integrator_scalar_mass)
MK_scalar_mass_2 = integrator_scalar_mass.assembly_cell_matrix_fast(trialspace=space, testspace=space, coefspace=space)
print("MK_scalar_mass_2:\n", MK_scalar_mass_2.shape, "\n", MK_scalar_mass_2[0])
bform_scalar_mass_2.fast_assembly(trialspace=space, testspace=space, coefspace=space)
M_scalar_mass_2 = bform_scalar_mass_2.get_matrix()
#print("M_scalar_mass_2:\n", M_scalar_mass_2.shape, "\n", M_scalar_mass_2.toarray())

local_matrices_equal = MK_scalar_mass_1 - MK_scalar_mass_2
print("local_matrices_equal:\n", local_matrices_equal.round(4))
global_matricee_equal = M_scalar_mass_1 - M_scalar_mass_2
#print("global_matricee_equal:\n", global_matricee_equal.toarray().round(4))
#local_matrices_equal = np.allclose(MK_scalar_mass_1, MK_scalar_mass_2, rtol=1e-05, atol=1e-08)
#print("local_matrices_equal:\n", local_matrices_equal)
#
#global_matrices_equal = np.allclose(M_scalar_mass_1.toarray(), M_scalar_mass_2.toarray(), rtol=1e-05, atol=1e-08)
#print("global_matrices_equal:\n", global_matrices_equal)
