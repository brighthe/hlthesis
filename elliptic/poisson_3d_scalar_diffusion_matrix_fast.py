import argparse
import numpy as np

from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.mesh import TetrahedronMesh

from poisson_3d import CosCosCosData

from fealpy.fem import ScalarDiffusionIntegrator
from fealpy.fem import BilinearForm


# Argument Parsing
parser = argparse.ArgumentParser(description=
        """
        Finite element method on a TetrahedronMesh of arbitrary order.
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Degree of the Lagrange finite element space. Default is 1.')

parser.add_argument('--GD',
        default=3, type=int,
        help='模型问题的维数, 默认求解 3 维问题.')

parser.add_argument('--nx',
        default=2, type=int,
        help='Number of initial mesh divisions along x.')

parser.add_argument('--ny',
        default=2, type=int,
        help='Number of initial mesh divisions along y.')

parser.add_argument('--nz',
        default=2, type=int,
        help='Number of initial mesh divisions along z.')

args = parser.parse_args()

p = args.degree
GD = args.GD
nx = args.nx
ny = args.ny
nz = args.nz

# Initialize the problem with given true solution
pde = CosCosCosData()
domain = pde.domain()

# Create the initial triangle mesh
mesh = TetrahedronMesh.from_box(box = domain, nx = nx, ny = ny, nz = nz)
NN = mesh.number_of_nodes()
NC = mesh.number_of_cells()
print("NC:", NC)
node = mesh.entity('node')
cell = mesh.entity('cell')

space = Space(mesh, p=p)
gdof = space.number_of_global_dofs()
ldof = space.number_of_local_dofs()
print("ldof", ldof)

scalar_coef = 2
arr_coef = np.full(NC, 2)
from fealpy.decorator import cartesian
@cartesian
def func_coef(p):
    x = p[..., 0]
    y = p[..., 1]
    z = p[..., 2]
    return x + y + z

#integrator_scalar_diffusion = ScalarDiffusionIntegrator(q=p+2)
#integrator_scalar_diffusion = ScalarDiffusionIntegrator(c=scalar_coef, q=p+2)
#integrator_scalar_diffusion = ScalarDiffusionIntegrator(c=arr_coef, q=p+2)
integrator_scalar_diffusion = ScalarDiffusionIntegrator(c=func_coef, q=p+2)

bform_scalar_diffusion_1 = BilinearForm(space)
bform_scalar_diffusion_1.add_domain_integrator(integrator_scalar_diffusion)
MK_scalar_diffusion_1 = integrator_scalar_diffusion.assembly_cell_matrix(space=space)
print("MK_scalar_diffusion_1:\n", MK_scalar_diffusion_1.shape, "\n", MK_scalar_diffusion_1[0])
bform_scalar_diffusion_1.assembly()
M_scalar_diffusion_1 = bform_scalar_diffusion_1.get_matrix()
print("M_scalar_diffusion_1:\n", M_scalar_diffusion_1.shape, "\n", M_scalar_diffusion_1.toarray())

bform_scalar_diffusion_2 = BilinearForm(space)
bform_scalar_diffusion_2.add_domain_integrator(integrator_scalar_diffusion)
MK_scalar_diffusion_2 = integrator_scalar_diffusion.assembly_cell_matrix_fast(space=space)
print("MK_scalar_diffusion_2:\n", MK_scalar_diffusion_2.shape, "\n", MK_scalar_diffusion_2[0])
bform_scalar_diffusion_2.fast_assembly()
M_scalar_diffusion_2 = bform_scalar_diffusion_2.get_matrix()
print("M_scalar_diffusion_2:\n", M_scalar_diffusion_2.shape, "\n", M_scalar_diffusion_2.toarray())

local_matrices_equal = np.allclose(MK_scalar_diffusion_1, MK_scalar_diffusion_2, rtol=1e-05, atol=1e-08)
print("local_matrices_equal:\n", local_matrices_equal)

global_matrices_equal = np.allclose(M_scalar_diffusion_1.toarray(), M_scalar_diffusion_2.toarray(), rtol=1e-05, atol=1e-08)
print("global_matrices_equal:\n", global_matrices_equal)

