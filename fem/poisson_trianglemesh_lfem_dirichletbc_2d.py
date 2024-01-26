import argparse

import numpy as np

import matplotlib.pyplot as plt

from scipy.sparse import spdiags

from poisson_2d import CosCosData

from fealpy.mesh.triangle_mesh import TriangleMesh 

from fealpy.functionspace.lagrange_fe_space import LagrangeFESpace

from fealpy.fem.diffusion_integrator import DiffusionIntegrator 
from fealpy.fem.scalar_source_integrator import ScalarSourceIntegrator
from fealpy.fem.bilinear_form import BilinearForm
from fealpy.fem.linear_form import LinearForm
from fealpy.fem.dirichlet_bc import DirichletBC

from fealpy.solver import GAMGSolver

from fealpy.tools.show import showmultirate


# Argument Parsing
parser = argparse.ArgumentParser(description=
        """
        Finite element method on a TriangleMesh of arbitrary order.
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Degree of the Lagrange finite element space. Default is 1.')

parser.add_argument('--nx',
        default=10, type=int,
        help='Number of initial mesh divisions along x.')

parser.add_argument('--ny',
        default=10, type=int,
        help='Number of initial mesh divisions along y.')

parser.add_argument('--maxit',
        default=4, type=int,
        help='Number of times to refine the mesh and solve. Default is 4.')

args = parser.parse_args()

p = args.degree
nx = args.nx
ny = args.ny
maxit = args.maxit

# Initialize the problem with given true solution
pde = CosCosData()
domain = pde.domain()

# Create the initial triangle mesh
mesh = TriangleMesh.from_box(box = domain, nx = nx, ny = ny)

errorType = ['$|| u - u_h||_{\\Omega,0}$', 
        '$||\\nabla u - \\nabla u_h||_{\\Omega, 0}$']
errorMatrix = np.zeros((2, maxit), dtype=np.float64)
NDof = np.zeros(maxit, dtype=np.int_)

# Main loop for mesh refinement and solution
for i in range(maxit):
    print("The {}-th computation:".format(i))
    # Create Lagrange finite element space
    space = LagrangeFESpace(mesh, p = p, spacetype = 'C', doforder = 'vdims')
    NDof[i] = space.number_of_global_dofs()

    # Assemble the stiffness matrix
    bform = BilinearForm(space)
    bform.add_domain_integrator(DiffusionIntegrator(q = p+2))
    A = bform.assembly()

    # Assemble the load vector
    lform = LinearForm(space)
    lform.add_domain_integrator(ScalarSourceIntegrator(f = pde.source, q = p+2))
    F = lform.assembly()

    # Apply Dirichlet boundary conditions
    uh = space.function() 
    isDDof = space.boundary_interpolate(gD=pde.dirichlet, uh=uh, threshold=pde.is_dirichlet_boundary)

    F = F - A@uh.reshape(-1)

    bdIdx = np.zeros(A.shape[0], dtype=np.int_)
    bdIdx[isDDof.reshape(-1)] = 1
    D0 = spdiags(1-bdIdx, 0, A.shape[0], A.shape[0])
    D1 = spdiags(bdIdx, 0, A.shape[0], A.shape[0])
    A = D0@A@D0 + D1

    F[isDDof.reshape(-1)] = uh[isDDof].reshape(-1)
    #bc = DirichletBC(space = space, gD = pde.dirichlet) 
    #A, F = bc.apply(A, F, uh)

    # Solve the linear system using GAMG solver
    solver = GAMGSolver(ptype='W', sstep=2)
    solver.setup(A)
    uh[:] = solver.solve(F)

    # Compute the errors
    errorMatrix[0, i] = mesh.error(pde.solution, uh)
    errorMatrix[1, i] = mesh.error(pde.gradient, uh.grad_value)

    # Refine the mesh for the next iteration
    if i < maxit-1:
        mesh.uniform_refine()

# Print error with its type
for i, errType in enumerate(errorType):
    print(errType)
    print(errorMatrix[i])
    print('------')

# Compute error convergence rates
ratios = errorMatrix[:, 0:-1] / errorMatrix[:, 1:]
convergence_rates = np.log2(ratios)

# Error
showmultirate(plt, 1, NDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)

# Numerical solution
fig = plt.figure()
axes = fig.add_subplot(1, 2, 1, projection = '3d')
uh.add_plot(axes, cmap = 'rainbow')

plt.show()

