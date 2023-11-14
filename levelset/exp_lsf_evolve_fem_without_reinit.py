import argparse 
import os
import numpy as np

from fealpy.mesh.triangle_mesh import TriangleMesh
from fealpy.timeintegratoralg import UniformTimeLine
from fealpy.functionspace import LagrangeFESpace
from fealpy.levelset.ls_fem_solver import LSFEMSolver, LSSolver
from fealpy.decorator import cartesian


# Command line argument parser
parser = argparse.ArgumentParser(description=
        """
        Finite element method to solve the level set evolution equation with Crank-Nicholson time discretization.
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Degree of the Lagrange finite element space. Default is 1.')

parser.add_argument('--ns',
        default=100, type=int,
        help='Number of spatial divisions in each direction. Default is 100.')

parser.add_argument('--nt',
        default=100, type=int,
        help='Number of time divisions. Default is 100.')

parser.add_argument('--T',
        default=1, type=float,
        help='End time of the evolution. Default is 1.')

parser.add_argument('--output',
        default='./results/', type=str,
        help='Output directory for the results. Default is ./results/')
        
parser.add_argument('--step',
        default=10, type=int,
        help='')

args = parser.parse_args()

# Check if the directory exists, if not, create it
if not os.path.exists(args.output):
    os.makedirs(args.output)

degree = args.degree
nt = args.nt
ns = args.ns
T = args.T
output = args.output

# Define the velocity field $u$ for the evolution
@cartesian
def velocity_field(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = np.sin((np.pi*x))**2 * np.sin(2*np.pi*y)
    u[..., 1] = -np.sin((np.pi*y))**2 * np.sin(2*np.pi*x)
    return u

# Initial level set function $\phi0$ representing the circle
@cartesian
def circle(p):
    x = p[...,0]
    y = p[...,1]
    val = np.sqrt((x-0.5)**2+(y-0.75)**2)-0.15
    return val

# Define the domain and generate the triangular mesh
domain = [0, 1, 0, 1]
mesh = TriangleMesh.from_box(domain, nx=ns, ny=ns)
cellmeasure = mesh.entity_measure('cell')

# Generate the uniform timeline
timeline = UniformTimeLine(0, T, nt)
dt = timeline.dt

# Define the finite element space
space = LagrangeFESpace(mesh, p=degree)

# Initialize the level set function $phi0$ and velocity field $u$ on the mesh nodes
phi0 = space.interpolate(circle)
u = space.interpolate(velocity_field, dim=2)

lsfemsolver = LSFEMSolver(space = space, u = u)

lssolver = LSSolver(space = space, phi0 = phi0, u = u)

# If output is enabled, save the initial state
lssolver.output(timestep = 0, output_dir = output, filename_prefix = 'lsf_without_reinit')

diff_avg, diff_max = lssolver.check_gradient_norm(phi = phi0)
print(f"Average diff: {diff_avg:.4f}, Max diff: {diff_max:.4f}")

# Time iteration
for i in range(nt):
    t1 = timeline.next_time_level()
    print("t1=", t1)

    phi0[:] = lsfemsolver.solve(phi0 = phi0, dt = dt)
    print("phi0", phi0)

    # Save the current state if output is enabled
    lssolver.output(timestep = i+1, output_dir = output, filename_prefix = 'lsf_without_reinit')

    # Move to the next time level
    timeline.advance()

diff_avg, diff_max = lssolver.check_gradient_norm(phi = phi0)
print(f"Average diff: {diff_avg:.4f}, Max diff: {diff_max:.4f}")

