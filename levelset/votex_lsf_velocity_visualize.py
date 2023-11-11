import argparse 
import numpy as np

from functools import partial

from fealpy.mesh.triangle_mesh import TriangleMesh
from fealpy.functionspace import LagrangeFESpace
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

args = parser.parse_args()


degree = args.degree
nt = args.nt
ns = args.ns
T = args.T


# Define the velocity field $u$ for the evolution
@cartesian
def velocity_field(p, t):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = 2 * np.sin(2 * np.pi * y) * np.sin(np.pi * x)**2 * np.cos(np.pi * t)
    u[..., 1] = -2 * np.sin(2 * np.pi * x) * np.sin(np.pi * y)**2 * np.cos(np.pi * t)
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

# Define the finite element space
space = LagrangeFESpace(mesh, p=degree)

# Initialize the level set function $phi0$ and velocity field $u$ on the mesh nodes
phi0 = space.interpolate(circle)
u = space.interpolate(partial(velocity_field, t=0.5), dim=2)

import matplotlib.pyplot as plt

# Extract the coordinates from the mesh
X = mesh.node[:, 0].reshape((ns+1, ns+1))
Y = mesh.node[:, 1].reshape((ns+1, ns+1))

# Reshape phi0 and u for plotting
phi0_reshaped = phi0.reshape(ns+1, ns+1)
u_reshaped = u.reshape(ns+1, ns+1, 2)

# Plotting the level set phi0
fig, ax = plt.subplots(figsize=(10, 8))

# Plotting the velocity field u using quiver (arrows)
skip = 5  # Adjust as needed
quiver = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip], np.array(u_reshaped[::skip, ::skip, 0]), np.array(u_reshaped[::skip, ::skip, 1]), color='blue')

# Add dashed lines from circle center to axes
circle_center = (0.5, 0.75)
ax.axvline(x=circle_center[0], ymin=0, ymax=circle_center[1], color='black', linestyle='--', linewidth=0.5)
ax.axhline(y=circle_center[1], xmin=0, xmax=circle_center[0], color='black', linestyle='--', linewidth=0.5)

# Add circle
circle = plt.Circle((0.5, 0.75), 0.15, color='black', fill=False)
ax.add_patch(circle)
ax.annotate('r=0.15', xy=circle_center, fontsize=12, weight='bold', ha='right', va='bottom')

ax.set_title("Initial configuration and domain for the vortex test(t=0.5)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
plt.tight_layout()
plt.show()




