import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from truss_model import Truss_3d

from fealpy.functionspace.lagrange_fe_space import LagrangeFESpace

from fealpy.fem.bilinear_form import BilinearForm
from fealpy.fem.truss_structure_integrator import TrussStructureIntegrator

from scipy.sparse.linalg import spsolve
from scipy.sparse import spdiags

# Argument Parsing
parser = argparse.ArgumentParser(description=
    """
    Finite element method on simplex meshes (triangles, tetrahedra)
    """)

parser.add_argument('--degree',
        default = 1, type = int,
        help = 'Degree of the Lagrange finite element space, default is 1.')

parser.add_argument('--GD',
        default = 3, type = int,
        help = 'Dimension of the model problem, default is 3D.')

parser.add_argument('--nrefine',
        default = 2, type = int,
        help = 'Number of times the initial mesh is refined, default is 2.')

parser.add_argument('--scale',
        default = 100, type = float,
        help = 'Mesh deformation factor, default is 100.')

parser.add_argument('--doforder',
        default = 'vdims', type = str,
        help = 'Convention for ordering of degrees of freedom, default is vdims.')

args = parser.parse_args()

# Assigning the parsed arguments to local variables for clarity
p = args.degree
GD = args.GD
n = args.nrefine
scale = args.scale
doforder = args.doforder

# Checking the validity of the parsed arguments
if GD not in [2, 3] or p < 0 or n < 0 or scale < 0 or doforder not in ['vdims', 'sdofs']:
    parser.print_help()
    exit(1)

# Initializing the truss model
pde = Truss_3d()
# Initializing the mesh for the problem
mesh = pde.init_mesh()

# Constructing LagrangeFESpace with scalar basis
space = LagrangeFESpace(mesh, p = p, doforder = doforder)
vspace = GD*(space, )

# Retrieving cross-sectional area and Young's modulus
A = pde.A
E = pde.E

# Constructing bilinear form representing the differential formulation of the problem
bform = BilinearForm(vspace)
bform.add_domain_integrator(TrussStructureIntegrator(E=E, A=A, q=p+2))
K = bform.assembly()
print("K:", K.shape, "\n", K.toarray().round(4))

# Load and boundary conditions
uh = space.function(dim = GD)
print("uh:", uh.shape, "\n", uh)
F = np.zeros((uh.shape[0], GD), dtype = np.float64)
print("F:", F.shape, "\n", F)
idx_f, f = mesh.meshdata['force_bc']
F[idx_f] = f
idx_disp, disp = mesh.meshdata['disp_bc']

# Applying Dirichlet boundary conditions
dflag = vspace[0].boundary_interpolate(gD=disp, uh=uh, threshold=idx_disp)
F = F.flat- K@uh.flat
bdIdx = np.zeros(K.shape[0], dtype=np.int_)
bdIdx[dflag.flat] = 1
D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
K = D0@K@D0 + D1
F[dflag.flat] = uh.ravel()[dflag.flat]
#bc = DirichletBC(vspace, disp, threshold = idx_disp)
#A, F = bc.apply(K, F.flat, uh)

# Solving the linear system
uh.flat[:] = spsolve(K, F)

# Visualization of the original and deformed truss
fig = plt.figure()
axes = fig.add_subplot(1, 1, 1, projection='3d')

# Adding a title to the plot
axes.set_title("Original vs. Deformed Truss Structure")

# Plotting the original truss
mesh.add_plot(axes, nodecolor='k', cellcolor='k', linewidth=0.5)

# Applying deformation to the mesh nodes
mesh.node += scale*uh

# Plotting the deformed truss
mesh.add_plot(axes, nodecolor='m', cellcolor='m', linewidth=1.0)

legend_elements = [Line2D([0], [0], color='k', lw=2, label='Original Truss'),
                   Line2D([0], [0], color='m', lw=2, label='Deformed Truss')]
axes.legend(handles=legend_elements, loc='upper right')
plt.show()

