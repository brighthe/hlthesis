import argparse

import numpy as np

from poisson import SinSinData

from fealpy.mesh.triangle_mesh import TriangleMesh
from fealpy.mesh.polygon_mesh import PolygonMesh

from fealpy.functionspace.lagrange_fe_space import LagrangeFESpace

from fealpy.fem.scalar_mass_integrator import ScalarMassIntegrator
from fealpy.fem.bilinear_form import BilinearForm

from scipy.sparse import csr_matrix

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
pde = SinSinData()
domain = pde.domain()

# Create the initial triangle mesh
# mesh = TriangleMesh.from_box(box = domain, nx = nx, ny = ny)
node = np.array([[0.0, 0.0], [0.0, 0.5], [0.0, 1.0],
                 [0.5, 0.0], [0.5, 0.5], [0.5, 1.0],
                 [1.0, 0.0], [1.0, 0.5], [1.0, 1.0]], dtype=np.float64)
cell = np.array([0, 3, 4, 1, 3, 6, 7, 4, 1, 4, 5, 2, 4, 7, 8, 4, 8, 5], dtype=np.int_)
cellLocation = np.array([0, 4, 8, 12, 15, 18], dtype=np.int_)
mesh = PolygonMesh(node=node, cell=cell, cellLocation=cellLocation)
edge = mesh.entity('edge')
# mesh = PolygonMesh.from_box(box = domain, nx = nx, ny = ny)


import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=8, fontsize=16, fontcolor='r')
mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=8, fontsize=16, fontcolor='b')
mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=8, fontsize=16, fontcolor='g')
plt.show()

printt("sadasd")




space = LagrangeFESpace(mesh, p = p, spacetype = 'C', doforder = 'vdims')
#space = LagrangeFESpace(mesh, p = p, spacetype = 'C', doforder = 'sdofs')
gdof = space.number_of_global_dofs()

# Assemble the mass matrix
integrator = ScalarMassIntegrator(c=1, q=p+1)
bform = BilinearForm(space)
bform.add_domain_integrator(integrator)
MK = integrator.assembly_cell_matrix(space=space)
print("MK:", MK.shape, "\n", MK.round(4))
bform.assembly()
M = bform.get_matrix()
print("M:", M.shape, "\n", M.toarray().round(4))

qf = mesh.integrator(p+1) 
bcs, ws = qf.get_quadrature_points_and_weights()
print("ws:", ws.shape)
#import ipdb
#ipdb.set_trace()
phi = space.basis(bcs)
gphi = space.grad_basis(bcs)
print("phi", phi.shape)
print("gard_phi", gphi.shape)
cm = mesh.entity_measure('cell')
print("cm:", cm.shape)
MK_1 = np.einsum('q, qcj, qck, c -> cjk', ws, phi, phi, cm)
print("MK_1:", MK_1.shape, "\n", MK_1.round(4))
cell2dof = space.cell_to_dof()
I = np.broadcast_to(cell2dof[:, :, None], shape=MK_1.shape)
J = np.broadcast_to(cell2dof[:, None, :], shape=MK_1.shape)
M_1 = csr_matrix((MK_1.flat, (I.flat, J.flat)), shape=(gdof, gdof))
print("M_1:", M_1.shape, "\n", M_1.toarray().round(4))

#A = np.einsum('i, ijkl, ijml, j -> jkm', ws, gphi, gphi, cm)


