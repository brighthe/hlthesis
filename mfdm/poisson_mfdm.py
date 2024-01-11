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

node = mesh.entity('node')
NN = mesh.number_of_nodes()
print("node:", NN, ":\n", node)
cell = mesh.entity('cell')
NC = mesh.number_of_cells()
print("cell:", NC, ":\n", cell)
edge = mesh.entity('edge')
NE = mesh.number_of_edges()
print("edge:", NE, ":\n", edge)

edge_lengths = mesh.entity_measure('edge')
print("edge_lengths:\n", edge_lengths)

cell2edge = mesh.ds.cell_to_edge()
print("cell2edge:\n", cell2edge)
cell2edge_sorted = [np.array(sorted(arr)) for arr in cell2edge]
print("cell2edge_sorted:\n", cell2edge_sorted)

#vertices = mesh.ds.cell_to_node()
#print("vertices:\n", vertices)

area_of_p = mesh.entity_measure('cell')
print("area_of_p:\n", area_of_p)

edge_normals = mesh.edge_unit_normal()
edge_normals = -edge_normals
print("edge_unit_normal:\n", edge_normals)

edge_centers = mesh.entity_barycenter('edge')
print("edge_centers:\n", edge_centers)

cell2edge_sign = mesh.ds.cell_to_edge_sign(return_sparse=False)
print("cell2edge_sign:\n", cell2edge_sign)
Alphas = []
start = 0
for c in cell:
    end = start + len(c)
    Alphas.append(cell2edge_sign[start:end])
    start = end
print("Alphas:\n", Alphas)

barycenters = mesh.entity_barycenter('cell')
print("barycenters:\n", barycenters)

Mf = np.zeros((NE, NE))
# TODO: sparse

for i in range(NC):
    num_edges = len(cell2edge[i])
    R = np.zeros((num_edges, 2))
    N = np.zeros((num_edges, 2))
    for j, edge_index in enumerate(cell2edge_sorted[i]):
        R[j, :] = Alphas[i][j] * (edge_centers[edge_index] - barycenters[i,:]) * edge_lengths[edge_index]
        N[j, :] = edge_normals[edge_index]

    print("R:\n", R)
    print("N:\n", N)

    size = 1 / num_edges
    M0 = R @(np.linalg.pinv(R.T @ N) @ R.T)
    print("M0:\n", M0)
    M1 = size * np.trace(M0) * (np.eye(num_edges) - N @ np.linalg.pinv(N.T @ N) @ N.T)
    print("M1:\n", M1)
    M = M0 + M1
    print("M:\n", M)

    assembly = np.zeros((num_edges, NE))
    for k, edge in enumerate(cell2edge_sorted[i]):
        #print("k:", k)
        #print("edge:", edge)
        assembly[k, edge] = 1

    #print("assembly:\n", assembly)
    Mf += assembly.T @ M @ assembly
    print("Mf:\n", Mf)

divh = np.zeros((NC, NE))
for i in range(NC):
    for j, edge_index in enumerate(cell2edge_sorted[i]):
        divh[i, edge_index] = Alphas[i][j] * edge_lengths[edge_index] / area_of_p[i]
print("divh:\n", divh)

Mc = np.diag(area_of_p)
print("Mc:\n", Mc)

A = np.block([
    [Mf, -divh.T @ Mc],
    [Mc @ divh, np.zeros((NC, NC))]
])
print("A:", A.shape, "\n", A)


import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=8, fontsize=16, fontcolor='r')
mesh.find_cell(axes, showindex=True, color='b', marker='o', markersize=8, fontsize=16, fontcolor='b')
mesh.find_edge(axes, showindex=True, color='g', marker='o', markersize=8, fontsize=16, fontcolor='g')
plt.show()

import os
output = './csm_results/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'polygon_mesh.vtu')
mesh.to_vtk(fname=fname)

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


