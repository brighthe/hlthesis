import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from mimetic_solver import Mimetic
from lsf_model import ClassicalLsfData

pde = ClassicalLsfData()
#mesh = pde.polygon_mesh()
ns=5
mesh = pde.polygon_mesh_2(n=ns)
#mesh = pde.triangle_mesh(domain=[0, 1, 0, 1], nx=10, ny=10)


maxit = 5
errorType = ['$|| u - u_h||_{\\Omega,0}$', 
        '$||\\nabla u - \\nabla u_h||_{\\Omega, 0}$']
errorMatrix = np.zeros((1, maxit), dtype=np.float64)
NDof = np.zeros(maxit, dtype=np.int_)

for iter in range(maxit):

    NN = mesh.number_of_nodes()
    print("NN:", NN)
    NE = mesh.number_of_edges()
    print("NE:", NE)
    NC = mesh.number_of_cells()
    print("NC:", NC)

    solver = Mimetic(mesh)
    M_c = solver.M_c()
    M_f = solver.M_f()
    div_operator = solver.div_operator()
    gradh = solver.gard_operator()

    node = mesh.entity('node')
    cell = mesh.entity('cell')
    edge = mesh.entity('edge')

    edge_norm = mesh.edge_unit_normal()
    cell_measure = mesh.entity_measure('cell') 

    A10 = solver.u_M_f(velocity=pde.velocity_field)
    #print("A10:", A10.shape, "\n", A10)

    exact_values = mesh.integral(u=pde.scalar_product, q=3, celltype=True)
    #print("exact_values:", exact_values.shape, "\n", exact_values)

    K_nodes = pde.grad_circle(node)
    cell2edge = mesh.ds.cell_to_edge()
    flag = np.where(mesh.ds.cell_to_edge_sign().toarray(), 1, -1)
    result = []
    K_h = np.zeros(NE, )
    for i in range(NC):
        edge_c = edge[cell[i]] # (LNE, 2)
        K_edges_c = (K_nodes[edge_c[:, 0]] + K_nodes[edge_c[:, 1]]) / 2 # (LNE, 2)
        tmp1 = edge_norm[cell2edge[i]]
        tmp2 = flag[i, cell2edge[i]].reshape(-1, 1)
        edge_norm_c = tmp1 * tmp2
        K_c = np.einsum('ij, ij -> i', K_edges_c, edge_norm_c) # (LNE, )
        result.append(K_c)
        K_h[cell2edge[i]] = K_c

    #print("result:\n", result)
    #print("test:", test.shape, "\n", test)

    #approx_values = np.zeros(NC, )

    #for i, cell_edges in enumerate(cell):
    #    A_row = A10[i]
    #    result_values = result[i]
    #    multiply_result = sum(A_row[cell_edges] * result_values)
    #    approx_values[i] = multiply_result
    ##print("approx_values:", approx_values.shape, "\n", approx_values)

    approx_values_new = A10 @ K_h
    #print("approx_values_new:", approx_values_new.shape, "\n", approx_values_new)

    errorMatrix[0, iter] = np.max(np.abs(exact_values - approx_values_new))
    print("errorMatrix:", errorMatrix)

    #error = np.max(np.abs(exact_values - approx_values))
    #print("error:", error)

    if iter < maxit-1:
        #mesh.uniform_refine()
        print("iter:", iter)
        ns = ns*2
        mesh = pde.polygon_mesh_2(n=ns)

    NDof[iter] = NN

import matplotlib.pyplot as plt
from fealpy.tools.show import showmultirate

showmultirate(plt, 2, NDof, errorMatrix, errorType, propsize=20, lw=2, ms=4)
plt.show()
