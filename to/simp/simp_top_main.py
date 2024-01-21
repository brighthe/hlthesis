import numpy as np

from simp_top import TopSimp

ts = TopSimp()

# Initialize optimization parameterse
nelx, nely, volfrac, penal, rmin = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

# Initialize design variable field to the volume fraction
x = np.full((nely, nelx), volfrac)
mesh.celldata['x'] = x.flatten('F') # 按列增加

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'quad_mesh.vtu')
mesh.to_vtk(fname=fname)

#KE = ts.lk()
#print("KE:", KE.shape, "\n", KE.round(4))

#from fealpy.functionspace import LagrangeFESpace as Space
#from fealpy.fem import LinearElasticityOperatorIntegrator
#from fealpy.fem import BilinearForm
#space = Space(mesh, p=1, doforder='vdims')
#vspace = 2*(space, )
#E, nu = 1.0, 0.3
#mu = E / (2 *(1 + nu))
#lambda_ = ( E * nu) / ((1 + nu) * (1 - 2 * nu))
#integrator1 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=5)
#bform = BilinearForm(vspace)
#bform.add_domain_integrator(integrator1)
#KK = integrator1.assembly_cell_matrix(space=vspace)
#print("KK:", KK.shape, "\n", KK[0].round(4))

loop = 0 # Iteration counter
change = 1.0 # Maximum change in design variables between iterations

    # Optimization loop, runs until the change is less than 1%
while change > 0.01:
    loop += 1
    xold = np.copy(x)

    # FE-Analysis: perform finite element analysis on the current design
    U = ts.FE(nelx, nely, penal, x)

    # Objective Function And Sensitivity Analysis
    KE = ts.lk() # Retrieve element stiffness matrix
    c = 0 # Initialize objective (compliance) to zero
    dc = np.zeros((nely, nelx)) # Initialize sensitivity array to zero

    # Loop over every element to calculate the objective and sensitivity
    for elx in range(nelx):
        for ely in range(nely):
            # Global node numbers for the upper left and upper right nodes of the element
            n1 = (nely+1) * elx + ely
            n2 = (nely+1) * (elx+1) + ely
            # Degrees of freedom for the element
            edof = np.array([2*n1, 2*n1 + 1, 2*n2, 2*n2 + 1, 2*n2 + 2, 2*n2 + 3, 2*n1 + 2, 2*n1 + 3])
            # Extract element displacements
            Ue = U[edof]
            # Update objective (compliance) and its sensitivity
            c += x[ely, elx]**penal * Ue.T @ KE @ Ue 
            dc[ely, elx] = -penal * x[ely, elx]**(penal - 1) * Ue.T @ KE @ Ue

    # Filtering of Sensitivity: apply mesh-independent filter to the sensitivities
    dc = ts.check(nelx, nely, rmin, x, dc)

    # Design Update By The Optimality Criteria Method
    x = ts.OC(nelx, nely, volfrac, x, dc)

    # Print Results: output the current iteration results
    change = np.max(np.abs(x - xold))
    print(f' Iter.: {loop:4d} Objective.: {c:10.4f} Volfrac.: {np.sum(x)/(nelx*nely):6.3f} change.: {change:6.3f}')




