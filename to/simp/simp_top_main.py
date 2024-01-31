import numpy as np

from simp_top import TopSimp

from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

nelx = 2
nely = 2
volfrac = 0.5
penal = 3.0
ts = TopSimp(nelx=nelx, nely=nely, volfrac=volfrac, penal=penal)

# 初始化优化参数
nelx, nely, volfrac, penal, rmin = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin
mesh = ts._mesh


node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

fname = os.path.join(output, 'simp_quad_mesh.vtu')
mesh.to_vtk(fname=fname)

# 根据体积分数 volfrac 初始化设计变量场
x = np.full((nely, nelx), volfrac)
print("x:", x.shape, "\n", x)

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

E0 = 1.0
nu = 0.3

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
from fealpy.fem import BilinearForm
from fealpy.functionspace import LagrangeFESpace as Space

p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
uh = space.function(dim=GD)
print("uh:", uh.shape)
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
ldof = vspace[0].number_of_local_dofs()
vldof = ldof * GD
print("vgdof", vgdof)
print("vldof", vldof)

integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, penal=penal, x=x)
bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator)
KK = integrator.assembly_cell_matrix(space=vspace)
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))

# 定义荷载
F = np.zeros(vgdof)
F[1] = -1
print("F:", F.shape, "\n", F.round(4))

# 定义支撑(边界处理)
fixeddofs = np.union1d( np.arange(0, 2*(nely+1), 2), np.array([2*(nelx+1)*(nely+1) - 1]) )
dflag = fixeddofs
F = F - K@uh.flat
bdIdx = np.zeros(K.shape[0], dtype=np.int_)
bdIdx[dflag.flat] = 1
D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
K = D0@K@D0 + D1
F[dflag.flat] = uh.ravel()[dflag.flat]

# 线性方程组求解
uh.flat[:] = spsolve(K, F)
print("uh:", uh.shape, "\n", uh)

import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
plt.show()

#loop = 0 # Iteration counter
#change = 1.0 # Maximum change in design variables between iterations
#
## Optimization loop, runs until the change is less than 1%
#while change > 0.01:
#    loop += 1
#    xold = np.copy(x)
#
#    # FE-Analysis: perform finite element analysis on the current design
#    U = ts.FE(nelx, nely, penal, x)
#
#    # Objective Function And Sensitivity Analysis
#    KE = ts.lk() # Retrieve element stiffness matrix
#    c = 0 # Initialize objective (compliance) to zero
#    dc = np.zeros((nely, nelx)) # Initialize sensitivity array to zero
#
#    # Loop over every element to calculate the objective and sensitivity
#    for elx in range(nelx):
#        for ely in range(nely):
#            # Global node numbers for the upper left and upper right nodes of the element
#            n1 = (nely+1) * elx + ely
#            n2 = (nely+1) * (elx+1) + ely
#            # Degrees of freedom for the element
#            edof = np.array([2*n1, 2*n1 + 1, 2*n2, 2*n2 + 1, 2*n2 + 2, 2*n2 + 3, 2*n1 + 2, 2*n1 + 3])
#            # Extract element displacements
#            Ue = U[edof]
#            # Update objective (compliance) and its sensitivity
#            c += x[ely, elx]**penal * Ue.T @ KE @ Ue 
#            dc[ely, elx] = -penal * x[ely, elx]**(penal - 1) * Ue.T @ KE @ Ue
#
#    # Filtering of Sensitivity: apply mesh-independent filter to the sensitivities
#    dc = ts.check(nelx, nely, rmin, x, dc)
#
#    # Design Update By The Optimality Criteria Method
#    x = ts.OC(nelx, nely, volfrac, x, dc)
#
#    # Print Results: output the current iteration results
#    change = np.max(np.abs(x - xold))
#    print(f' Iter.: {loop:4d} Objective.: {c:10.4f} Volfrac.: {np.sum(x)/(nelx*nely):6.3f} change.: {change:6.3f}')
#
#    mesh.celldata['x'] = np.flipud(x).flatten('F') # 因为网格的左下角为 0 号单元，所以先将 x 翻转，再按列展平
#
#    fname = os.path.join(output, f'quad_mesh_{loop:04d}.vtu')
#    mesh.to_vtk(fname=fname)
#
