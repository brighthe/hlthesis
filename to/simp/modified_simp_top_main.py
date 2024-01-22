import numpy as np

from modified_simp_top import TopModifiedSimp

nelx = 4
nely = 3
ts = TopModifiedSimp(nelx=nelx, nely=nely)

# Initialize optimization parameterse
nelx, nely, volfrac, penal, rmin, ft = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin, ts._ft
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针，单元从下往上
#cell = np.flipud(cell)
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

fname = os.path.join(output, 'modified_simp_quad_mesh.vtu')
mesh.to_vtk(fname=fname)

# Initialize design variable field to the volume fraction
x = np.full((nely, nelx), volfrac)
xPhys = x

loop = 0 # Iteration counter
change = 1.0 # Maximum change in design variables between iterations

E0 = 1.0
Emin = 1e-9
nu = 0.3

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
from fealpy.fem import BilinearForm
from fealpy.functionspace import LagrangeFESpace as Space

space = Space(mesh, p=1, doforder='vdims')
GD = 2
uh = space.function(dim=GD)
vspace = GD*(space, )
gdof = vspace[0].number_of_global_dofs()
vgdof = gdof * GD
ldof = vspace[0].number_of_local_dofs()
vldof = ldof * GD
print("vgdof", vgdof)
print("vldof", vldof)

integrator = MbbBeamOperatorIntegrator(nu=nu, nelx=nelx, nely=nely, xPhys=xPhys, penal=penal, E0=E0, Emin=Emin)
bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator)
sK = integrator.assembly_cell_matrix(space=vspace)
print("sK", sK.shape, "\n", sK)
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))


#nu = 0.3
#U = ts.FE(nelx=nelx, nely=nely, nu=nu)

# Optimization loop, runs until the change is less than 1%
#while change > 0.01:
#    loop += 1

