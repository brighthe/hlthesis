import numpy as np

from lsf_top import TopLsf

from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

nelx = 2
nely = 2
volReq = 0.4
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq)

# 初始化优化参数
nelx, nely, volReq = ts._nelx, ts._nely, ts._volReq
mesh = ts._mesh
#nelx, nely, volReq, stepLength, numReinit, topWeight = ts._nelx, ts._nely, ts._volReq, ts._stepLength, ts._numReinit, ts._topWeight

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

# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))
print("struc:", struc.shape, "\n", struc)

#strucFull = np.zeros((struc.shape[0] + 2, struc.shape[1] + 2))
strucFull = np.zeros((nely + 2, nelx + 2))
strucFull[1:-1, 1:-1] = struc

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

integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator)
KK = integrator.assembly_cell_matrix(space=vspace)
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))

F = np.zeros(vgdof)
F[vgdof-1] = 1
print("F:", F.shape, "\n", F.round(4))

# 定义支撑(边界处理)
fixeddofs = np.arange(0, 2*(nely+1), 1)
dflag = fixeddofs
print("dflag:", dflag)
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

