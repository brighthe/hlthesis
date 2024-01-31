import numpy as np
import matplotlib.pyplot as plt

from modified_simp_top import TopModifiedSimp

from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

nelx = 4
nely = 3
volfrac = 0.5
penal = 3.0
ts = TopModifiedSimp(nelx=nelx, nely=nely, volfrac=volfrac, penal=penal)

# 初始化优化参数
nelx, nely, volfrac, penal, rmin, ft = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin, ts._ft
mesh = ts._mesh

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
plt.show()
node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针，单元从下往上
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

fname = os.path.join(output, 'modified_simp_quad_mesh.vtu')
mesh.to_vtk(fname=fname)

# 根据体积分数 volfrac 初始化设计变量场
x = np.full((nely, nelx), volfrac)
xPhys = x

E0 = 1.0
Emin = 1e-9
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

integrator1 = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, 
                                    penal=penal, xPhys=xPhys, Emin=Emin)
bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator1)
KK = integrator1.assembly_cell_matrix(space=vspace)
#print("sK", sK.shape, "\n", sK[0].round(4))
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))

# 定义载荷
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
