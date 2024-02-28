from truss_model import Truss_2d

from fealpy.functionspace import LagrangeFESpace

from fealpy.fem import BilinearForm
from fealpy.fem import TrussStructureIntegrator

from scipy.sparse.linalg import spsolve
from scipy.sparse import spdiags

import matplotlib.pyplot as plt
import numpy as np

pde = Truss_2d()

mesh = pde.init_mesh()
node = mesh.entity('node')
NN = mesh.number_of_nodes()
cell = mesh.entity('cell') 
NC = mesh.number_of_cells()
print("NN:", NN)
print("NC:", NC)

GD = mesh.geo_dimension()
print("GD:", GD)

p = 1
space = LagrangeFESpace(mesh, p=p, spacetype='C', doforder='vdims')

gdof = space.number_of_global_dofs() 
ldof = space.number_of_local_dofs()
print("gdof:", gdof)
print("ldof:", ldof)

vspace = GD*(space, )

uh = vspace[0].function(dim=GD) 

bform = BilinearForm(vspace)

E0 = pde.E # 杨氏模量
A0 = pde.A # 横截面积
bform.add_domain_integrator(TrussStructureIntegrator(E=E0, A=A0, q=p+2))
K = bform.assembly()
print("K:", K.shape, "\n", K.toarray().round(4))

F = np.zeros((uh.shape[0], GD), dtype=np.float64)
idx_f, f = mesh.meshdata['force_bc']
F[idx_f] = f
print("F:", F.shape, "\n", F)

idx_disp, disp = mesh.meshdata['disp_bc']
# 按分量处理自由度索引
dflag = idx_disp
print("dflag:", dflag.shape, "\n", dflag)

F = F.flat- K@uh.flat
bdIdx = np.zeros(K.shape[0], dtype=np.int_)
bdIdx[dflag.flat] = 1
D0 = spdiags(1-bdIdx, 0, K.shape[0], K.shape[0])
D1 = spdiags(bdIdx, 0, K.shape[0], K.shape[0])
K = D0@K@D0 + D1
F[dflag.flat] = uh.ravel()[dflag.flat]

# Solving the linear system
uh.flat[:] = spsolve(K, F)
print("uh:", uh.shape, "\n", uh)

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='g')
plt.show()

