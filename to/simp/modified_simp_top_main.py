import numpy as np

from modified_simp_top import TopModifiedSimp

nelx = 60
nely = 20
ts = TopModifiedSimp(nelx=nelx, nely=nely)

# Initialize optimization parameterse
nelx, nely, volfrac, penal, rmin, ft = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin, ts._ft
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
#cell = np.flipud(cell)
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

fname = os.path.join(output, 'modified_simp_quad_mesh.vtu')
mesh.to_vtk(fname=fname)

nu = 0.3
U = ts.FE(nelx=nelx, nely=nely, nu=nu)


