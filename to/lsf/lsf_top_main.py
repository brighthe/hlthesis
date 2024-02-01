import numpy as np

from lsf_top import TopLsf

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


# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))
print("struc:", struc.shape, "\n", struc)

strucFull = np.zeros((nely + 2, nelx + 2))
strucFull[1:-1, 1:-1] = struc

# 初始化水平集函数
lsf = ts.reinit(strucFull = strucFull)
print("lsf:", lsf.shape, "\n", lsf.round(4))

# 初始化灵敏度
shapeSens = np.zeros((nely, nelx))
topSens = np.zeros((nely, nelx))

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
E0 = 1.0
nu = 0.3
integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
KE = integrator.matrix()

# 优化循环
num = 1
for iterNum in range(num):
    U, Ue = ts.FE(mesh=mesh, struc=struc)

    temp1 = -np.maximum(struc, 0.0001)
    temp2 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nely, nelx).T
    shapeSens[:] = np.einsum('ij, ij -> ij', temp1, temp2)
    print("shapeSens:", shapeSens.shape, "\n", shapeSens)


# 可视化
import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'lsf_quad_mesh.vtu')
mesh.celldata['strucFull'] = strucFull.flatten('F') # 按列增加
mesh.celldata['lsf'] = lsf.flatten('F') # 按列增加
mesh.to_vtk(fname=fname)



import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
plt.show()

