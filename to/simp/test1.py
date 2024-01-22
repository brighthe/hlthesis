#from fealpy.mesh import QuadrangleMesh
#nelx = 4
#nely = 3
#mesh = QuadrangleMesh.from_box(box = [0, nelx+1, 0, nely+1], nx = nelx, ny = nely)
#
#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_cell(axes, showindex=True, color='k', marker='s', markersize=6, fontsize=15, fontcolor='k')
#mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=6, fontsize=15, fontcolor='r')
#plt.show()

import numpy as np

# 定义变量
Emin = 1e-9  # 假设的最小杨氏模量值
E0 = 1.0    # 假设的最大杨氏模量值
penal = 3.0 # 惩罚因子
nelx, nely = 4, 3  # 单元的数量
KE = np.array([  # 给定的8x8单元刚度矩阵
    [ 0.4945,  0.1786, -0.3022, -0.0137, -0.2473, -0.1786,  0.0549,  0.0137],
    [ 0.1786,  0.4945,  0.0137,  0.0549, -0.1786, -0.2473, -0.0137, -0.3022],
    [-0.3022,  0.0137,  0.4945, -0.1786,  0.0549, -0.0137, -0.2473,  0.1786],
    [-0.0137,  0.0549, -0.1786,  0.4945,  0.0137, -0.3022,  0.1786, -0.2473],
    [-0.2473, -0.1786,  0.0549,  0.0137,  0.4945,  0.1786, -0.3022, -0.0137],
    [-0.1786, -0.2473, -0.0137, -0.3022,  0.1786,  0.4945,  0.0137,  0.0549],
    [ 0.0549, -0.0137, -0.2473,  0.1786, -0.3022,  0.0137,  0.4945, -0.1786],
    [ 0.0137, -0.3022,  0.1786, -0.2473, -0.0137,  0.0549, -0.1786,  0.4945]
])
xPhys = np.full((nely, nelx), 0.5)  # 物理密度矩阵

# 计算每个单元的杨氏模量
E = Emin + xPhys.ravel()**penal * (E0 - Emin)

# 将每个单元的杨氏模量乘以单元刚度矩阵
sK = np.einsum('i, jk -> ijk', E, KE)

# sK现在是一个三维数组，其中sK[i]是第i个单元的刚度矩阵
# 您可以通过sK[i]来访问第i个单元的刚度矩阵，它保持了8x8的形状

# 如果您想要将这些矩阵堆叠成一个大矩阵（可能是为了与全局刚度矩阵组装），可以这样做：
# sK_stacked = np.vstack(sK)

# 检查sK的形状是否符合预期
print("sK:", sK.shape, "\n", sK)

from fealpy.fem.mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
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

integrator1 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+1)

bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator1)
KK = integrator1.assembly_cell_matrix(space=vspace)
print("KK", KK.shape)
bform.assembly()
K = bform.get_matrix()
print("K:", K.shape, "\n", K.toarray().round(4))

