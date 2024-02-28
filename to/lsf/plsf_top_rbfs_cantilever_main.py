from plsf_top_rbfs_cantilever import TopRBFPlsm

# Cantilever
nelx = 6
nely = 3
volfrac = 0.5
ts = TopRBFPlsm(nelx=nelx, nely=nely, volfrac=volfrac)

# 初始化优化参数
nelx, nely, volfrac = ts._nelx, ts._nely, ts._volfrac
mesh = ts._mesh

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
print("cell:", cell.shape, "\n", cell)

Phi = ts.lsf_init(mesh = mesh)
print("Phi:", Phi.shape, "\n", Phi.round(4))

A, G, pGpX, pGpY, Alpha = ts.rbf_init(mesh = mesh, Phi = Phi)
#print("G:", G.shape, "\n", G.round(4))
#print("pGpX:", pGpX.shape, "\n", pGpX.round(4))
#print("pGpY:", pGpY.shape, "\n", pGpY.round(4))
#print("Alpha:", Alpha.shape, "\n", Alpha.round(4))



import matplotlib.pyplot as plt
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
plt.show()
