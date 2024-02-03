import numpy as np

from lsf_top import TopLsf

nelx = 32
nely = 20
volReq = 0.4
stepLength = 2;
topWeight = 2;
numReinit = 3
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq, stepLength=stepLength, topWeight=topWeight, numReinit=3)

# 初始化优化参数
nelx, nely, volReq = ts._nelx, ts._nely, ts._volReq
#mesh = ts._mesh
mesh = ts._mesh_top2
#mesh2 = ts._mesh2
#mesh3 = ts._mesh_top2

node = mesh.entity('node') # 按列增加
#node2 = mesh2.entity('node') # 按列增加
#sorted_indices = np.lexsort((-node2[:, 1], node2[:, 0]))
#node3 = node2[sorted_indices]
#NN = mesh2.number_of_nodes()
#cell = mesh2.entity('cell')
#mesh2.ds.reinit(NN=NN, cell=cell)


#import matplotlib.pyplot as plt
#fig = plt.figure()
#axes = fig.gca()
#mesh.add_plot(axes)
#mesh.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
#mesh.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
#
#fig = plt.figure()
#axes = fig.gca()
#mesh2.add_plot(axes)
#mesh2.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
#mesh2.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
#
#fig = plt.figure()
#axes = fig.gca()
#mesh3.add_plot(axes)
#mesh3.find_node(axes, showindex=True, fontsize=12, fontcolor='r')
#mesh3.find_cell(axes, showindex=True, fontsize=12, fontcolor='b')
#plt.show()



cell = mesh.entity('cell') # 左下角逆时针
print("node:", node.shape, "\n", node)
#print("node2:", node2.shape, "\n", node2)
#print("node3:", node3.shape, "\n", node3)
print("cell:", cell.shape, "\n", cell)

#def reorder_nodes(node_data):
#    # The sorting is done first by the second column (y values) in descending order
#    # and then by the first column (x values) in ascending order
#    sorted_indices = np.lexsort((-node_data[:, 1], node_data[:, 0]))
#    return node_data[sorted_indices]



# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))
#print("struc:", struc.shape, "\n", struc)

# 初始化水平集函数
lsf = ts.reinit(struc = struc)
#print("lsf0:", lsf.shape, "\n", lsf.round(4))

# 初始化灵敏度
shapeSens = np.zeros((nely, nelx))
topSens = np.zeros((nely, nelx))

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
E0 = 1.0
nu = 0.3
integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, struc=struc)
KE = integrator.stiff_matrix()
KTr = integrator.trace_matrix()
lambda_, mu = integrator.lame()

# 优化循环
num = 6
# 初始化 compliance objective value
objective = np.zeros(num)
for iterNum in range(num):
    # 计算全局位移和局部单元位移
    print("iterNum:", iterNum)
    print("struc:\n", struc.round(4))
    U, Ue = ts.FE(mesh=mesh, struc=struc)
    print("U:\n", U.round(4))

    # 计算每个单元的柔度的形状灵敏度
    temp1 = -np.maximum(struc, 0.0001)
    #print("temp1:", temp1.shape, "\n", temp1.round(4))
    #temp2 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nely, nelx).T
    temp2 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nelx, nely).T
    #print("temp2:", temp2.shape, "\n", temp2.round(4))
    shapeSens[:] = np.einsum('ij, ij -> ij', temp1, temp2)
    #print("shapeSens1:", shapeSens.shape, "\n", shapeSens.round(4))

    # 计算每个单元的柔度的拓扑灵敏度
    coef = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
    #temp3 = (4 * mu) * np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nely, nelx).T
    #temp4 = (lambda_ - mu) * np.einsum('ij, jk, ki -> i', Ue, KTr, Ue.T).reshape(nely, nelx).T
    temp3 = (4 * mu) * np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nelx, nely).T
    temp4 = (lambda_ - mu) * np.einsum('ij, jk, ki -> i', Ue, KTr, Ue.T).reshape(nelx, nely).T
    topSens[:] = np.einsum('ij, ij -> ij', coef*struc, (temp3+temp4))
    #print("topSens1:", topSens.shape, "\n", topSens.round(4))

    #if iterNum == 3:
    #    print("U:", U.round(4))
        #print("shapeSens:", shapeSens.round(4))
    # 存储当前迭代的 compliance objective
    objective[iterNum] = -np.sum(shapeSens)

    # 计算当前的 volume fraction
    volCurr = np.sum(struc) / (nelx*nely)

    # 打印当前迭代的结果
    #print(f'Iter: {iterNum}, \
    #      Objective Change: {np.all(abs(objective[-1]-objective[-6:]) < 0.01*abs(objective[-1]))}, \
    #      Volume Change: {abs(volCurr-volReq) < 0.005}')
    #print("struc:", np.sum(struc))

    print(f'Iter: {iterNum}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}')

    # 绘制结果图
    import matplotlib.pyplot as plt
    plt.imshow(-struc, cmap='gray', vmin=-1, vmax=0)
    plt.axis('off')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-5)

    # 五次迭代后执行收敛性检查
    if iterNum > 5 and (abs(volCurr-volReq) < 0.005) and \
        np.all( abs(objective[-1]-objective[-6:-1]) < 0.01*abs(objective[-1]) ):
        break

    # 设置  augmented Lagrangian parameters
    if iterNum == 0:
        la = -0.01
        La = 1000
        alpha = 0.9
    else:
        # TODO 与理论不一致
        la = la - 1/La * (volCurr - volReq)
        La = alpha * La

    # Update the sensitivities with augmented Lagrangian terms
    #print("test2:", volCurr-volReq)
    #print("test2:", -la + 1/La * (volCurr-volReq))
    shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
    #print("shapeSens2:", shapeSens.shape, "\n", shapeSens.round(4))
    topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )
    #print("topSens2:", topSens.shape, "\n", topSens.round(4))

    # Smooth the sensitivities
    shapeSens = ts.smooth_sens(sens=shapeSens)
    #print("shapeSens3:", shapeSens.shape, "\n", shapeSens.round(4))
    topSens = ts.smooth_sens(sens=topSens)
    #print("topSens3:", topSens.shape, "\n", topSens.round(4))

    # 执行设计更新
    struc, lsf = ts.updateStep(lsf=lsf, shapeSens=shapeSens, topSens=topSens, stepLength=stepLength, topWeight=topWeight)
    #print("struc:", struc.shape, "\n", struc)

    # Reinitialize the level set function at specified iterations
    if (iterNum+1) % numReinit == 0:
        #print("lsf:\n", lsf.round(4))
        lsf = ts.reinit(struc)
        #print("struc2:\n", struc.round(4))
        #print("lsf:\n", lsf.round(4))
    #if iterNum == 2:
    #    print("struc:\n", struc.round(4))

    #print("shapeSens:", shapeSens.shape, "\n", shapeSens.round(4))
    #print("topSens:", topSens.shape, "\n", topSens.round(4))
    #print("lsf:\n", lsf.round(4))


plt.ioff()
plt.show()






# 可视化
import os
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'lsf_quad_mesh.vtu')
#mesh.celldata['strucFull'] = strucFull.flatten('F') # 按列增加
#mesh.celldata['lsf'] = lsf.flatten('F') # 按列增加
mesh.to_vtk(fname=fname)




