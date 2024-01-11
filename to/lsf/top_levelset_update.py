import os
import numpy as np

from fealpy.mesh import QuadrangleMesh

def initializeLevelSet(nelx, nely, X, Y):
    # 初始孔洞的半径
    r = nely * 0.1
    # hX 和 hY 分别是初始孔洞的中心处的 x 和 y 坐标
    hX = nelx * np.concatenate([np.tile([1/6, 5/6], 3), np.tile([0, 1/3, 2/3, 1], 2), [1/2]])
    hY = nely * np.concatenate([np.repeat([0, 1/2, 1], 2), np.repeat([1/4, 3/4], 4), [1/2]])

    # dX 和 dY 分别是所有网格点在 x 方向上和 y 方向上与初始孔洞的中心之间的距离差
    dX = X[:, :, np.newaxis] - hX
    dY = Y[:, :, np.newaxis] - hY

    # 计算网格点到最近孔洞附近的欧氏距离，并限制在 -3 到 3 之间
    Phi = np.sqrt(dX**2 + dY**2) - r
    Phi = np.min(Phi, axis=2)
    Phi = np.clip(Phi, -3, 3)

    return Phi

def initializeRBF(mesh, X, Y, Phi):
    cRBF = 1e-4  # RBF 参数
    nNode = mesh.number_of_nodes() # 节点总数

    # 计算 MQ 样条
    # 计算所有节点间的 X 和 Y 方向距离差
    Ax = np.subtract.outer(X.flatten('F'), X.flatten('F'))
    Ay = np.subtract.outer(Y.flatten('F'), Y.flatten('F'))
    # 构建矩阵 A
    A = np.sqrt(Ax**2 + Ay**2 + cRBF**2)
    print("A:", A.shape)

    # 构建矩阵 G
    P = np.vstack((np.ones(nNode), X.flatten('F'), Y.flatten('F'))).T
    I = np.zeros((3, 3))
    G_upper = np.hstack((A, P))
    G_lower = np.hstack((P.T, I))
    G = np.vstack((G_upper, G_lower))
    print("G:", G.shape)

    # 计算 MQ 样条在 x 方向和 y 方向上的偏导数
    # 构建矩阵 pGpX
    pGpX_upper = np.hstack((Ax / A, np.tile(np.array([0, 1, 0]), (nNode, 1))))
    pGpX_lower = np.hstack((np.tile(np.array([[0], [1], [0]]), (1, nNode)), np.zeros((3, 3))))
    pGpX = np.vstack((pGpX_upper, pGpX_lower))
    print("pGpX:", pGpX.shape)

    # 构建矩阵 pGpY
    pGpY_upper = np.hstack((Ay / A, np.tile(np.array([0, 0, 1]), (nNode, 1))))
    pGpY_lower = np.hstack((np.tile(np.array([[0], [0], [1]]), (1, nNode)), np.zeros((3, 3))))
    pGpY = np.vstack((pGpY_upper, pGpY_lower))
    print("pGpY:", pGpY.shape)

    # 计算 Alpha
    Phi_flat = Phi.flatten('F')
    Alpha = np.linalg.solve(G, np.hstack((Phi_flat, np.zeros(3))))
    print("Alpha:", Alpha)

    return G, pGpX, pGpY, Alpha

if __name__ == '__main__':
    nelx = 60
    nely = 30
    mesh = QuadrangleMesh.from_box(box = [0, nelx, 0, nely], nx = nelx, ny = nely)

    output = './mesh/'
    if not os.path.exists(output):
        os.makedirs(output)
    fname = os.path.join(output, 'quad_mesh.vtu')

    node = mesh.entity('node') # 按列增加
    #print("node:", node)
    # 网格中点的 x 坐标
    X = node[:, 0].reshape(nelx+1, nely+1).T
    #print(X.shape)
    #print("X:", X)
    # 网格中点的 y 坐标
    Y = node[:, 1].reshape(nelx+1, nely+1).T
    #print(Y.shape)
    #print("Y:", Y)

    Phi = initializeLevelSet(nelx = nelx, nely = nely, X = X, Y = Y)
    print("Phi:", Phi.shape)
    print("Initialized Level Set Function:\n", Phi)

    mesh.nodedata['phi'] = Phi.flatten('F') # 按列增加

    G, pGpX, pGpY, Alpha = initializeRBF(mesh = mesh, X = X, Y = Y, Phi = Phi)
    nNode = mesh.number_of_nodes() # 节点总数
    A = G[:nNode, :nNode]

    # 找到中心节点的索引，调整为从 0 开始的索引
    centerNodeIndex = (nely // 2) + (nelx // 2) * (nely + 1)

    # 提取与中心节点相关的 RBF 值
    gi_x = A[centerNodeIndex, :]

    # 将 RBF 值重新整形为网格的形状以便绘制
    gi_x_reshaped = gi_x.reshape((nely+1, nelx+1), order='F')
    print("gi_x_reshaped:", gi_x_reshaped.shape)
    print("MQ Spline at (0, 0):\n", gi_x_reshaped)

    mesh.nodedata['g_i(x)'] = gi_x.flatten('F') # 按列增加

    # 提取与中心节点相关的偏导数
    dgi_dx = pGpX[centerNodeIndex, :]
    dgi_dy = pGpY[centerNodeIndex, :]
    
    # 将 pGpX 和 pGpY 值重新整形为网格的形状以便绘制
    dgi_dx_reshaped = dgi_dx[:-3].reshape((nely+1, nelx+1), order='F')
    dgi_dy_reshaped = dgi_dy[:-3].reshape((nely+1, nelx+1), order='F')
    print("dgi_dx_reshaped:", dgi_dx_reshaped.shape)
    print("The partial derivates of MQ Spline in x direction:", dgi_dx_reshaped)
    print("dgi_dy_reshaped:", dgi_dy_reshaped.shape)
    print("The partial derivates of MQ Spline in y direction:", dgi_dy_reshaped)

    mesh.nodedata['dgi_dx'] = dgi_dx_reshaped.flatten('F') # 按列增加
    mesh.nodedata['dgi_dy'] = dgi_dy_reshaped.flatten('F') # 按列增加


    mesh.to_vtk(fname=fname)

