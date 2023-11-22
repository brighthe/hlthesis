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

if __name__ == '__main__':
    nelx = 60
    nely = 30
    mesh = QuadrangleMesh.from_box(box = [0, nelx, 0, nely], nx = nelx, ny = nely)
    node = mesh.entity('node') # 按列增加
    print("node:", node)
    # 网格中点的 x 坐标
    X = node[:, 0].reshape(nelx+1, nely+1).T
    print(X.shape)
    print("X:", X)
    # 网格中点的 y 坐标
    Y = node[:, 1].reshape(nelx+1, nely+1).T
    print(Y.shape)
    print("Y:", Y)

    Phi = initializeLevelSet(nelx = nelx, nely = nely, X = X, Y = Y)
    print(Phi.shape)
    print("Initialized Level Set Function:\n", Phi)

    output = './mesh/'
    if not os.path.exists(output):
        os.makedirs(output)
    fname = os.path.join(output, 'quad_mesh.vtu')
    mesh.nodedata['phi'] = Phi.flatten('F') # 按列增加
    mesh.to_vtk(fname=fname)

