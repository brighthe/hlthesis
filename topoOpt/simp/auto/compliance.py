from jax import value_and_grad

import jax
import jax.numpy as jnp

import numpy as np

# 定义全局变量 mesh
mesh = None

# 投影滤波器
def projectionFilter(rho):
    projection = {'isOn': True, 'beta': 4, 'c0': 0.5}
    if(projection['isOn']):
        v1 = jnp.tanh(projection['c0'] * projection['beta'])
        nm = v1 + jnp.tanh(projection['beta'] * (rho - projection['c0']))
        dnm = v1 + jnp.tanh(projection['beta'] * (1. - projection['c0']))
        return nm / dnm
    else:
        return rho

# SIMP 材料插值模型
def materialModel(rho):
    Emin = 1e-3
    Emax = 1
    penal = 3
    E = Emin + (Emax - Emin) * (rho + 0.01) ** penal
    return E

def getK0():
    E = 1.
    nu = 0.3
    k = jnp.array([1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, \
                -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8])
    KE = E / (1 - nu**2) *\
    jnp.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
                [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
                [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
                [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
                [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
                [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
                [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
                [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ])
    return KE

def getMeshStructure():
    global mesh
    nelx, nely = mesh['nelx'], mesh['nely']
    edofMat = jnp.zeros((nelx*nely, 8), dtype=int)
    for elx in range(nelx):
        for ely in range(nely):
            el = ely + elx * nely
            n1 = (nely + 1) * elx + ely
            n2 = (nely + 1) * (elx + 1) + ely
            edofMat = edofMat.at[el, :].set(jnp.array([2*n1+2, 2*n1+3, 2*n2+2,
                                                       2*n2+3, 2*n2, 2*n2+1, 2*n1, 2*n1+1]))
    iK = tuple(jnp.kron(edofMat, jnp.ones((8, 1))).flatten().astype(int))
    jK = tuple(jnp.kron(edofMat, jnp.ones((1, 8))).flatten().astype(int))
    idx = (iK, jK)
    return edofMat, idx

# 组装全局刚度矩阵
def assembleK(E):
    global mesh
    nelx, nely, ndof = mesh['nelx'], mesh['nely'], mesh['ndof']
    K_asm = jnp.zeros((ndof, ndof))
    K0 = getK0()
    print("K0:", K0.shape, "\n", K0.round(4))
    K_elem = (K0.flatten()[jnp.newaxis]).T

    K_elem = (K_elem * E).T.flatten()
    idx = getMeshStructure()[1]
    K_asm = K_asm.at[(idx)].add(K_elem)
    return K_asm

# 直接法求解线性方程组
def solveKuf(K):
    global mesh
    nelx, nely, ndof = mesh['nelx'], mesh['nely'], mesh['ndof']
    force = jnp.zeros((ndof, 1))
    dofs = jnp.arange(ndof)
    fixed = dofs[0:2 * (nely + 1):1]
    free = jnp.setdiff1d(jnp.arange(ndof), fixed)
    force = force.at[2 * (nelx + 1) * (nely + 1) - 2 * nely + 1, 0].set(-1)

    u_free = jax.scipy.linalg.solve(K[free, :][:, free], force[free], check_finite=False)
    u = jnp.zeros((ndof))
    u = u.at[free].set(u_free.reshape(-1))
    return u

def computeCompliance(rho):
    global mesh
    nelx, nely, ndof = mesh['nelx'], mesh['nely'], mesh['ndof']
    force = jnp.zeros((ndof, 1))
    force = force.at[2 * (nelx + 1) * (nely + 1) - 2 * nely + 1, 0].set(-1)
    print("force:", force.shape, "\n", force)

    print("rho:", rho.shape, "\n", rho)
    rho = projectionFilter(rho)
    print("rho_projection:", rho.shape, "\n", rho)
    E = materialModel(rho)
    K = assembleK(E)
    u = solveKuf(K)
    J = jnp.dot(force.T, u)[0]
    
    return J

# 计算体积约束
def computeGlobalVolumeConstraint(rho):
    vf = 0.5
    g = jnp.mean(rho) / vf - 1.
    return g

def computeFilter(mesh, rmin):
    nelx, nely = mesh['nelx'], mesh['nely']
    H = np.zeros((nelx*nely, nelx*nely))

    for i1 in range(nelx):
        for j1 in range(nely):
            e1 = (i1) * nely + j1
            imin = max(i1 - (np.ceil(rmin) - 1), 0.)
            imax = min(i1 + (np.ceil(rmin)), nelx)
            for i2 in range(int(imin), int(imax)):
                jmin = max(j1 - (np.ceil(rmin) - 1), 0.)
                jmax = min(j1 + (np.ceil(rmin)), nely)
                for j2 in range(int(jmin), int(jmax)):
                    e2 = i2*nely+j2;
                    H[e1, e2] = max(0., rmin - np.sqrt((i1-i2)**2 + (j1-j2)**2))

    Hs = np.sum(H, 1)
    return H, Hs

def applySensitivityFilter(ft, x, dc, dv):
    if (ft['type'] == 1):
        dc = np.matmul(ft['H'], np.multiply(x, dc) / ft['Hs'] / np.maximum(1e-3, x))
    elif (ft['type'] == 2):
        dc = np.matmul(ft['H'], (dc / ft['Hs']))
        dv = np.matmul(ft['H'], (dv / ft['Hs']))
    return dc, dv

# 执行计算并打印调试信息
def main():
    global mesh
    nelx, nely = 6, 3
    elemSize = jnp.array([1., 1.])
    mesh = {'nelx': nelx, 'nely': nely, 'elemSize': elemSize,
            'ndof': 2 * (nelx + 1) * (nely + 1), 'numElems': nelx * nely}
    filterRadius = 1.5
    H, Hs = computeFilter(mesh, filterRadius)
    print("H:", H.shape, "Hs:", Hs.shape)

    rho = jnp.ones((nelx * nely))
    
    # 计算材料模型
    E = materialModel(rho)
    print("E:", E.shape, "\n", E.round(4))

    # 组装全局刚度矩阵
    K = assembleK(E)
    print("K:", K.shape, "\n", K.round(4))

    # 计算位移
    u = solveKuf(K)
    print("u:", u.shape, "\n", u)

    # 计算目标函数及其梯度
    J, dJ = value_and_grad(computeCompliance)(rho)
    print("J:", J)
    print("dJ:", dJ.shape, "\n", dJ)

    # 计算体积约束及其梯度
    vc, dvc = value_and_grad(computeGlobalVolumeConstraint)(rho)
    print("vc:", vc)
    print("dvc:", dvc.shape, "\n", dvc)
    vc, dvc = vc.reshape((1, 1)), dvc.reshape((1, -1))
    

    ft = {'type': 1, 'H': H, 'Hs': Hs}
    dJ, dvc = applySensitivityFilter(ft, rho, dJ, dvc)
    print("dJ_filter:", dJ.shape, "\n", dJ)
    print("dvc_filter:", dvc.shape, "\n", dvc)

if __name__ == "__main__":
    main()