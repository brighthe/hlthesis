from jax import value_and_grad

import jax
import jax.numpy as jnp

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
    nelx = 6
    nely = 3
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
    nelx = 6
    nely = 3
    ndof = 2 * (nelx+1) * (nely+1)
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
    nelx = 6
    nely = 3
    ndof = 2 * (nelx+1) * (nely+1)
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
    nelx = 6
    nely = 3
    ndof = 2 * (nelx+1) * (nely+1)
    force = jnp.zeros((ndof, 1))
    force = force.at[2 * (nelx + 1) * (nely + 1) - 2 * nely + 1, 0].set(-1)
    print("force:", force.shape, "\n", force)

    E = materialModel(rho)
    K = assembleK(E)
    u = solveKuf(K)
    J = jnp.dot(force.T, u)[0]
    
    return J

# 执行计算并打印调试信息
def main():
    nelx = 6
    nely = 3
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
    J, gradJ = value_and_grad(computeCompliance)(rho)
    print("J:", J)
    print("gradJ:", gradJ.shape, "\n", gradJ)

if __name__ == "__main__":
    main()