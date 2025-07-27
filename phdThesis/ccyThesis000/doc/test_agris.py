import sympy as sp
from sympy import Matrix, symbols, diff, simplify, lambdify
import numpy as np

def agris_basis_functions_sympy(vertices):
    """
    使用sympy计算Agris元的基函数（符号计算版本）
    
    参数:
        vertices: 三角形三个顶点的坐标，形状为(3,2)的符号矩阵
    
    返回:
        包含基函数和必要几何量的字典，可以用于后续计算
    """
    # 定义符号变量
    x, y = symbols('x y')
    lambda0, lambda1, lambda2 = symbols('lambda0 lambda1 lambda2')
    
    # 顶点坐标
    x0, y0 = vertices[0, 0], vertices[0, 1]
    x1, y1 = vertices[1, 0], vertices[1, 1]
    x2, y2 = vertices[2, 0], vertices[2, 1]

    v0 = vertices[0]
    v1 = vertices[1]
    v2 = vertices[2]

    # 计算重心坐标
    A = Matrix([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]])
    Ainv = A.inv().T
    lambda0 = Ainv[0, 0] * x + Ainv[0, 1] * y + Ainv[0, 2]
    lambda1 = Ainv[1, 0] * x + Ainv[1, 1] * y + Ainv[1, 2]
    lambda2 = Ainv[2, 0] * x + Ainv[2, 1] * y + Ainv[2, 2]
    lambdas = [lambda0, lambda1, lambda2]

    lam0 = lambdify((x, y), lambda0, 'numpy')
    lam1 = lambdify((x, y), lambda1, 'numpy')
    lam2 = lambdify((x, y), lambda2, 'numpy')
    lam = [lam0, lam1, lam2]

    de = np.zeros([3, 3], dtype=np.float64)
    for i in range(3):
        for j in range(3):
            de[i, j] = lam[i](vertices[j, 0], vertices[j, 1])
    print(de)
    
    # 计算必要的几何量
    t01 = np.array([x1 - x0, y1 - y0])
    t02 = np.array([x2 - x0, y2 - y0])
    t12 = np.array([x2 - x1, y2 - y1])
    
    # 计算边的法向量 (逆时针旋转90度)
    n0 = -np.array([t12[1], -t12[0]])  # 边e0 (x1-x2) 的法向量
    n1 = -np.array([-t02[1], t02[0]])  # 边e1 (x2-x0) 的法向量
    n2 = -np.array([t01[1], -t01[0]])  # 边e2 (x0-x1) 的法向量

    # 计算面积
    area = (t01[0] * t02[1] - t01[1] * t02[0]) / 2
    
    # 计算梯度 (lambda是线性函数，梯度为常数)
    grad_lambda0 = -np.array([t12[1], -t12[0]]) / (2 * area)
    grad_lambda1 = -np.array([-t02[1], t02[0]]) / (2 * area)
    grad_lambda2 = -np.array([t01[1], -t01[0]]) / (2 * area)

    # 构建Λ矩阵
    Lambda = Matrix.zeros(3, 3)
    Lambda[0, 0] = grad_lambda0.dot(n0)
    Lambda[0, 1] = grad_lambda0.dot(n1)
    Lambda[0, 2] = grad_lambda0.dot(n2)
    Lambda[1, 0] = grad_lambda1.dot(n0)
    Lambda[1, 1] = grad_lambda1.dot(n1)
    Lambda[1, 2] = grad_lambda1.dot(n2)
    Lambda[2, 0] = grad_lambda2.dot(n0)
    Lambda[2, 1] = grad_lambda2.dot(n1)
    Lambda[2, 2] = grad_lambda2.dot(n2)
    
    # 定义边基函数
    phi_e0 = (16 / Lambda[0, 0]) * (lambda1**2) * (lambda2**2) * lambda0
    phi_e1 = (16 / Lambda[1, 1]) * (lambda0**2) * (lambda2**2) * lambda1
    phi_e2 = (16 / Lambda[2, 2]) * (lambda0**2) * (lambda1**2) * lambda2
    phi_e = [phi_e0, phi_e1, phi_e2]
    phi_v = []
    
    # 计算顶点v=0的间接函数 (i=1, j=2)
    vstar = [[1, 2], [0, 2], [0, 1]]
    for v in range(3):
        i, j = vstar[v]
        lambda0 = lambdas[v]
        lambda1 = lambdas[i]
        lambda2 = lambdas[j]

        phi_e1 = phi_e[i]
        phi_e2 = phi_e[j]

        #print(Lambda[i, j]/16, i , j)

        #print(v)
        #print("lambda0:", lambda0)
        #print("lambda1:", lambda1)
        #print("lambda2:", lambda2)
    
        tphi0 = sp.Rational(1,2) * (lambda0**3) * (lambda1**2) - (
                sp.Rational(3,32)*Lambda[v,j] + sp.Rational(1,16)*Lambda[i,j]) * phi_e2
        tphi1 = (lambda0**3) * lambda1 * lambda2 - (sp.Rational(1,16)*Lambda[j,j]) * phi_e2 - (
                sp.Rational(1,16)*Lambda[i,i]) * phi_e1
        tphi2 = sp.Rational(1,2) * (lambda0**3) * (lambda2**2) - (
                sp.Rational(3,32)*Lambda[v,i] + sp.Rational(1,16)*Lambda[j,i]) * phi_e1
        tphi3 = (lambda0**4) * lambda1 + 8 * tphi0 + 4 * tphi1 - (
                sp.Rational(1,16)*Lambda[i,i]) * phi_e1 - (sp.Rational(1,4)*Lambda[v,j] + sp.Rational(1,16)*Lambda[i,j]) * phi_e2
        tphi4 = (lambda0**4) * lambda2 + 8 * tphi2 + 4 * tphi1 - (
                sp.Rational(1,16)*Lambda[j,j]) * phi_e2 - (sp.Rational(1,4)*Lambda[v,i] + sp.Rational(1,16)*Lambda[j,i]) * phi_e1
        tphi5 = (lambda0**5) - 20 * tphi0 - 20 * tphi1 - 20 * tphi2 + 5 * tphi3 + 5 * tphi4 - (sp.Rational(5,16)*Lambda[v,i]) * phi_e1 - (sp.Rational(5,16)*Lambda[v,j]) * phi_e2
    
        # 构建系数矩阵 C^{v,2} 和 C^{v,1}
        t_vi = vertices[i] - vertices[v]  # t01
        t_vj = vertices[j] - vertices[v]  # t02
        
        # C^{v,2} 矩阵
        term1 = t_vi[:, np.newaxis] @ t_vi[np.newaxis, :]  # outer product
        term2 = 0.5*(t_vi[:, np.newaxis]@t_vj[np.newaxis, :] + t_vj[:,
                                    np.newaxis]@t_vi[np.newaxis, :])
        term3 = t_vj[:, np.newaxis] @ t_vj[np.newaxis, :]  # outer product
        
        vec_term1 = np.array([term1[0,0], 2*term1[0,1], term1[1,1]])
        vec_term2 = np.array([term2[0,0], 2*term2[0,1], term2[1,1]])
        vec_term3 = np.array([term3[0,0], 2*term3[0,1], term3[1,1]])

        C_v2 = np.array([vec_term1, vec_term2, vec_term3]).T
        
        # C^{v,1} 矩阵
        C_v1 = np.array([t_vi, t_vj]).T
        
        # 计算顶点基函数
        phi0 = tphi0*C_v2[0, 0] + tphi1*C_v2[0, 1] + tphi2*C_v2[0, 2]
        phi1 = tphi0*C_v2[1, 0] + tphi1*C_v2[1, 1] + tphi2*C_v2[1, 2]
        phi2 = tphi0*C_v2[2, 0] + tphi1*C_v2[2, 1] + tphi2*C_v2[2, 2]

        phi3 = tphi3*C_v1[0, 0] + tphi4*C_v1[0, 1]
        phi4 = tphi3*C_v1[1, 0] + tphi4*C_v1[1, 1]
        phi5 = tphi5
        phi_v.append([phi0, phi1, phi2, phi3, phi4, phi5])
    
    # 返回所有基函数和几何量
    return phi_e, phi_v

def compute_dof(phi, vertices):
    x, y = symbols('x y')

    phi_x = diff(phi, x)
    phi_y = diff(phi, y)
    phi_xx = diff(phi_x, x)
    phi_xy = diff(phi_x, y)
    phi_yy = diff(phi_y, y)

    phi   = lambdify((x, y), phi, 'numpy')
    phi_x = lambdify((x, y), phi_x, 'numpy')
    phi_y = lambdify((x, y), phi_y, 'numpy')
    phi_xx = lambdify((x, y), phi_xx, 'numpy')
    phi_xy = lambdify((x, y), phi_xy, 'numpy')
    phi_yy = lambdify((x, y), phi_yy, 'numpy')

    dof = []
    for v in range(3):
        phixxval = phi_xx(vertices[v, 0], vertices[v, 1])
        phixyval = phi_xy(vertices[v, 0], vertices[v, 1])
        phiyyval = phi_yy(vertices[v, 0], vertices[v, 1])
        phixval = phi_x(vertices[v, 0], vertices[v, 1])
        phiyval = phi_y(vertices[v, 0], vertices[v, 1])
        phival = phi(vertices[v, 0], vertices[v, 1])
        dof += [phixxval, phixyval, phiyyval, phixval, phiyval, phival]

    es = [[1, 2], [2, 0], [0, 1]]
    R = np.array([[0, 1], [-1, 0]])
    for e in range(3):
        i, j = es[e]
        me = 0.5*(vertices[i] + vertices[j])
        ne = (vertices[j] - vertices[i])@R
        phixval = phi_x(*me)
        phiyval = phi_y(*me)
        val = phixval * ne[0] + phiyval * ne[1]
        dof += [val]
    for i, val in enumerate(dof):
        dof[i] = float(val)
    return dof

def c1_basis_functions_sympy(vertices):
    """
    使用sympy计算 C1 元的基函数（符号计算版本）
    
    参数:
        vertices: 三角形三个顶点的坐标，形状为(3,2)的符号矩阵
    
    返回:
        包含基函数和必要几何量的字典，可以用于后续计算
    """
    # 定义符号变量
    x, y = symbols('x y')
    lambda0, lambda1, lambda2 = symbols('lambda0 lambda1 lambda2')
    
    # 顶点坐标
    x0, y0 = vertices[0, 0], vertices[0, 1]
    x1, y1 = vertices[1, 0], vertices[1, 1]
    x2, y2 = vertices[2, 0], vertices[2, 1]

    v0 = vertices[0]
    v1 = vertices[1]
    v2 = vertices[2]

    # 计算重心坐标
    A = Matrix([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]])
    Ainv = A.inv().T
    lambda0 = Ainv[0, 0] * x + Ainv[0, 1] * y + Ainv[0, 2]
    lambda1 = Ainv[1, 0] * x + Ainv[1, 1] * y + Ainv[1, 2]
    lambda2 = Ainv[2, 0] * x + Ainv[2, 1] * y + Ainv[2, 2]
    lambdas = [lambda0, lambda1, lambda2]

    lam0 = lambdify((x, y), lambda0, 'numpy')
    lam1 = lambdify((x, y), lambda1, 'numpy')
    lam2 = lambdify((x, y), lambda2, 'numpy')
    lam = [lam0, lam1, lam2]

    de = np.zeros([3, 3], dtype=np.float64)
    for i in range(3):
        for j in range(3):
            de[i, j] = lam[i](vertices[j, 0], vertices[j, 1])
    print(de)
    
    # 计算必要的几何量
    t01 = np.array([x1 - x0, y1 - y0])
    t02 = np.array([x2 - x0, y2 - y0])
    t12 = np.array([x2 - x1, y2 - y1])
    
    # 计算边的法向量 (逆时针旋转90度)
    n0 = -np.array([t12[1], -t12[0]])  # 边e0 (x1-x2) 的法向量
    n1 = -np.array([-t02[1], t02[0]])  # 边e1 (x2-x0) 的法向量
    n2 = -np.array([t01[1], -t01[0]])  # 边e2 (x0-x1) 的法向量

    # 计算面积
    area = (t01[0] * t02[1] - t01[1] * t02[0]) / 2
    
    # 计算梯度 (lambda是线性函数，梯度为常数)
    grad_lambda0 = -np.array([t12[1], -t12[0]]) / (2 * area)
    grad_lambda1 = -np.array([-t02[1], t02[0]]) / (2 * area)
    grad_lambda2 = -np.array([t01[1], -t01[0]]) / (2 * area)

    # 构建Λ矩阵
    Lambda = Matrix.zeros(3, 3)
    Lambda[0, 0] = grad_lambda0.dot(n0)
    Lambda[0, 1] = grad_lambda0.dot(n1)
    Lambda[0, 2] = grad_lambda0.dot(n2)
    Lambda[1, 0] = grad_lambda1.dot(n0)
    Lambda[1, 1] = grad_lambda1.dot(n1)
    Lambda[1, 2] = grad_lambda1.dot(n2)
    Lambda[2, 0] = grad_lambda2.dot(n0)
    Lambda[2, 1] = grad_lambda2.dot(n1)
    Lambda[2, 2] = grad_lambda2.dot(n2)
    
    # 定义边基函数
    phi_e0 = (16 / Lambda[0, 0]) * (lambda1**2) * (lambda2**2) * lambda0
    phi_e1 = (16 / Lambda[1, 1]) * (lambda0**2) * (lambda2**2) * lambda1
    phi_e2 = (16 / Lambda[2, 2]) * (lambda0**2) * (lambda1**2) * lambda2
    phi_e = [phi_e0, phi_e1, phi_e2]
    phi_v = []
    
    # 计算顶点v=0的间接函数 (i=1, j=2)
    vstar = [[1, 2], [0, 2], [0, 1]]
    for v in range(3):
        i, j = vstar[v]
        lambda0 = lambdas[v]
        lambda1 = lambdas[i]
        lambda2 = lambdas[j]

        phi_e1 = phi_e[i]
        phi_e2 = phi_e[j]

        #print(Lambda[i, j]/16, i , j)

        #print(v)
        #print("lambda0:", lambda0)
        #print("lambda1:", lambda1)
        #print("lambda2:", lambda2)
    
        tphi0 = sp.Rational(1,2) * (lambda0**3) * (lambda1**2) - (
                sp.Rational(3,32)*Lambda[v,j] + sp.Rational(1,16)*Lambda[i,j]) * phi_e2
        tphi1 = (lambda0**3) * lambda1 * lambda2 - (sp.Rational(1,16)*Lambda[j,j]) * phi_e2 - (
                sp.Rational(1,16)*Lambda[i,i]) * phi_e1
        tphi2 = sp.Rational(1,2) * (lambda0**3) * (lambda2**2) - (
                sp.Rational(3,32)*Lambda[v,i] + sp.Rational(1,16)*Lambda[j,i]) * phi_e1
        tphi3 = (lambda0**4) * lambda1 + 8 * tphi0 + 4 * tphi1 - (
                sp.Rational(1,16)*Lambda[i,i]) * phi_e1 - (sp.Rational(1,4)*Lambda[v,j] + sp.Rational(1,16)*Lambda[i,j]) * phi_e2
        tphi4 = (lambda0**4) * lambda2 + 8 * tphi2 + 4 * tphi1 - (
                sp.Rational(1,16)*Lambda[j,j]) * phi_e2 - (sp.Rational(1,4)*Lambda[v,i] + sp.Rational(1,16)*Lambda[j,i]) * phi_e1
        tphi5 = (lambda0**5) - 20 * tphi0 - 20 * tphi1 - 20 * tphi2 + 5 * tphi3 + 5 * tphi4 - (sp.Rational(5,16)*Lambda[v,i]) * phi_e1 - (sp.Rational(5,16)*Lambda[v,j]) * phi_e2
    
        # 构建系数矩阵 C^{v,2} 和 C^{v,1}
        t_vi = vertices[i] - vertices[v]  # t01
        t_vj = vertices[j] - vertices[v]  # t02
        
        # C^{v,2} 矩阵
        term1 = t_vi[:, np.newaxis] @ t_vi[np.newaxis, :]  # outer product
        term2 = 0.5*(t_vi[:, np.newaxis]@t_vj[np.newaxis, :] + t_vj[:,
                                    np.newaxis]@t_vi[np.newaxis, :])
        term3 = t_vj[:, np.newaxis] @ t_vj[np.newaxis, :]  # outer product
        
        vec_term1 = np.array([term1[0,0], 2*term1[0,1], term1[1,1]])
        vec_term2 = np.array([term2[0,0], 2*term2[0,1], term2[1,1]])
        vec_term3 = np.array([term3[0,0], 2*term3[0,1], term3[1,1]])

        C_v2 = np.array([vec_term1, vec_term2, vec_term3]).T
        
        # C^{v,1} 矩阵
        C_v1 = np.array([t_vi, t_vj]).T
        
        # 计算顶点基函数
        phi0 = tphi0*C_v2[0, 0] + tphi1*C_v2[0, 1] + tphi2*C_v2[0, 2]
        phi1 = tphi0*C_v2[1, 0] + tphi1*C_v2[1, 1] + tphi2*C_v2[1, 2]
        phi2 = tphi0*C_v2[2, 0] + tphi1*C_v2[2, 1] + tphi2*C_v2[2, 2]

        phi3 = tphi3*C_v1[0, 0] + tphi4*C_v1[0, 1]
        phi4 = tphi3*C_v1[1, 0] + tphi4*C_v1[1, 1]
        phi5 = tphi5
        phi_v.append([phi0, phi1, phi2, phi3, phi4, phi5])
    
    # 返回所有基函数和几何量
    return phi_e, phi_v

def test_agris_basis_functions_sympy():
    """
    测试Agris基函数的计算
    """
    # 定义三角形的顶点坐标
    vertices = np.array([[0, 0], [2.12312414, 1.244512], [0.214124, 4.1251124]])

    phi_e, phi_v = agris_basis_functions_sympy(vertices)

    M = np.zeros([21, 21], dtype=np.float64)
    idx = 0

    print("\n顶点基函数:")
    for i, phi in enumerate(phi_v):
        for j, p in enumerate(phi):
            dof = compute_dof(p, vertices)
            M[idx] = dof
            idx += 1
            print(f"phi_v{i}_{j}", dof)
    
    # 计算基函数

    # 打印结果
    print("边基函数:")
    for i, phi in enumerate(phi_e):
        dof = compute_dof(phi, vertices)
        M[idx] = dof
        idx += 1
        print(f"phi_e{i}", dof)

    # 打印矩阵
    error = np.max(np.abs(M-np.eye(M.shape[0])))
    print("误差:", error)


# 测试代码
if __name__ == "__main__":
    test_agris_basis_functions_sympy()




