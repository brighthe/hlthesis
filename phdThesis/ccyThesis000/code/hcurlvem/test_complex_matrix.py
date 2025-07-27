import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix, bmat

from mumps import DMumpsContext
from scipy.sparse.linalg import minres, gmres, cg
def Solve(A, b):

    N = A.shape[0]
    A = bmat([[A.real, -A.imag], [A.imag, A.real]])
    b = np.r_[b.real, b.imag]

    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_icntl(7, 3)
    ctx.set_icntl(14, 0)
    ctx.set_centralized_sparse(A)

    ctx.set_rhs(b)
    ctx.run(job=6)
    ctx.destroy() # Cleanup
    '''
    x, _ = minres(A.real, b.real, tol=1e-10)
    #x, _ = gmres(A, b, tol=1e-10)
    '''
    return b[:N] + b[N:]*(1j)

def solve_complex_sparse(A, b):
    # Convert A and b to complex csc_matrix and array types
    A = csc_matrix(A, dtype=np.complex128)
    b = np.array(b, dtype=np.complex128)

    # Create a DMumpsContext object and set options
    ctx = DMumpsContext()
    ctx.set_silent()
    ctx.set_centralized_sparse(A)

    # Set the right-hand side vector
    ctx.set_rhs(b)

    # Set ICNTL(7) to 3 to enable complex MUMPS solver
    ctx.set_icntl(7, 3)

    # Run MUMPS solver and destroy context
    ctx.run(job=6)
    ctx.destroy()

    # Solve the system using minres method
    x, _ = minres(A, b, tol=1e-10)

    return x

# 构造稀疏矩阵
data = np.array([1+1j, 2+2j, 3+3j, 4+4j, 5+5j, 6+6j])
row = np.array([0, 0, 1, 2, 2, 3])
col = np.array([0, 1, 1, 2, 3, 3])
A = csc_matrix((data, (row, col)), shape=(4, 4))

A = A.real
A = A*(1+1j)

# 构造右侧向量
b = np.array([1+1j, 2+2j, 3+3j, 4+4j])

# 求解线性方程组
x = spsolve(A, b)
#x = Solve(A, b) 

print("稀疏矩阵A：\n", A.toarray())
print("右侧向量b：\n", b)
print("解x：\n", x)

