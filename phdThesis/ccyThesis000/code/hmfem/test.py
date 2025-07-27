import itertools
import string
import time
import numpy as np
from scipy.special import factorial
from fealpy.mesh import TriangleMesh

def comb(gamma):
    N = np.sum(gamma)
    #array = np.zeros(N, dtype=np.float_)
    #i = 0
    #for k in gamma:
    #    array[i:i+k] = np.arange(1, k+1)
    #    i += k
    #res = np.prod(array/np.arange(1, N+1))
    res = 1/factorial(N)
    for v in gamma:
        res *= factorial(v)
    return res

def custom_next_permutation(arr, compare_function):
    n = len(arr)
    i = n - 2
    while i >= 0 and compare_function(arr[i], arr[i + 1]) >= 0:
        i -= 1
    if i == -1:
        # 如果没有找到降序的元素，说明当前排列已经是最大的排列
        return False
    # 从右向左查找第一个大于arr[i]的元素
    j = n - 1
    while compare_function(arr[j], arr[i]) <= 0:
        j -= 1
    # 交换arr[i]和arr[j]
    arr[i], arr[j] = arr[j], arr[i]
    # 反转arr[i+1:]，使其成为升序
    arr[i + 1:] = arr[i + 1:][::-1]
    return True

def symmetry_array(arr):
    idx0 = np.arange(len(arr.shape))
    ret = np.zeros_like(arr)
    for idx in itertools.permutations(idx0):
        ret += np.transpose(arr, idx) 
    ret /= factorial(len(idx0))  
    return ret

def span_array(arr, alpha):
    """
    arr : (NC, l, d) 
    alpha : (l, )
    """
    N = np.sum(alpha)
    s = string.ascii_lowercase[:N]
    ss = 'i'+',i'.join(s)
    s = ss+'->i'+s

    tup = (s, )
    for i in range(len(alpha)):
        for j in range(alpha[i]):
            tup = tup + (arr[:, i], )
    return np.einsum(*tup, optimize=False)

def symmetry_span_array(arr, alpha):
    M = span_array(arr, alpha)

    N = np.sum(alpha)
    idx = np.arange(N)
    idx1 = []
    for count, value in enumerate(alpha):
        idx1.extend([count] * value)
    ret = np.zeros_like(M)
    count = 0
    while True:
        ret += np.transpose(M, (0, ) + tuple([i+1 for i in idx])) 
        count += 1
        sss = custom_next_permutation(idx, lambda x, y : idx1[x]-idx1[y])
        if not sss:
            # ret *= np.prod(factorial(alpha))/factorial(np.sum(alpha))
            ret /= count
            break
    return ret

    # 调用C++函数
    arr_ptr = (ctypes.c_int * len(arr))(*arr)
    arr1_ptr = (ctypes.c_int * len(arr1))(*arr1)

    permutations_lib.next_permutation(arr_ptr, arr1_ptr, len(arr))
    pass


def test_formula0():
    """
    @brief 测试张量对称化公式
    """
    mesh = TriangleMesh.from_box([0, 1, 0, 1], 2, 2)
    node = mesh.entity('node')
    c = mesh.entity('cell')[0]
    c2e = mesh.ds.cell_to_edge()[0]
    glambda = mesh.grad_lambda()[0]
    print("\033[91m这是红色的文本\033[0m")
    t = np.zeros((3, 2), dtype=np.float_)
    t[0] = node[c[2]]-node[c[1]]
    t[1] = node[c[0]]-node[c[2]]
    t[2] = node[c[1]]-node[c[0]]

    n = mesh.edge_unit_normal()[c2e]
    N = np.einsum('ik, jk->ij', glambda, n)
    tdelta = np.zeros((3, 2, 2), dtype=np.float_)
    tdelta[0, 0] = t[2]
    tdelta[0, 1] = -t[1]
    tdelta[1, 0] = -t[2]
    tdelta[1, 1] = t[0]
    tdelta[2, 0] = t[1]
    tdelta[2, 1] = -t[0]

    for i in range(3):
        t = tdelta[i]
        gamma = np.array([3, 5])
        tg = span_array(t[None, :], gamma)[0]
        stg0 = symmetry_array(tg)            # t^{gamma} 的对称化
        stg1 = symmetry_span_array(t[None, :], gamma)[0] # 快速计算 t^{gamma} 的对称化
        print(np.max(np.abs(stg0-stg1))<1e-15)

        beta = np.array([1, 3, 4])
        sLg = symmetry_span_array(glambda[None, :], beta)[0]
        t0 = time.time()
        invo0 = np.sum(sLg*tg)
        t1 = time.time()
        invo1 = np.sum(sLg*stg1)
        t2 = time.time()

        flag = np.ones(3, dtype=np.bool_)
        flag[i] = False
        if np.all(gamma-beta[flag]>=0):
            invo2 = (-1)**(beta[i])*comb(gamma)/comb(gamma-beta[flag])
        else:
            invo2 = 0
        t3 = time.time()
        print("t3", t3-t2)
        print("t2", t2-t1)
        print("t1", t1-t0)
        print("invo0: ", invo0)
        print("invo1: ", invo1)
        print("invo2: ", invo2)

        ni = np.array([n[i]])
        nn = span_array(ni[None, :], [sum(beta)])[0]
        t0 = time.time()
        invon0 = np.sum(nn*sLg)
        t1 = time.time()
        invon1 = np.prod(N[:, i]**beta)
        t2 = time.time()
        print("inin : ", invon0)
        print("ininini : ", invon1)
        print(t2-t1)
        print(t1-t0)

def test_flatten_symtery_tensor(n=2):
    """
    T : (..., 2, 2, ..., 2) 有 n 个 2, 代表 2 为 n 阶张量
    """
    shape=(-1, )+tuple([2 for i in range(n)])
    T = np.random.rand(10*2**n)

    s = T.reshape(-1, 2**n)
    f = lambda x: np.array([int(ss) for ss in np.binary_repr(x, n)], dtype=np.int_)
    idx = np.array(list(map(f, np.arange(2**n))))
    print(idx)
    flag = np.ones(len(idx), dtype=np.bool_)
    for i in range(n)[1:]:
        flag = flag & (idx[:, i] <= idx[:, i-1])
    print(idx[flag])

def test_apply():
    from scipy.linalg import solve_triangular

    b = np.array([[1, 2, 3], [0, 4, 5], [0, 0, 6]], dtype=np.float_)
    b = np.tile(b, (3, 1, 1))
    def fff(x):
        print('x : ', x)
        #aaa = solve_triangular(x[y], np.identity(3))
        aaa = np.zeros((1, 1, 1))
        return aaa
    np.apply_along_axis(fff, (0, 1), b)




    
#test_formula0()
#test_flatten_symtery_tensor(n=4)
test_apply()
