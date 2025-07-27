import ctypes

# 加载共享库
permutations_lib = ctypes.CDLL('./libpermutation.so')  # 指定库的路径

# 声明函数参数类型
permutations_lib.next_permutation.argtypes = [ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int), ctypes.c_int]

# 创建一个示例数组
arr = [0, 1, 2]
arr1 = [0, 0, 1]

# 调用C++函数
arr_ptr = (ctypes.c_int * len(arr))(*arr)
arr1_ptr = (ctypes.c_int * len(arr1))(*arr1)

permutations_lib.next_permutation(arr_ptr, arr1_ptr, len(arr))
print(list(arr_ptr))
permutations_lib.next_permutation(arr_ptr, arr1_ptr, len(arr))
print(list(arr_ptr))
permutations_lib.next_permutation(arr_ptr, arr1_ptr, len(arr))
print(list(arr_ptr))
permutations_lib.next_permutation(arr_ptr, arr1_ptr, len(arr))
print(list(arr_ptr))
