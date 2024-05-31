import numpy as np

def applySensitivityFilter(ft, x, dc, dv):
    # 灵敏度滤波器
    if (ft['type'] == 1):
        dc = np.matmul(ft['H'], np.multiply(x, dc)/ft['Hs']/np.maximum(1e-3,x))
    # 密度滤波器
    elif (ft['type'] == 2):
        dc = np.matmul(ft['H'], (dc/ft['Hs']))
        dv = np.matmul(ft['H'], (dv/ft['Hs']))
    return dc, dv