import numpy as np

def computeFilter(mesh, rmin):
    nelx, nely = mesh['nelx'], mesh['nely']
    H = np.zeros((nelx*nely, nelx*nely));

    for i1 in range(nelx):
        for j1 in range(nely):
            e1 = (i1)*nely+j1
            imin = max(i1-(np.ceil(rmin)-1),0.)
            imax = min(i1+(np.ceil(rmin)),nelx)
            for i2 in range(int(imin), int(imax)):
                jmin = max(j1-(np.ceil(rmin)-1),0.)
                jmax = min(j1+(np.ceil(rmin)),nely)
                for j2 in range(int(jmin), int(jmax)):
                    e2 = i2*nely+j2;
                    H[e1, e2] = max(0., rmin-np.sqrt((i1-i2)**2+(j1-j2)**2))

    Hs = np.sum(H,1);
    return H, Hs;

def applySensitivityFilter(ft, x, dc, dv):
    # 灵敏度滤波器
    if (ft['type'] == 1):
        dc = np.matmul(ft['H'], np.multiply(x, dc)/ft['Hs']/np.maximum(1e-3, x))
    # 密度滤波器
    elif (ft['type'] == 2):
        dc = np.matmul(ft['H'], (dc/ft['Hs']))
        dv = np.matmul(ft['H'], (dv/ft['Hs']))
    return dc, dv