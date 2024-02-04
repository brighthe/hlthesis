
import numpy as np
np.set_printoptions(threshold=np.inf)

from lsf_top import TopLsf

nelx = 32
nely = 20
volReq = 0.4
stepLength = 2;
topWeight = 2;
numReinit = 3
ts = TopLsf(nelx=nelx, nely=nely, volReq=volReq, stepLength=stepLength, topWeight=topWeight, numReinit=3)

mesh = ts._mesh_top2

# 定义初始结构为 entirely solid
struc = np.ones((nely, nelx))
struc[0, -1] = 0

#print("struc:\n", struc.round(6))
U, Ue = ts.FE(mesh=mesh, struc=struc)
#print("U", U.shape)

import scipy.io as io
U2_0 = io.loadmat('U2.mat')['U']
#
#error = U.reshape(-1)-U2_all1.reshape(-1)
##print("error:", np.sum(np.abs(error)))
#
error2 = U.reshape(-1)-U2_0.reshape(-1)
print("error:", np.sum(np.abs(error2)))
print(error2)
