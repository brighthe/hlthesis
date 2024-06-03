import numpy as np
from TriTop_Class import TriTop

nx = 8
ny = 5
maxedge = 0.2;
minedge = 0.025;
# 矩形区域的边界
BDY = np.array([[-.5 * nx, -.5 * ny], [.5 * nx, .5 * ny]]) / 100
x = np.linspace(BDY[0,0], BDY[1,0], nx+1)
y = np.linspace(BDY[0,1], BDY[1,1], ny+1)
xn, yn = np.meshgrid(x, y)
# 网格节点的密度值
dN = np.sin(xn/BDY[1,0]*6*np.pi) * np.cos(yn/BDY[1,0]*6*np.pi) + 0.5

tt = TriTop()

for iterNum in range(1):
    C = tt.test(xn=xn, yn=yn, dN=dN)
    tt.generate_mesh(xn=xn, yn=yn, dN=dN, maxedge=maxedge, minedge=minedge, BDY=BDY, maxiter=80)