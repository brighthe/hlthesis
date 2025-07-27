import numpy as np
import sys
from fealpy.mesh import QuadrangleMesh
import matplotlib.pyplot as plt

mesh = QuadrangleMesh.from_box([0, 4, 0, 1], 8, 2)
NC = 16
puh = np.zeros([NC, 2], dtype=np.float_)
puh[:, 0] = np.linspace(float(sys.argv[1]), float(sys.argv[2]), NC, endpoint=True)
puh[:, 1] = np.linspace(float(sys.argv[3]), float(sys.argv[4]), NC, endpoint=True)
uh = np.linspace(float(sys.argv[5]), float(sys.argv[6]), NC, endpoint=True)
cuh = np.linspace(float(sys.argv[7]), float(sys.argv[8]), NC, endpoint=True)
fname = sys.argv[9]

uhs = '$u_h$'

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(plt, cellcolor=puh[:, 0], cmap=plt.cm.jet, linewidths=0, 
        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
plt.title('x component of the real part of ' + uhs)
plt.savefig(fname+'_fem_puhx.svg')

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(plt, cellcolor=puh[:, 1], cmap=plt.cm.jet, linewidths=0, 
        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
plt.title('y component of the real part of ' + uhs)
plt.savefig(fname+'_fem_puhy.svg')

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(plt, cellcolor=uh, cmap=plt.cm.jet, linewidths=0, 
        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
plt.title('magnitude the real part of ' + uhs)
plt.savefig(fname+'_fem_puhxy.svg')

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(plt, cellcolor=cuh, cmap=plt.cm.jet, linewidths=0, 
        showaxis=True, showcolorbar=True, aspect=1, colorbarshrink=0.4)
plt.title('rot of the real part of ' + uhs)
plt.savefig(fname+'_fem_cuh.svg')


