import numpy as np
import sys
from scipy.sparse.linalg import spsolve, cg

from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, MeshFactory
from fealpy.boundarycondition import DirichletBC #导入边界条件包
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d

def u(p):
    x = p[..., 0]
    y = p[..., 1]
    r = np.ones(p.shape, dtype=np.float_)
    r[..., 0] = np.sin(y)
    r[..., 1] = np.sin(x)
    return r

def rotu(p):
    x = p[..., 0]
    y = p[..., 1]
    r = np.cos(x)-np.cos(y)
    print("r : ", r.shape)
    return r  

def u(p):
    x = p[..., 0]
    y = p[..., 1]
    r = np.ones(p.shape, dtype=np.float_)
    r[..., 0] = y
    r[..., 1] = 2*x
    return r

def dir(p, t):
    x = p[..., 0]
    y = p[..., 1]
    r = np.ones(p.shape, dtype=np.float_)
    r[..., 0] = y
    r[..., 1] = 2*x
    return np.einsum('...ed, ed->...e', r, t)

def rotu(p):
    x = p[..., 0]
    y = p[..., 1]
    r = np.ones(x.shape, dtype=np.float_)
    return r  

N = int(sys.argv[1])
mesh = MeshFactory.boxmesh2d([0, 1, 0, 1], nx=N, ny=N, meshtype='noconvex')
mesh = HalfEdgeMesh2d.from_mesh(mesh)
space = HCurlVirtualElementSpace2d(mesh)

M = space.mass_matrix()
C = space.curl_matrix()
P = space.PI

bc = DirichletBC(space, dir) 

b = space.source_vector(u)
#print(b)
M, b = bc.apply(M, b)

uh = spsolve(M, b)
print("max uh : ", np.max(np.abs(uh)))
print("MAx m : ", np.max(np.abs(M)))
print("max b :", np.max(np.abs(b)))

cell2dof, cell2dofLoc = space.dof.cell_to_dof()
uhK = np.split(uh[cell2dof], cell2dofLoc[1:-1])

f = lambda x : x[0]@x[1]
uhK = np.concatenate(list(map(f, zip(space.PI, uhK)))).reshape(-1, 2)

e = space.L2_error(u, uh)
print("L2 error : ", e)
e = space.curl_error(rotu, uh)
print("curl error : ", e)










