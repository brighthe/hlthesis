import numpy as np

from fealpy.pde.MaxwellPDE2d import MaxwellPDE2d
from fealpy.decorator import cartesian, barycentric
from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh, MeshFactory, QuadrangleMesh

def boxmesh2d(box, nx=10, ny=10, eps=0.1):
    """

    Notes
    -----
    生成二维矩形区域上的网格，包括结构的三角形、四边形和三角形对偶的多边形网
    格. 
    """
    N = (nx+1)*(ny+1)
    NC = nx*ny
    node = np.zeros((N,2))
    X = np.linspace(box[0], box[1], nx+1)
    Y = np.linspace(box[2], box[3], ny+1)

    if nx > 2:
        X[nx//2+1] = eps
        #X[nx//2] = eps

    X, Y = np.meshgrid(X, Y)
    X = X.T
    Y = Y.T
    node[:, 0] = X.flatten()
    node[:, 1] = Y.flatten()
    NN = len(node)

    idx = np.arange(N).reshape(nx+1, ny+1)
    cell = np.zeros((NC,4), dtype=np.int_)
    cell[:,0] = idx[0:-1, 0:-1].flat
    cell[:,1] = idx[1:, 0:-1].flat
    cell[:,2] = idx[1:, 1:].flat
    cell[:,3] = idx[0:-1, 1:].flat
    return QuadrangleMesh(node, cell)

def plusmesh(box0, box1, nx=10, ny=10):
    mesh0 = MeshFactory.boxmesh2d(box0, nx=nx, ny=ny, meshtype='poly')
    mesh1 = MeshFactory.boxmesh2d(box1, nx=nx, ny=ny, meshtype='poly')

    node0 = mesh0.entity('node')
    cell0, cellloc0 = mesh0.entity('cell')

    node1 = mesh1.entity('node')
    cell1, cellloc1 = mesh1.entity('cell')

    ## 处理 node
    nidx0 = np.where(node0[:, 0]>-1e-16)[0]
    nidx1 = np.where(node1[:, 0]<1e-16)[0]

    nidx0 = nidx0[np.argsort(node0[nidx0, 1])]
    nidx1 = nidx1[np.argsort(node1[nidx1, 1])]

    idxmap = np.zeros(len(node1), dtype=np.int_)
    idxmap[nidx1] = nidx0
    idxmap[node1[:, 0]>1e-16] = np.arange(len(node1)-len(nidx1))+len(node0)

    node = np.r_[node0, node1[node1[:, 0]>1e-6]]

    cell1 = idxmap[cell1]
    cell = np.r_[cell0, cell1]

    cellloc1[1:] += cellloc0[-1]
    cellloc = np.r_[cellloc0, cellloc1[1:]]
    v = node[:, 0]
    #node[:, 0] = np.sign(v)*(np.abs(v))**(1/4)
    return PolygonMesh(node, cell, cellloc)
    #return mesh0 


class PDE0():
    def __init__(self, s):
        self.s = s

    def solution(self, p):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        theta = np.arctan2(y, x)
        theta = (theta >= 0)*theta + (theta < 0)*(theta+2*pi)
        r = np.sqrt(x**2 + y**2)
        val = np.zeros(p.shape, dtype=p.dtype)
        val[..., 0] = s*r**(s-1)*(sin(s*theta)*cos(theta) - cos(s*theta)*sin(theta))
        val[..., 1] = s*r**(s-1)*(sin(s*theta)*sin(theta) + cos(s*theta)*cos(theta))
        return val

    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x**2 + y**2 < 0.25
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val 

    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x**2 + y**2 < 0.25
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val 

    @cartesian
    def curl_solution(self, p):
        return np.zeros_like(p[..., 0])

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        beta = self.beta(p)[..., None]
        flag = x**2 + y**2 < 0.25
        val = -beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        mesh = MeshFactory.boxmesh2d(box, nx=nx, ny=ny, meshtype='poly')
        return HalfEdgeMesh2d.from_mesh(mesh)

class PDE1(MaxwellPDE2d):
    def __init__(self):
        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        super(PDE1, self).__init__(f)

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val

    @cartesian
    def source(self, p):
        x = p[..., 0, None]
        y = p[..., 1, None]

        alpha = self.alpha(p)[..., None]
        beta = self.beta(p)[..., None]

        ccFx = self.curlcurlFx(x, y)
        ccFy = self.curlcurlFy(x, y)
        if type(ccFx) is not np.ndarray:
            ccFx = np.ones(x.shape, dtype=np.float_)*ccFx
        if type(ccFy) is not np.ndarray:
            ccFy = np.ones(x.shape, dtype=np.float_)*ccFy
        ccf = np.c_[ccFx, ccFy] 
        return alpha*ccf - beta*self.solution(p)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        mesh = MeshFactory.boxmesh2d(box, nx=nx, ny=ny, meshtype='tri')
        return HalfEdgeMesh2d.from_mesh(mesh)

class PDE2():
    def __init__(self, s):
        self.s = s

    def solution(self, p):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        theta = np.arctan2(y, x)
        theta = (theta >= 0)*theta + (theta < 0)*(theta+2*pi)
        r = np.sqrt(x**2 + y**2)
        flag0 = x>0
        flag1 = x<0

        theta0 = theta[flag0]
        r0 = r[flag0]
        x0 = x[flag0]

        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag0, 0] = s*x0**(s-1)*sin(s*theta0) - s*x0**s*cos(s*theta0)*sin(theta0)/r0
        val[flag0, 1] = s*x0**s*cos(s*theta0)*cos(theta0)/r0

        theta1 = theta[flag1]
        r1 = r[flag1]
        x1 = np.abs(x[flag1])
        val[flag1, 0] = -s*x1**(s-1)*sin(s*theta1) - s*x1**s*cos(s*theta1)*sin(theta1)/r1
        val[flag1, 1] = s*x1**s*cos(s*theta1)*cos(theta1)/r1
        return val

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val

    @cartesian
    def curl_solution(self, p):
        return np.zeros_like(p[..., 0])

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        beta = self.beta(p)[..., None]
        val = -beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        mesh = MeshFactory.boxmesh2d(box, nx=nx, ny=ny, meshtype='quad')
        return HalfEdgeMesh2d.from_mesh(mesh)

class PDE3():
    def __init__(self, s, eps):
        self.s = s
        self.eps = eps

    def solution(self, p, *args):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0] - self.eps
        y = p[..., 1]

        flag0 = x>0
        flag1 = x<0

        x0 = x[flag0]
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag0, 0] = s*x0**(s-1)

        x1 = np.abs(x[flag1])
        val[flag1, 0] = -s*x1**(s-1)
        return val

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < self.eps
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < self.eps
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val

    @cartesian
    def curl_solution(self, p):
        return np.zeros_like(p[..., 0])

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        beta = self.beta(p)[..., None]
        val = -beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        eps = self.eps
        '''
        if meshtype=='quad':
            mesh = MeshFactory.boxmesh2d([-1, 1, -1, 1], nx=nx*2, ny=ny*2, meshtype='quad')
        elif meshtype=='poly':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx, ny=ny)
        elif meshtype=='quad0':
            mesh = boxmesh2d(box, nx=nx, ny=ny)
        elif meshtype=='poly0':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx*8, ny=ny)
        elif meshtype=='poly1':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx, ny=ny*4)
        '''
        mesh = boxmesh2d(box, nx=nx, ny=ny, eps=eps)
        return HalfEdgeMesh2d.from_mesh(mesh)

class PDE4():
    def __init__(self, s, eps):
        self.s = s
        self.eps = eps

    def solution(self, p, *args):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0] - self.eps
        y = p[..., 1]

        flag0 = x>0
        flag1 = x<0

        x0 = x[flag0]
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag0, 0] = s*x0**(s-1)

        x1 = np.abs(x[flag1])
        val[flag1, 0] = -s*x1**(s-1)
        return val

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < self.eps
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = x < self.eps
        val = np.ones_like(x)
        val[flag] = val[flag]*2
        return val

    @cartesian
    def curl_solution(self, p):
        return np.zeros_like(p[..., 0])

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        beta = self.beta(p)[..., None]
        val = -beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        eps = self.eps
        '''
        if meshtype=='quad':
            mesh = MeshFactory.boxmesh2d([-1, 1, -1, 1], nx=nx*2, ny=ny*2, meshtype='quad')
        elif meshtype=='poly':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx, ny=ny)
        elif meshtype=='quad0':
            mesh = boxmesh2d(box, nx=nx, ny=ny)
        elif meshtype=='poly0':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx*8, ny=ny)
        elif meshtype=='poly1':
            mesh = plusmesh([-1, 0, -1, 1], [0, 1, -1, 1], nx=nx, ny=ny*4)
        '''
        mesh = boxmesh2d(box, nx=nx, ny=ny, eps=eps)
        return HalfEdgeMesh2d.from_mesh(mesh)



