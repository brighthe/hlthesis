import numpy as np
import sympy as sym
from sympy.vector import CoordSys3D

from fealpy.pde.MaxwellPDE_2d import MaxwellPDE2d
from fealpy.decorator import cartesian, barycentric
from fealpy.mesh import TriangleMesh, PolygonMesh, QuadrangleMesh
from fealpy.mesh.halfedge_mesh import HalfEdgeMesh2d

from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve, DoubleBandY, Polygon, BandY 
from interface_mesh_2d import HalfEdgeMesh2dWithInterface, InterfaceMesher

def boxmesh2d(box, nx=10, ny=10, xlist=[], ylist=[]):
    X = np.linspace(box[0], box[1], nx+1)
    Y = np.linspace(box[2], box[3], ny+1)

    nx += len(xlist)
    ny += len(ylist)
    N = (nx+1)*(ny+1)
    NC = nx*ny
    node = np.zeros((N,2))

    X = np.sort(np.r_[X, np.array(xlist)])
    Y = np.sort(np.r_[Y, np.array(ylist)])

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

a = 1
class PDE0():
    """
    @brief 用于鲁棒性测试
    """
    def __init__(self, s, eps):
        self.s = s
        self.eps = eps

    @cartesian
    def solution(self, p, *args):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]

        val = np.zeros(p.shape, dtype=p.dtype)
        val[..., 0] = (1/s-1)*np.abs(x-self.eps)**s + cos(a*(x+y))
        val[..., 1] = sin(a*(x+y))
        return val

    @cartesian
    def rotrotsolution(self, p, *args):
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]

        val = np.zeros(p.shape, dtype=p.dtype)
        val[..., 1] =  a**2*(sin(a*(x+y)) - cos(a*(x+y)))
        val[..., 0] = - val[..., 1]
        return val

    def interface(self, p):
        return p[..., 0] - self.eps 

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
        s = self.s
        sin = np.sin
        cos = np.cos
        pi = np.pi
        x = p[..., 0]
        y = p[..., 1]
        return a*cos(a*(x+y)) + a*sin(a*(x+y))

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        alpha = self.alpha(p)[..., None]
        beta = self.beta(p)[..., None]
        val = alpha*self.rotrotsolution(p)-beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx=1, ny=1):
        box = [-1, 1, -1, 1]
        eps = self.eps
        mesh = boxmesh2d(box, nx=nx, ny=ny, xlist=[eps])
        #mesh = TriangleMesh.from_quadrangle_mesh(mesh)
        return HalfEdgeMesh2d.from_mesh(mesh)
        return mesh

class PDE1():
    """
    @brief 圆形界面的 maxwell 方程问题, 用于验证收敛性，
        其中的参数中 0: 内部; 1: 外部
    """
    def __init__(self, k1=20, r0=np.pi/5, r1=1, mu0=1, mu1=0.1):
        self.interface = CircleCurve([0.0, 0.0], r0)
        self.r0 = r0

        C = CoordSys3D('C')
        k0 = k1*(r1**2-r0**2)
        u0 = -mu0*k0*(r0**2 - C.x**2 - C.y**2)*C.y*C.i - mu0*k0*(r0**2 - 
                C.x**2 - C.y**2)*C.x*C.j
        u1 = -mu1*k1*(r1**2 - C.x**2 - C.y**2)*(r0**2 - 
                C.x**2 - C.y**2)*C.y*C.i -mu1*k1*(r1**2 - 
                        C.x**2 - C.y**2)*(r0**2 - C.x**2 - C.y**2)*C.x*C.j 

        self.pde0 = MaxwellPDE2d(u0)
        self.pde1 = MaxwellPDE2d(u1)

    def solution(self, p, *args):
        flag = self.interface(p)<0
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag] = self.pde0.solution(p[flag])
        val[~flag] = self.pde1.solution(p[~flag])
        return val

    def rotrotsolution(self, p, *args):
        flag = self.interface(p)<0
        val = np.zeros(p.shape, dtype=p.dtype)
        val[flag] = self.pde0.curl_curl_solution(p[flag])
        val[~flag] = self.pde1.curl_curl_solution(p[~flag])
        return val

    @cartesian
    def alpha(self, p):
        flag = self.interface(p)<0
        val = 10*np.ones(p.shape[:-1], dtype=p.dtype)
        val[flag] = 1
        return val 

    @cartesian
    def beta(self, p):
        flag = self.interface(p)<0
        val = 10*np.ones(p.shape[:-1], dtype=p.dtype)
        val[flag] = 1
        return val 

    @cartesian
    def curl_solution(self, p):
        flag = self.interface(p)<0
        val = np.zeros(p.shape[:-1], dtype=p.dtype)
        val[flag] = self.pde0.curl_solution(p[flag])
        val[~flag] = self.pde1.curl_solution(p[~flag])
        return val

    @cartesian
    def source(self, p):
        alpha = self.alpha(p)[..., None]
        beta = self.beta(p)[..., None]
        val = alpha*self.rotrotsolution(p)-beta*self.solution(p)
        #val =beta*self.solution(p)
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = self.solution(p)
        return np.einsum('...ed, ed->...e', val, t)

    def get_mesh(self, nx, ny):
        #mesh = HalfEdgeMesh2dWithInterface([-1, 1, -1, 1], 
        #        self.interface, nx, ny)
        mesh = HalfEdgeMesh2d.from_interface_cut_box(self.interface, [-0, 1,
            -0, 1], nx, ny)
        #mesh = HalfEdgeMesh2d.from_interface_cut_box_tri(self.interface, [0, 1, 0, 1], nx, ny)
        return mesh
        #mesh = interfacemesh2d([-1.1, 1, -1.1, 1], self.interface, nx)
        #return HalfEdgeMesh2d.from_mesh(mesh)

class PDE2(MaxwellPDE2d):
    def __init__(self, interface):
        self.interface = interface

        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        f = 0.1*sym.sin(2*C.y)*C.i + sym.sin(2*C.x)*C.j
        super(PDE2, self).__init__(f)

    @cartesian
    def alpha(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = self.interface(p) < 0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        x = p[..., 0]
        y = p[..., 1]
        flag = self.interface(p) < 0
        val = 1*np.ones_like(x)
        val[flag] = val[flag]*1
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
        #mesh = HalfEdgeMesh2d.from_interface_cut_box(self.interface, [0, 4, 0,
        #    1], nx*4, ny)
        mesh = PolygonMesh.nonconvex_octagonal_mesh([0, 4, 0, 1], nx*4, ny)
        mesh = HalfEdgeMesh2d.from_mesh(mesh)
        
        #mesh = HalfEdgeMesh2d.from_interface_cut_box_tri(self.interface, [0, 4,
        #    0, 1], nx*4, ny)
        return mesh

class PDE3():
    def __init__(self, om = 5, eps = 0.01, sig0 = 1, sig1=0.1, level = 7,
            interfacetype='double_circle'):
        self.om = om
        self.eps = eps
        self.sig0 = sig0
        self.sig1 = sig1

        if interfacetype=='double_circle':
            self.interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
            self.mesher = InterfaceMesher([0, 4, 0, 1], self.interface, 8, 2, level)
        elif interfacetype=='double_band':
            self.interface = DoubleBandY()
            self.mesher = InterfaceMesher([0, 4, 0, 1], self.interface, 16, 4, level)
        elif interfacetype=='five_band':
            self.interface = BandY([0.11, 0.13, 0.24, 0.26, 0.49, 0.51, 0.74, 0.76, 0.87, 0.89])
            #self.interface = BandY([0.09, 0.11, 0.24, 0.26, 0.49, 0.51, 0.74, 0.76, 0.89, 0.91])
            self.mesher = InterfaceMesher([0, 4, 0, 1], self.interface, 32, 8, level)
        elif interfacetype=='poly':
            points = np.array([[311+2/3, 50], [290+0.625, 47+2/3], [225, 68+1/39],
             [206.25, 19+1/3], [198.05, 66+2/3], [178.125, 15+1/3], [175+1/3, 62+1/3],
             [100+1/3, 69-1/4], [243.75, 82+1/3], [287.5, 52+1/3]])/100
            points[:, 1] = 1-points[:, 1]
            self.interface = Polygon(points)

    @cartesian
    def source(self, p):
        om = self.om
        eps = self.eps

        x = p[..., 0]
        y = p[..., 1]
        val = np.zeros_like(p, dtype=np.complex128)
        val[..., 1] = -om*(1j)*np.exp(-((x-3)**2)/(eps**2))
        #val[..., 1] = -om*(1j)*np.exp(-((y+5*x-15.5)**2/26)/(eps**2))
        return val

    @cartesian
    def alpha(self, p):
        val = np.ones_like(p[..., 0])
        return val 

    @cartesian
    def beta(self, p):
        flag = self.interface(p)<0

        om = self.om
        eps = self.eps
        sig0 = self.sig0
        sig1 = self.sig1

        val = np.zeros_like(p[..., 0], dtype=np.complex128)
        val[flag] = np.ones_like(p[..., 0])[flag]*(om**2*(eps + sig0/om * (1j)))
        val[~flag] = np.ones_like(p[..., 0])[~flag]*(om**2*(eps + sig1/om * (1j)))
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = np.zeros_like(p[..., 0])
        return val 

    def get_deepest_mesh(self):
        return self.mesher.mesh0

    def get_mesh(self, level):
        """
        level 要小于 init 里的 level
        """
        mesh = self.mesher.get_mesh(level)
        return mesh, self.mesher.cidx

    def get_mesh_0(self, nx, ny):
        mesh = HalfEdgeMesh2d.from_interface_cut_box(self.interface, [0, 4, 0, 1], nx, ny)
        #mesh = HalfEdgeMesh2d.from_interface_cut_box_tri(self.interface, [0, 4, 0, 1], nx, ny)
        return mesh

class PDE4(MaxwellPDE2d):
    def __init__(self, level = 7, interfacetype='double_circle'):
        if interfacetype=='double_circle':
            self.interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
            self.mesher = InterfaceMesher([0, 4, 0, 1], self.interface, 8, 2, level)
        elif interfacetype=='double_band':
            self.interface = DoubleBandY()
            self.mesher = InterfaceMesher([0, 4, 0, 1], interface, 16, 4, 5)

        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        super(PDE4, self).__init__(f)

    @cartesian
    def alpha(self, p):
        flag = x < self.interface(p)<0
        val = np.ones_like(x)
        val[flag] = val[flag]*1
        return val 

    @cartesian
    def beta(self, p):
        flag = x < self.interface(p)<0
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

    def get_deepest_mesh(self):
        return self.mesher.mesh0

    def get_mesh(self, level):
        """
        level 要小于 init 里的 level
        """
        mesh = self.mesher.get_mesh(level)
        return mesh, self.mesher.cidx

class PDE5():
    """
    @brief 五个层的情况, 不能计算阶数
    """
    def __init__(self, om = 5, eps = 0.01, sig0 = 1, sig1=0.1,
            ylist = [0.09, 0.11, 0.24, 0.26, 0.49, 0.51, 0.74, 0.76, 0.89, 0.91]):
        self.layers_list = np.array(ylist).reshape(-1, 2)
        self.om = om
        self.eps = eps
        self.sig0 = sig0
        self.sig1 = sig1

    def interface(self, p):
        N = len(self.layers_list)
        dist = np.ones(p.shape[:-1], dtype=np.float_) 
        for i in range(N):
            flag = (p[..., 1]>self.layers_list[i, 0])&(p[..., 1]<self.layers_list[i, 1])
            dist[flag] = -1
        return dist 

    @cartesian
    def source(self, p):
        om = self.om
        eps = self.eps

        x = p[..., 0]
        val = np.zeros_like(p, dtype=np.complex128)
        val[..., 1] = -om*(1j)*np.exp(-((x-3)**2)/(eps**2))
        return val

    @cartesian
    def alpha(self, p):
        val = np.ones_like(p[..., 0])
        return val 

    @cartesian
    def beta(self, p):
        flag = self.interface(p)<0

        om = self.om
        eps = self.eps
        sig0 = self.sig0
        sig1 = self.sig1

        val = np.zeros_like(p[..., 0], dtype=np.complex128)
        val[flag] = np.ones_like(p[..., 0])[flag]*(om**2*(eps + sig0/om * (1j)))
        val[~flag] = np.ones_like(p[..., 0])[~flag]*(om**2*(eps + sig1/om * (1j)))
        return val 

    @cartesian
    def dirichlet(self, p, t):
        val = np.zeros_like(p[..., 0])
        return val 

    def get_deepest_mesh(self):
        return self.mesher.mesh0

    def get_mesh(self, nx, ny, meshtype='quad'):
        """
        level 要小于 init 里的 level
        """
        box = [0, 4, 0, 1]
        mesh = boxmesh2d(box, nx, ny, [], list(self.layers_list.reshape(-1)))
        if meshtype=='tri':
            mesh = TriangleMesh.from_quadrangle_mesh(mesh)
            return mesh
        else:
            return HalfEdgeMesh2d.from_mesh(mesh)

    def get_mesh_0(self, nx, ny):
        mesh = HalfEdgeMesh2d.from_interface_cut_box(self.interface, [0, 4, 0, 1], nx, ny)
        return mesh

class PDE6(MaxwellPDE2d):
    def __init__(self, eps = 1e-7):
        self.epsilon = eps

        C = CoordSys3D('C')
        f = sym.sin(sym.pi*C.y)*C.i + sym.sin(sym.pi*C.x)*C.j
        super(PDE6, self).__init__(f)

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
        val[flag] = val[flag]*2000
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
        eps = self.epsilon
        mesh = boxmesh2d(box, nx=nx, ny=ny, xlist=[eps])
        return HalfEdgeMesh2d.from_mesh(mesh)


