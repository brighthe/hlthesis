import numpy as np
import matplotlib.pyplot as plt

from fealpy.geometry import CircleCurve, FoldCurve
from fealpy.mesh import QuadrangleMesh
from fealpy.mesh import PolygonMesh, HalfEdgeMesh2d 

from scipy.spatial import Delaunay

def find_cut_point(phi, p0, p1):

    """ Find cutted point between edge `(p0, p1)` and the curve `phi`
    
    Parameters
    ----------
    phi : function
        This is a Sign distance function.
    p0 : nd.ndarray, Nx2
        p0 is leftpotint of an edge.
    p1 : nd.ndarray, Nx2  
        p1 is rightpoint of an edge.

    Returns
    -------
    cutpoint : numpy.ndarray, Nx2 
        The return value 'cutpoint' of type is float and 'cutpoint' is a
        Intersection between edge '(p0,p1)' and the curve 'phi'
        
    Raises
    ------
        BadException
        'eps' is a very small number,This is done to prevent a
        situation equal to zero and ignore the point.

    """
    cutPoint = (p0+p1)/2.0
    phi0 = phi(p0)
    phi1 = phi(p1)
    phic = phi(cutPoint)

    isLeft = np.zeros(p0.shape[0], dtype=np.bool_)
    isRight = np.zeros(p0.shape[0], dtype=np.bool_)
    vec = p1 - p0
    h = np.sqrt(np.sum(vec**2, axis=1))

    eps = np.finfo(p0.dtype).eps
    tol = np.sqrt(eps)*h*h
    isNotOK = (h > tol) & (phic != 0)
    while np.any(isNotOK):
        cutPoint[isNotOK, :] = (p0[isNotOK, :] + p1[isNotOK,:])/2
        phic[isNotOK] = phi(cutPoint[isNotOK, :])
        isLeft[isNotOK] = phi0[isNotOK] * phic[isNotOK] > 0
        isRight[isNotOK] = phi1[isNotOK] * phic[isNotOK] > 0
        p0[isLeft, :] = cutPoint[isLeft, :]
        p1[isRight, :] = cutPoint[isRight, :]

        phi0[isLeft] = phic[isLeft]
        phi1[isRight] = phic[isRight]
        h[isNotOK] /= 2
        isNotOK[isNotOK] = (h[isNotOK] > tol[isNotOK]) & (phic[isNotOK] != 0)
        isLeft[:] = False
        isRight[:] = False 
    return cutPoint

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

class HalfEdgeMesh2dWithInterface(HalfEdgeMesh2d):
    def __init__(self, box, interface, nx=10, ny=10):
        backmesh = boxmesh2d(box, nx = nx, ny = ny)
        mesh = HalfEdgeMesh2d.from_mesh(backmesh)

        super(HalfEdgeMesh2dWithInterface, self).__init__(mesh.entity('node')[:], 
                mesh.ds.halfedge[:], mesh.ds.subdomain[:])
        
        self.interface = interface
        self.number_cut_cell = 0
        self.generate_mesh()

    def generate_mesh(self):
        NN = self.number_of_nodes()
        NC = self.number_of_all_cells()
        NE = self.number_of_edges()

        node = self.node
        cell = self.entity('cell')
        halfedge = self.entity('halfedge')
        cstart = self.ds.cellstart

        isMainHEdge = self.ds.main_halfedge_flag()

        phiValue = self.interface(node[:])
        #phiValue[np.abs(phiValue) < 0.1*h**2] = 0.0
        phiSign = np.sign(phiValue)

        # Step 1: find the nodes near interface
        edge = self.entity('edge')
        isCutHEdge = phiValue[halfedge[:, 0]]*phiValue[halfedge[halfedge[:, 4], 0]] < 0 

        cutHEdge, = np.where(isCutHEdge&isMainHEdge)
        cutEdge = self.ds.halfedge_to_edge(cutHEdge)

        e0 = node[edge[cutEdge, 0]]
        e1 = node[edge[cutEdge, 1]]
        cutNode = find_cut_point(self.interface, e0, e1)

        print("aaa", type(node))
        self.refine_halfedge(isCutHEdge, newnode = cutNode)

        newHE = np.where((halfedge[:, 0] >= NN)&(halfedge[:, 1]>=cstart))[0]
        ## 注意这里要把 newHE 区分为界面内部和界面外部的
        cen = self.entity_barycenter('halfedge', index=newHE)
        isinnewHE = self.interface(cen)<0
        newHEin = newHE[isinnewHE]
        newHEout = newHE[~isinnewHE]

        idx = np.argsort(halfedge[newHEout, 1])
        newHE[::2] = newHEout[idx]
        idx = np.argsort(halfedge[newHEin, 1])
        newHE[1::2] = newHEin[idx]
        newHE = newHE.reshape(-1, 2)
        ################################################

        ne = len(newHE)
        NE = len(halfedge)//2

        self.number_cut_cell = ne

        halfedgeNew = halfedge.increase_size(ne*2)
        halfedgeNew[:ne, 0] = halfedge[newHE[:, 1], 0]
        halfedgeNew[:ne, 1] = halfedge[newHE[:, 1], 1]
        halfedgeNew[:ne, 2] = halfedge[newHE[:, 1], 2]
        halfedgeNew[:ne, 3] = newHE[:, 0] 
        halfedgeNew[:ne, 4] = np.arange(NE*2+ne, NE*2+ne*2)

        halfedgeNew[ne:, 0] = halfedge[newHE[:, 0], 0]
        halfedgeNew[ne:, 1] = np.arange(NC, NC+ne)
        halfedgeNew[ne:, 2] = halfedge[newHE[:, 0], 2]
        halfedgeNew[ne:, 3] = newHE[:, 1] 
        halfedgeNew[ne:, 4] = np.arange(NE*2, NE*2+ne) 

        halfedge[halfedge[newHE[:, 0], 2], 3] = np.arange(NE*2+ne, NE*2+ne*2)
        halfedge[halfedge[newHE[:, 1], 2], 3] = np.arange(NE*2, NE*2+ne)
        halfedge[newHE[:, 0], 2] = np.arange(NE*2, NE*2+ne)
        halfedge[newHE[:, 1], 2] = np.arange(NE*2+ne, NE*2+ne*2)

        isNotOK = np.ones(ne, dtype=np.bool_)
        current = np.arange(NE*2+ne, NE*2+ne*2)
        while np.any(isNotOK):
            halfedge[current[isNotOK], 1] = np.arange(NC, NC+ne)[isNotOK]
            current[isNotOK] = halfedge[current[isNotOK], 2]
            isNotOK = current != np.arange(NE*2+ne, NE*2+ne*2)


        #增加主半边
        self.ds.hedge.extend(np.arange(NE*2, NE*2+ne))

        #更新subdomain
        subdomainNew = self.ds.subdomain.increase_size(ne)
        subdomainNew[:] = self.ds.subdomain[halfedge[newHE[:, 0], 1]]

        #更新起始边
        self.ds.hcell.increase_size(ne)
        self.ds.hcell[halfedge[:, 1]] = np.arange(len(halfedge)) # 的编号

        self.ds.NN = self.node.size
        self.ds.NC += ne 
        self.ds.NE = halfedge.size//2

        self.init_level_info()


class InterfaceMesher():
    def __init__(self, box, interface, nx=5, ny=5, level=5):

        self.box = box 
        self.interface = interface
        self.top_nx = nx
        self.top_ny = ny
        nx0 = (2**level)*nx
        ny0 = (2**level)*ny

        self.level0 = level
        self.mesh0 = HalfEdgeMesh2dWithInterface(self.box, interface, nx0, ny0)
        self.cidx = np.zeros(self.mesh0.number_of_cells(), dtype=np.int_)

    def get_mesh(self, level):
        nx = self.top_nx*(2**level)
        ny = self.top_ny*(2**level)
        mesh = HalfEdgeMesh2dWithInterface(self.box, self.interface, nx, ny) 

        cstart = mesh.ds.cellstart

        level0 = self.level0
        mesh0 = self.mesh0
        cidx = self.cidx

        N = 2**(level0-level)
        X = np.arange(self.top_nx*(2**level0))//N
        Y = np.arange(self.top_ny*(2**level0))//N

        cidx[:self.top_nx*self.top_ny*(4**level0)] = (X[:, None]*ny + Y[None, :]).reshape(-1)

        halfedge  = mesh.entity('halfedge')
        halfedge0 = mesh0.entity('halfedge')

        ncc0 = mesh0.number_cut_cell
        ncc = mesh.number_cut_cell
        cidx[-ncc0:] = cidx[halfedge0[-2*ncc0:-ncc0, 1]-cstart]

        NC  = mesh.number_of_cells()
        NC0 = mesh0.number_of_cells()

        iscutcell = np.zeros(NC, dtype=np.bool_) # 小网格上的cutcell
        iscutcell[halfedge[-2*ncc:, 1]-cstart] = True

        cutcellmap = np.zeros(nx*ny, dtype=np.int_)
        cutcellmap[halfedge[-2*ncc:-ncc, 1]-cstart] = halfedge[-ncc:, 1]-cstart

        cutcell0 = np.where(iscutcell[cidx])[0]
        cutcellbary0 = mesh0.entity_barycenter('cell')[cutcell0]
        flag = self.interface(cutcellbary0)<0

        cidx[cutcell0[flag]] = cutcellmap[cidx[cutcell0[flag]]]

        self.cidx = cidx 
        return mesh




