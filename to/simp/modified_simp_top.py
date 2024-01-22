import numpy as np
from numpy.lib import fix

from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

from fealpy.mesh import QuadrangleMesh

class TopModifiedSimp:
    def __init__(self, nelx: int = 60, nely: int = 20, volfrac: float = 0.5, penal: float = 3.0,
                rmin: float = 1.5, ft: int = 1):
        '''
        Parameters:
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 20.
        - volfrac (float): Volume fraction, representing the desired fraction of
        the design space to be occupied by material. Defaults to 0.5.
        - penal (float): Penalization power, controlling the penalization of intermediate
        densities in the SIMP method. Defaults to 3.0.
        - rmin (float): Filter radius (divided by the element size), used to achieve
        mesh-independence in the design. Defaults to 1.5.
        - ft (1 or 2):
        '''
        self._nelx = nelx
        self._nely = nely
        self._volfrac = volfrac
        self._penal = penal
        self._rmin = rmin
        self._ft = ft
        self._mesh = QuadrangleMesh.from_box(box = [0, nelx+1, 0, nely+1], nx = nelx, ny = nely)

    def FE(self, nelx, nely, nu):
        """
        Finite Element Analysis in Python for Topology Optimization.

        Args:
        nelx (int): Number of elements along the x-axis.
        nely (int): Number of elements along the y-axis.
        nu (float): Poisson's ratio.

        Returns:
        U (numpy array): Displacement vector.
        """

        # Material properties and stiffness matrix
        A11 = np.array([[12, 3, -6, -3], [3, 12, 3, 0], [-6, 3, 12, -3], [-3, 0, -3, 12]])
        A12 = np.array([[-6, -3, 0, 3], [-3, -6, -3, -6], [0, -3, -6, 3], [3, -6, 3, -6]])
        B11 = np.array([[-4, 3, -2, 9], [3, -4, -9, 4], [-2, -9, -4, -3], [9, 4, -3, -4]])
        B12 = np.array([[2, -3, 4, -9], [-3, 2, 9, -2], [4, 9, 2, 3], [-9, -2, 3, 2]])

        KE = 1/(1-nu**2)/24 * (np.block([[A11, A12], [A12.T, A11]]) +
             nu * np.block([[B11, B12], [B12.T, B11]]))
        print("KE:", KE.shape, "\n", KE.round(4))



        #from fealpy.functionspace import LagrangeFESpace
        #GD = 2
        #space = LagrangeFESpace(self._mesh, p=1, doforder='vdims')
        #vspace = GD * (space, )
        #cell2dof = vspace[0].cell_to_dof()
        #print("cell2dof:", cell2dof.shape, "\n", cell2dof)
        # Node numbers and DOF indices
        #nodenrs = np.arange(1, (1+nelx)*(1+nely)+1).reshape((1+nely, 1+nelx), order='F')
        nodenrs = np.arange(0, (nelx+1)*(nely+1)).reshape((nely+1, nelx+1), order='F')
        print("nodenrs:", nodenrs.shape, "\n", nodenrs)
        #edofVec = np.reshape(2*nodenrs[:-1, :-1]+1, nelx*nely, order='F')
        edofVec = np.reshape(2*nodenrs[:-1, :-1]+1, nelx*nely, order='F')
        print("edofVec:", edofVec.shape, "\n", edofVec)
        edofMat = np.tile(edofVec, (8, 1)).T
        offsets = np.array([0, 1, 2*nely+2, 2*nely+3, 2*nely, 2*nely+1, -2, -1])
        offsets_mat = np.tile(offsets, (nelx*nely, 1))
        edofMat += offsets_mat
        print("edofMat:", edofMat.shape, "\n", edofMat)

        # Global stiffness matrix assembly
        iK = np.kron(edofMat, np.ones((8, 1))).flatten()
        print("iK:", iK.shape, "\n", iK)
        jK = np.kron(edofMat, np.ones((1, 8))).flatten()
        print("jK:", jK.shape, "\n", jK)
        
        # Define loads and supports (Half MBB-Beam)
        F = np.zeros((2*(nely+1)*(nelx+1), 1))
        F[2, 0] = -1
        fixeddofs = np.union1d(np.arange(1, 2*(nely+1), 2), np.array([2*(nelx+1)*(nely+1)]))
        print("fixeddofs:", fixeddofs.shape, "\n", fixeddofs)
        alldofs = np.arange(2*(nely+1)*(nelx+1))
        freedofs = np.setdiff1d(alldofs, fixeddofs)

        return None
