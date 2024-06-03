import numpy as np

from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

from fealpy.mesh import QuadrangleMesh

class TopSimp:
    def __init__(self, nelx: int = 60, nely: int = 20, volfrac: float = 0.5, penal: float = 3.0, rmin: float = 1.5):
        '''
        Parameters:
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 20.
        - volfrac (float): Volume fraction, representing the desired fraction of the design space to be occupied by material. Defaults to 0.5.
        - penal (float): Penalization power, controlling the penalization of intermediate densities in the SIMP method. Defaults to 3.0.
        - rmin (float): Filter radius (divided by the element size), used to achieve mesh-independence in the design. Defaults to 1.5.
        '''
        self._nelx = nelx
        self._nely = nely
        self._volfrac = volfrac
        self._penal = penal
        self._rmin = rmin
        self._mesh = QuadrangleMesh.from_box(box = [0, nelx+1, 0, nely+1], nx = nelx, ny = nely)

    def lk(self):
        """
        Compute the local stiffness matrix for a plane stress problem.

        Returns:
            np.array: 8x8 local stiffness matrix.
        """
        # Young's module and Poisson's ratio
        E, nu = 1.0, 0.3
        k = [1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3*nu/8, -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3*nu/8]
        KE = E / (1 - nu**2) * np.array([
            [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
            [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
            [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
            [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
            [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
            [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
            [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
            [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
        ])
        
        return KE

    def FE(self, nelx, nely, penal, x):
        """
        Compute the global displacement vector by assembling the global stiffness matrix
        and solving the system of equations.

        Parameters:
            nelx (int): Number of elements in the x direction.
            nely (int): Number of elements in the y direction.
            penal (float): Penalization power.
            x (np.array): Density distribution matrix.

        Returns:
            np.array: Global displacement vector.
        """
        # Get local stiffness matrix
        KE = self.lk()

        # Initialize the global stiffness matrix, load matrix and global displacement vector
        K = lil_matrix( ( 2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1) ) )
        F = lil_matrix( ( 2*(nelx+1)*(nely+1), 1) )
        U = np.zeros( 2*(nely+1)*(nelx+1) )
        
        # Assembly of the global stiffness matrix
        for elx in range(nelx):
            for ely in range(nely):
                # Node numbers for the element
                n1 = (nely+1) * elx + ely
                n2 = (nely+1) * (elx+1) + ely

                # DOF mapping
                edof = np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2, 2*n1+3])

                # Insert the local stiffness matrix into the global matrix
                K[np.ix_(edof, edof)] += x[ely, elx] ** penal * KE

        # Define loads and supports (Half MBB-Beam)
        F[1] = -1
        fixeddofs = np.union1d( np.arange(0, 2*(nely+1), 2), np.array([2*(nelx+1)*(nely+1) - 1]) )
        alldofs = np.arange(2 * (nely+1) * (nelx+1))
        freedofs = np.setdiff1d(alldofs, fixeddofs)
        
        # Solve the system of equations
        U[freedofs] = spsolve(csc_matrix(K[np.ix_(freedofs, freedofs)]), F[freedofs])
        U[fixeddofs] = 0
        
        return U

    def optimize(self):
        # Initialize optimization parameterse
        nelx, nely, volfrac = self._nelx, self._nely, self._volfrac
        mesh = self._mesh
        penal = self._penal
        node = mesh.entity('node') # 按列增加
        cell = mesh.entity('cell') # 左下角逆时针
        print("node:", node.shape, "\n", node)
        print("cell:", cell.shape, "\n", cell)
        # Initialize design variable field to the volume fraction
        x = np.full((nely, nelx), volfrac)
        mesh.celldata['x'] = x.flatten('F') # 按列增加

        import os
        output = './mesh/'
        if not os.path.exists(output):
            os.makedirs(output)
        fname = os.path.join(output, 'quad_mesh.vtu')
        mesh.to_vtk(fname=fname)

        KE = self.lk()
        print("KE:", KE.shape, "\n", KE.round(4))

        from fealpy.functionspace import LagrangeFESpace as Space
        from fealpy.fem import LinearElasticityOperatorIntegrator
        from fealpy.fem import BilinearForm
        mesh = self._mesh
        space = Space(mesh, p=1, doforder='vdims')
        vspace = 2*(space, )
        E, nu = 1.0, 0.3
        mu = E / (2 *(1 + nu))
        lambda_ = ( E * nu) / ((1 + nu) * (1 - 2 * nu))
        integrator1 = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=5)
        bform = BilinearForm(vspace)
        bform.add_domain_integrator(integrator1)
        KK = integrator1.assembly_cell_matrix(space=vspace)
        print("KK:", KK.shape, "\n", KK[0].round(4))


if __name__ == "__main__":

    ts = TopSimp()
    ts.optimize()
