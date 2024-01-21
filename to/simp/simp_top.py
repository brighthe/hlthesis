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
        - volfrac (float): Volume fraction, representing the desired fraction of
        the design space to be occupied by material. Defaults to 0.5.
        - penal (float): Penalization power, controlling the penalization of intermediate
        densities in the SIMP method. Defaults to 3.0.
        - rmin (float): Filter radius (divided by the element size), used to achieve
        mesh-independence in the design. Defaults to 1.5.
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
        k = [1/2 - nu/6, 1/8 + nu/8, -1/4 - nu/12, -1/8 + 3*nu/8,
            -1/4 + nu/12, -1/8 - nu/8, nu/6, 1/8 - 3*nu/8]
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

    def check(self, nelx, nely, rmin, x, dc):
        """
        Apply the mesh-independency filter to modify the element sensitivities.

        Parameters:
            nelx (int): Number of elements in the x direction.
            nely (int): Number of elements in the y direction.
            rmin (float): Filter radius.
            x (np.array): Density distribution matrix of shape (nely, nelx).
            dc (np.array): Original sensitivities matrix of shape (nely, nelx).

        Returns:
            np.array: Filtered sensitivities matrix of shape (nely, nelx).

        Notes:
            The convolution operator is defined such that it decays linearly with 
            distance from the considered element, and is zero outside the filter region.
        """
        # Assert to ensure all dc values are non-positive
        # assert np.all(dc <= 0), "dc should be non-positive (usually negative)"

        # Initialize the modified sensitivities matrix
        dcn = np.zeros((nely, nelx))

        # Loop over all elements in the design region
        for i in range(nelx):
            for j in range(nely):
                sum_val = 0.0

                # Loop over the square region surrounding the element (i,j) 
                # with a side length of twice round(rmin) to find the elements 
                # that lie within the filter radius
                for k in range( max( i - int(rmin), 0 ), min( i + int(rmin) + 1, nelx ) ):
                    for l in range( max( j - int(rmin), 0 ), min( j + int(rmin) + 1, nely ) ):

                        # Calculate convolution operator value for the element (k,l) with respect to (i,j)
                        fac = rmin - np.sqrt((i - k)**2 + (j - l)**2)

                        # Accumulate the convolution sum
                        sum_val += max(0, fac)

                        # Modify the sensitivity for element (i,j) based on the value of the convolution operator
                        dcn[j, i] += max(0, fac) * x[l, k] * dc[l, k]

                # Normalize the modified sensitivity for element (i,j)
                dcn[j, i] /= (x[j, i] * sum_val)

        return dcn

    def OC(self, nelx, nely, volfrac, x, dc):
        """
        Perform topology optimization using the Optimality Criteria method.

        Parameters:
            nelx (int): Number of elements in the x direction.
            nely (int): Number of elements in the y direction.
            volfrac (float): Volume fraction, representing the desired fraction of
            the design space to be occupied by material.
            x (np.array): Density distribution matrix of shape (nely, nelx).
            dc (np.array): Original sensitivities matrix of shape (nely, nelx).

        Returns:
            xnew (np.array)-- updated design variable vector after the optimization step
        """
        # Ensure that the sensitivity values are non-positive as required for the OC method
        assert np.all(dc <= 0), "dc should be non-positive (usually negative)"

        # Initialize Lagrange multipliers for bi-section algorithm
        l1, l2 = 0, 1e5

        # Define the maximum change allowed in the design variables
        move = 0.2  

        # Create a copy of the design variable vector to be updated
        xnew = np.copy(x)  

        # Bi-section loop to find the correct Lagrange multiplier that satisfies the volume constraint
        while (l2 - l1) > 1e-4:
            # Calculate the midpoint of the current Lagrange multiplier interval
            lmid = 0.5 * (l2 + l1)
            # Lower limit move restriction
            tmp0 = x - move
            # Upper limit move restriction
            tmp1 = x + move
            # Design variable update (intermediate step) using OC update scheme
            tmp2 = x * np.sqrt(-dc / lmid)

            tmp3 = np.minimum(tmp1, tmp2)
            tmp4 = np.minimum(1, tmp3) 
            tmp5 = np.maximum(tmp0, tmp4)
            xnew = np.maximum(0.001, tmp5)
            # Check if the current design violates the volume constraint
            if np.sum(xnew) - volfrac * nelx * nely > 0:
                l1 = lmid
            else:
                l2 = lmid

        # Return the updated design variable vector
        return xnew
