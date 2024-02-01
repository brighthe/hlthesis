import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.signal import convolve2d

class TopLevelSet:

    def __init__(self, nelx: int = 60, nely: int = 30, volReq: float = 0.3, 
                stepLength: int = 3, numReinit: int = 2, topWeight: int = 2):
        '''
        初始化拓扑优化问题

        Parameters: 
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 30.
        - volReq : 最终设计所需的体积分数. Defaults to 0.3.
        - stepLength (int): The CFL time step size used in each iteration of the evolution equation. Defaults to 3.
        - numReinit (int): The number of algorithm iterations before performing level set reinitialization. Defaults to 2.
        - topWeight (int): Weight of the topology derivative term in the evolution equation. Defaults to 2.
        '''

        self._nelx = nelx
        self._nely = nely
        self._volReq = volReq
        self._stepLength = stepLength
        self._numReinit = numReinit
        self._topWeight = topWeight


    def reinit(self, strucFull):
        """
        根据给定的结构重置化水平集函数

        This function extends the input structure by adding a boundary of void cells,
        computes the Euclidean distance to the nearest solid and void cells,
        and computes the level set function which is negative inside the solid phase
        and positive inside the void phase.

        Parameters:
        - struc (numpy.ndarray): A 2D array representing the solid (1) and void (0) cells of the structure.

        Returns:
        numpy.ndarray: A 2D array of the same shape as 'struc', 表示重置化后的水平集函数
        """
        # strucFull = np.zeros((struc.shape[0] + 2, struc.shape[1] + 2))
        # strucFull[1:-1, 1:-1] = struc

        # Compute the distance to the nearest void (0-valued) cells.
        dist_to_0 = ndimage.distance_transform_edt(strucFull)

        # Compute the distance to the nearest solid (1-valued) cells.
        dist_to_1 = ndimage.distance_transform_edt(strucFull - 1)

        # Offset the distances by 0.5 to center the level set function on the boundaries.
        temp_1 = dist_to_1 - 0.5
        temp_2 = dist_to_0 - 0.5

        # Calculate the level set function, ensuring the correct sign inside and outside the structure.
        lsf = (~strucFull.astype(bool)).astype(int) * temp_1 - strucFull * temp_2

        return lsf


    def stiffnessMatrix(self, k):
        """
        Constructs an 8x8 elemental stiffness matrix for plane stress problems.

        Parameters:
        - k (numpy.ndarray): An array of material properties used to construct the stiffness matrix.

        Returns:
        - numpy.ndarray: An 8x8 elemental stiffness matrix.
        """
        # Element stiffness matrix symmetry is exploited for efficient assembly
        K = np.array([
            [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
            [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
            [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
            [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
            [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
            [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
            [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
            [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
        ])

        return K


    def materialInfo(self):
        """
        Calculates and returns material property information necessary for 
        stress analysis and topology derivative calculations.

        Returns:
        - KE: Element stiffness matrix for stress analysis.
        - KTr: Element stiffness matrix for topology derivative calculations.
        - lambda_, mu: Lame parameters for the material.
        """
        E = 1.0
        nu = 0.3
        lambda_ = E * nu / ((1 + nu) * (1 - nu))
        mu = E / (2 * (1 + nu))
        k = np.array([1/2 - nu/6,   1/8 + nu/8,   -1/4 - nu/12, -1/8 + 3 * nu/8,
                    -1/4 + nu/12,  -1/8 - nu/8,    nu/6,         1/8 - 3 * nu/8])

        KE = E / (1 - nu**2) * self.stiffnessMatrix(k)

        k = np.array([1/3, 1/4, -1/3, 1/4, -1/6, -1/4, 1/6, -1/4])
        KTr = E / (1 - nu) * self.stiffnessMatrix(k)

        return KE, KTr, lambda_, mu

    def FE(self, nelx, nely, KE, struc):
        """
        Finite element analysis performed for each optimization iteration.

        Parameters:
        - nelx, nely: Number of elements in x and y directions.
        - KE: Element stiffness matrix.
        - struc: Current structure as a binary matrix.

        Returns:
        - U: Displacement vector for the entire structure.
        """
        K = lil_matrix( (2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1)) )
        F = lil_matrix( (2*(nelx+1)*(nely+1), 1) )
        U = np.zeros( 2*(nelx+1)*(nely+1) )

        for elx in range(nelx):
            for ely in range(nely):
                n1 = (nely+1) * elx + ely
                n2 = (nely+1) * (elx+1) + ely
                edof = [2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2, 2*n1+3]

                K[np.ix_(edof, edof)] += max(struc[ely, elx], 0.0001) * KE

        # Define loads and supports (Bridge)
        F[2 * (round(nelx/2)+1) * (nely+1) - 1] = 1
        fixeddofs = np.concatenate( [np.arange( 2*(nely+1)-2, 2*(nely+1) ), 
                                     np.arange( 2*(nelx+1)*(nely+1)-2, 2*(nelx+1)*(nely+1) )] )
        alldofs = np.arange( 2*(nely+1)*(nelx+1) )
        freedofs = np.setdiff1d(alldofs, fixeddofs)

        U[freedofs] = spsolve(csc_matrix(K[np.ix_(freedofs, freedofs)]), F[freedofs])

        return U


    def smooth_sensitivity(self, sensitivity, kernel_size=3, padding_mode='edge'):
        """
        Smooth the sensitivity matrix using convolution with a predefined kernel.

        Parameters:
        - sensitivity (numpy.ndarray): The sensitivity matrix to be smoothed.
        - kernel_size (int): The size of the convolution kernel. Default is 3.
        - padding_mode (str): The mode used for padding. Default is 'edge' which pads with the edge values of the array.

        Returns:
        - numpy.ndarray: The smoothed sensitivity matrix.
        """
        # Define the convolution kernel
        kernel_value = 1 / (2*kernel_size)
        kernel = kernel_value * np.array([[0, 1, 0], 
                                        [1, 2, 1], 
                                        [0, 1, 0]])

        # Apply padding to the sensitivity array
        padded_sensitivity = np.pad(sensitivity, ((1, 1), (1, 1)), mode=padding_mode)

        # Perform the convolution using the padded array and the kernel
        smoothed_sensitivity = convolve2d(padded_sensitivity, kernel, mode='valid')

        return smoothed_sensitivity


    def evolve(self, v, g, lsf, stepLength, w):
        """
        Evolves the level set function to simulate the movement of interfaces over a grid.
        
        Parameters:
        - v (numpy.ndarray): A 2D array representing the speed field (velocity) at each grid point.
        - g (numpy.ndarray): A 2D array representing the force term at each grid point.
        - lsf (numpy.ndarray): A 2D array representing the level set function, which defines the interface.
        - stepLength (float): The total time for which the level set function should be evolved.
        - w (float): A weighting parameter for the influence of the force term on the evolution.

        Returns:
        - struc (numpy.ndarray): The updated structure after evolution, represented as a 2D array.
        - lsf (numpy.ndarray): The evolved level set function as a 2D array.
        """
        # Pad the velocity field and force term with a border of zeros to handle boundary conditions
        vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
        gFull = np.pad(g, ((1,1),(1,1)), mode='constant', constant_values=0)

        # Calculate the time step based on the maximum velocity present and a CFL condition factor
        dt = 0.1 / np.max(np.abs(v))

        # Iteratively update the level set function based on the evolution equation
        for _ in range(int(10 * stepLength)):
            # Compute forward and backward differences in the x and y directions
            dpx = np.roll(lsf, shift=(0, -1), axis=(0, 1)) - lsf # forward differences in x directions
            dmx = lsf - np.roll(lsf, shift=(0, 1), axis=(0, 1)) # backward differences in x directions
            dpy = np.roll(lsf, shift=(-1, 0), axis=(0, 1)) - lsf
            dmy = lsf - np.roll(lsf, shift=(1, 0), axis=(0, 1))
            
            # Use the upwind scheme to update the level set function
            lsf = lsf - dt * np.minimum(vFull, 0) * np.sqrt( np.minimum(dmx, 0)**2 + np.maximum(dpx, 0)**2 + np.minimum(dmy, 0)**2 + np.maximum(dpy, 0)**2 ) \
                    - dt * np.maximum(vFull, 0) * np.sqrt( np.maximum(dmx, 0)**2 + np.minimum(dpx, 0)**2 + np.maximum(dmy, 0)**2 + np.minimum(dpy, 0)**2 ) \
                    - dt*w*gFull

        # Derive the new structure based on the zero level set
        strucFULL = (lsf < 0).astype(int)
        struc = strucFULL[1:-1, 1:-1]
        
        return struc, lsf

    def updateStep(self, lsf, shapeSens, topSens, stepLength, topWeight):
        """
        Performs a design update using shape and topological sensitivity information.
        
        Parameters:
        - lsf (numpy.ndarray): The level set function, which describes the interface of the current structure.
        - shapeSens (numpy.ndarray): The shape sensitivity of the current design.
        - topSens (numpy.ndarray): The topological sensitivity of the current design.
        - stepLength (float): The step length parameter controlling the extent of the evolution.
        - topWeight (float): The weighting factor for the topological sensitivity in the evolution.

        Returns:
        - struc (numpy.ndarray): The updated structure after the design step, represented as a 2D array.
        - lsf (numpy.ndarray): The updated level set function after the design step.
        """
        # Smooth the sensitivities to avoid drastic changes between adjacent elements
        shapeSens_smoothed = self.smooth_sensitivity(shapeSens)
        topSens_smoothed = self.smooth_sensitivity(topSens)

        # Load bearing pixels must remain solid - Bridge
        # Ensure load-bearing elements remain solid by setting their sensitivity to zero
        # The last row's first, middle two, and last elements are made zero
        _, cols = shapeSens.shape
        shapeSens_smoothed[-1, [0, cols//2 - 1, cols//2, -1]] = 0
        topSens_smoothed[-1, [0, cols//2 - 1, cols//2, -1]] = 0

        # Evolve the level set function to find the new structure
        # The negative shape sensitivity is used as the velocity field
        # The topological sensitivity is scaled by the topWeight factor and applied only to the solid part of the structure
        struc, lsf = self.evolve(-shapeSens_smoothed, topSens_smoothed*(lsf[1:-1, 1:-1] < 0), lsf, stepLength, topWeight)

        return struc, lsf


    def optimize(self, Num: int = 200):
        '''
        Perform the topology optimization process.

        Parameters:
        - Num (int): Maximum number of iterations of the optimization algorithm.

        Returns:
        - None
        '''
        # Initialization of parameters and variables
        nelx, nely, volReq, stepLength, numReinit, topWeight = self._nelx, self._nely, self._volReq, self._stepLength, self._numReinit, self._topWeight

        struc = np.ones((nely, nelx))

        shapeSens = np.zeros((nely, nelx))
        topSens = np.zeros((nely, nelx))

        # Retrieve material properties for later calculations
        KE, KTr, lambda_, mu = self. materialInfo()

        # Initialize the level set function based on the initial structure
        lsf = self.reinit(struc)

        # Allocate space for the objective function history
        objective = np.zeros(Num)

        # Initialize the augmented Lagrangian parameters for volume constraint
        la = -0.01 # Lagrange multiplier
        La = 1000 # Lagrange multiplier
        alpha = 0.9 # Reduction rate for the penalty parameter

        # Start the optimization loop
        for iterNum in range(1, Num + 1):
            # Perform finite element analysis and get global displacement vector
            U = self.FE(nelx, nely, KE, struc)

            # Loop over each element to compute sensitivities
            for elx in range(nelx):
                for ely in range(nely):
                    # Global indices of the nodes of the element
                    n1 = (nely + 1) * elx + ely
                    n2 = (nely + 1) * (elx + 1) + ely
                    print("n1:", n1)
                    print("n2:", n2)
                    # Local displacement vector for the element
                    Ue = U[np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2, 2*n1+3])]
                    print("Ue:", Ue.shape, "\n", Ue)

                    # Compute shape sensitivity for compliance
                    shapeSens[ely, elx] = -max(struc[ely, elx], 0.0001) * Ue.T @ KE @ Ue
                    
                    # Compute topological sensitivity for compliance
                    coeff = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
                    UeT_KE_Ue = 4 * mu * Ue.T @ KE @ Ue
                    additional_term = (lambda_ - mu) * Ue.T @ KTr @ Ue
                    topSens[ely, elx] = struc[ely, elx] * coeff * UeT_KE_Ue * (UeT_KE_Ue + additional_term)

            # Store the compliance objective for current iteration
            objective[iterNum - 1] = -np.sum(shapeSens)

            # Calculate the current volume fraction
            volCurr = np.sum(struc) / (nelx * nely)

            # Check for convergence after a certain number of iterations
            if iterNum > 5 and (abs(volCurr-volReq) < 0.005) and np.all( abs(objective[-1]-objective[-6:-1]) < 0.01*abs(objective[-1]) ):
                break

            # Update the augmented Lagrangian parameters for the next iteration
            if iterNum > 1:
                la = la - 1/La * (volCurr - volReq)
                La = alpha * La

            # Update the sensitivities with augmented Lagrangian terms
            shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
            topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )

            # Perform the design update step
            struc, lsf = self.updateStep(lsf, shapeSens, topSens, stepLength, topWeight)

            # Reinitialize the level set function at specified iterations
            if iterNum % numReinit == 0:
                lsf = self.reinit(struc)

            print(f'Iter: {iterNum}, Objective Change: {np.all(abs(objective[-1]-objective[-6:]) < 0.01*abs(objective[-1]))}, Volume Change: {abs(volCurr-volReq) < 0.005}')

            # Print the current iteration's results to the console
            print(f'Iter: {iterNum}, Compliance.: {objective[iterNum - 1]:.4f}, Volfrac.: {volCurr:.3f}, la: {la:.3f}, La: {La:.3f}')

            plt.imshow(-struc, cmap='gray', vmin=-1, vmax=0)
            plt.axis('off')
            plt.axis('equal')
            plt.draw()
            plt.pause(1e-5)

        plt.ioff()
        plt.show()


if __name__ == "__main__":
    from fealpy.mesh import QuadrangleMesh
    nelx = 60
    nely = 30
    mesh = QuadrangleMesh.from_box(box = [0, nelx+2, 0, nely+2], nx = nelx+2, ny = nely+2)
    node = mesh.entity('node') # 按列增加
    print(node.shape)
    print("node:", node)
    # 网格中点的 x 坐标
    # X = node[:, 0].reshape(nelx+1, nely+1).T
    # print(X.shape)
    # print("X:", X)
    # 网格中点的 y 坐标
    # Y = node[:, 1].reshape(nelx+1, nely+1).T
    # print(Y.shape)
    # print("Y:", Y)

    struc = np.ones((nely, nelx))
    print(struc.shape)
    print("struc:\n", struc)
    strucFull = np.zeros((struc.shape[0] + 2, struc.shape[1] + 2))
    strucFull[1:-1, 1:-1] = struc
    print(strucFull.shape)
    print("strucFull:\n", strucFull)

    tls = TopLevelSet()

    lsf = tls.reinit(strucFull = strucFull)
    print(lsf.shape)
    print("lsf:\n", lsf)

    shapeSens = np.zeros((nely, nelx))
    topSens = np.zeros((nely, nelx))
    KE, KTr, lambda_, mu = tls.materialInfo()
    print("KE:\n", KE, "\n", "KTr:\n", KTr, "\n", "lambda_:", lambda_, "mu:", mu)

    # Allocate space for the objective function history
    objective = np.zeros(200)

    volReq = 0.3

    # Initialize the augmented Lagrangian parameters for volume constraint
    la = -0.01 # Lagrange multiplier
    La = 1000 # Lagrange multiplier
    alpha = 0.9 # Reduction rate for the penalty parameter


    # Start the optimization loop
    for iterNum in range(1):
        # Perform finite element analysis and get global displacement vector
        U = tls.FE(nelx, nely, KE, struc)
        
        # Loop over each element to compute sensitivities
        for elx in range(nelx):
            for ely in range(nely):
                # Global indices of the nodes of the element
                n1 = (nely + 1) * elx + ely
                n2 = (nely + 1) * (elx + 1) + ely
                # Local displacement vector for the element
                Ue = U[np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2, 2*n1+3])]

                # Compute shape sensitivity for compliance
                shapeSens[ely, elx] = -max(struc[ely, elx], 0.0001) * Ue.T @ KE @ Ue
                
                # Compute topological sensitivity for compliance
                coeff = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
                UeT_KE_Ue = 4 * mu * Ue.T @ KE @ Ue
                additional_term = (lambda_ - mu) * Ue.T @ KTr @ Ue
                topSens[ely, elx] = struc[ely, elx] * coeff * UeT_KE_Ue * (UeT_KE_Ue + additional_term)

        # print("shapeSens:", shapeSens.shape, "\n", shapeSens)
        # print("topSens:", topSens.shape, "\n", topSens)

        # Store the compliance objective for current iteration
        objective[iterNum] = -np.sum(shapeSens)
        # print("objective:", objective)

        # Calculate the current volume fraction
        volCurr = np.sum(struc) / (nelx * nely)
        # print("volCurr:", volCurr)
        
        # print("iterNum:", iterNum)
        # Check for convergence after a certain number of iterations
        if iterNum > 4 and (abs(volCurr-volReq) < 0.005) and np.all( abs(objective[-1]-objective[-6:-1]) < 0.01*abs(objective[-1]) ):
            break

        # Update the augmented Lagrangian parameters for the next iteration
        if iterNum > 0:
            la = la - 1/La * (volCurr - volReq)
            La = alpha * La
        #print("la:", la)
        #print("La:", La)

        # Update the sensitivities with augmented Lagrangian terms
        shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
        topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )
        print("shapeSens:", shapeSens.shape, "\n", shapeSens)
        print("topSens:", topSens.shape, "\n", topSens)





    import os
    output = './mesh/'
    if not os.path.exists(output):
        os.makedirs(output)
    fname = os.path.join(output, 'quad_mesh_2.vtu')
    mesh.celldata['strucFull'] = strucFull.flatten('F') # 按列增加

    mesh.celldata['lsf'] = lsf.flatten('F') # 按列增加
    mesh.to_vtk(fname=fname)
