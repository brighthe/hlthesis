import numpy as np
import matplotlib.pyplot as plt

from scipy import ndimage
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.signal import convolve2d

from visualizer import Visualizer

class TopLevelSet:

    def __init__(self, nelx: int = 60, nely: int = 30, volReq: float = 0.3, 
                 stepLength: int = 3, numReinit: int = 2, topWeight: int = 2):

        # Structure and Optimization
        self._nelx = nelx
        self._nely = nely
        self._struc = np.ones((nely, nelx))
        self._volReq = volReq

        self._shapeSens = np.zeros((nely, nelx))
        self._topSens = np.zeros((nely, nelx))

        self._stepLength = stepLength
        self._numReinit = numReinit
        self._topWeight = topWeight


    def reinit(self, struc):
        """
        Reinitialize the level set function for a given structure.

        This function extends the input structure by adding a boundary of void cells,
        computes the Euclidean distance to the nearest solid and void cells,
        and computes the level set function which is negative inside the solid phase
        and positive inside the void phase.

        Parameters:
        - struc (numpy.ndarray): A 2D array representing the solid (1) and void (0) cells of the structure.

        Returns:
        numpy.ndarray: A 2D array of the same shape as 'struc', representing the reinitialized level set function.
        """
        strucFull = np.zeros((struc.shape[0] + 2, struc.shape[1] + 2))
        strucFull[1:-1, 1:-1] = struc

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
        vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
        gFull = np.pad(g, ((1,1),(1,1)), mode='constant', constant_values=0)

        dt = 0.1 / np.max(np.abs(v))
        for _ in range(int(10 * stepLength)):
            dpx = np.roll(lsf, shift=(0, -1), axis=(0, 1)) - lsf
            dmx = lsf - np.roll(lsf, shift=(0, 1), axis=(0, 1))
            dpy = np.roll(lsf, shift=(-1, 0), axis=(0, 1)) - lsf
            dmy = lsf - np.roll(lsf, shift=(1, 0), axis=(0, 1))
            
            lsf = lsf - dt * np.minimum(vFull, 0) * np.sqrt( np.minimum(dmx, 0)**2 + np.maximum(dpx, 0)**2 + np.minimum(dmy, 0)**2 + np.maximum(dpy, 0)**2 ) \
                    - dt * np.maximum(vFull, 0) * np.sqrt( np.maximum(dmx, 0)**2 + np.minimum(dpx, 0)**2 + np.maximum(dmy, 0)**2 + np.minimum(dpy, 0)**2 ) \
                    - dt*w*gFull

        strucFULL = (lsf < 0).astype(int)
        struc = strucFULL[1:-1, 1:-1]
        
        return struc, lsf

    def updateStep(self, lsf, shapeSens, topSens, stepLength, topWeight):

        # Smooth the sensitivities
        shapeSens_smoothed = self.smooth_sensitivity(shapeSens)
        topSens_smoothed = self.smooth_sensitivity(topSens)

        # Load bearing pixels must remain solid - Bridge
        key_positions = [0, round((shapeSens_smoothed.shape[1]-1)/2), round((shapeSens_smoothed.shape[1]-1)/2) + 1, -1]

        shapeSens_smoothed[-1, key_positions] = 0
        topSens_smoothed[-1, key_positions] = 0

        struc, lsf = self.evolve(-shapeSens_smoothed, topSens_smoothed*(lsf[1:-1, 1:-1] < 0), lsf, stepLength, topWeight)

        return struc, lsf


    def optimize(self, Num: int = 1):
        '''
        Num : Maximum number of iterations of the optimization algorithm
        '''
        nelx, nely, struc, volReq, stepLength, numReinit, topWeight = self._nelx, self._nely, self._struc, self._volReq, self._stepLength, self._numReinit, self._topWeight
        shapeSens, topSens = self._shapeSens, self._topSens

        KE, KTr, lambda_, mu = self. materialInfo()
        lsf = self.reinit(struc)
        
        objective = np.zeros(Num)

        la = -0.01
        La = 1000
        alpha = 0.9

        visualize = Visualizer()

        for iterNum in range(Num):
            U = self.FE(nelx, nely, KE, struc)

            for elx in range(nelx):
                for ely in range(nely):
                    n1 = (nely + 1) * elx + ely
                    n2 = (nely + 1) * (elx + 1) + ely
                    Ue = U[np.array([2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n2+2, 2*n2+3, 2*n1+2, 2*n1+3])]

                    shapeSens[ely, elx] = -max(struc[ely, elx], 0.0001) * Ue.T @ KE @ Ue
                    
                    coeff = np.pi/2 * (lambda_ + 2*mu) / mu / (lambda_ + mu)
                    UeT_KE_Ue = 4 * mu * Ue.T @ KE @ Ue
                    additional_term = (lambda_ - mu) * Ue.T @ KTr @ Ue
                    topSens[ely, elx] = struc[ely, elx] * coeff * UeT_KE_Ue * (UeT_KE_Ue + additional_term)

            objective[iterNum] = -np.sum(shapeSens)
            volCurr = np.sum(struc) / (nelx * nely)

            if iterNum > 4 and (abs(volCurr-volReq) < 0.005) and np.all( abs(objective[-1]-objective[-6:-1]) < 0.01*abs(objective[-1]) ):
                break

            if iterNum > 0:
                la = la - 1/La * (volCurr - volReq)
                La = alpha * La

            shapeSens = shapeSens - la + 1/La * (volCurr - volReq)
            topSens = topSens + np.pi * ( la - 1/La * (volCurr - volReq) )
            print("topSens:\n", topSens)

            smoothed_shapeSens = self.smooth_sensitivity(sensitivity = shapeSens)
            smoothed_topSens = self.smooth_sensitivity(sensitivity = topSens)
            print("smoothed_topSens:\n", smoothed_topSens)
            visualize.plot_matrices(shapeSens, smoothed_shapeSens, 
                  titles=['Original shapeSens', 'Smoothed shapeSens'],
                  fmt=".1e", annot_kws={"size": 6})
            visualize.plot_matrices(topSens, smoothed_topSens, 
                  titles=['Original topSens', 'Smoothed topSens'],
                  fmt=".1e", annot_kws={"size": 6})


            struc, lsf = self.updateStep(lsf, shapeSens, topSens, stepLength, topWeight)

            if iterNum % numReinit == 0:
                lsf = self.reinit(struc)
        
            print(f'Iter: {iterNum+1}, Compliance.: {objective[iterNum]:.4f}, Volfrac.: {volCurr:.3f}, la: {la:.3f}, La: {La:.3f}')
        plt.tight_layout()
        plt.show()



if __name__ == "__main__":
    tls = TopLevelSet()
    print(tls.optimize())
