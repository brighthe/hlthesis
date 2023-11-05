import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from top_levelset import TopLevelSet

class Visualizer:
    def __init__(self):
        pass


    def plot_matrices(self, *args, titles=None, cmap="viridis", annot=None, fmt=".2f", annot_kws={"size": 8}):
        """
        Visualizes matrices using seaborn's heatmap.

        Parameters:
        - args: A series of numpy arrays that represent the matrices to be visualized.
        - titles: List of titles corresponding to the matrices. If provided, should match the number of matrices.
        - cmap: Colormap for the heatmap.
        - annot: Whether to annotate each cell with the numeric value. If set to None, automatic based on matrix size.
        - fmt: String formatting code to use when adding annotations.
        - annot_kws: Keyword arguments for heatmap annotation.
        """
        num_matrices = len(args)
        if titles and len(titles) != num_matrices:
            raise ValueError("Number of titles must match the number of matrices.")
        
        plt.figure(figsize=(7 * num_matrices, 7))
        
        for i, matrix in enumerate(args):
            if annot is None:  # Automatic annotation based on matrix size
                # Don't annotate if any dimension is larger than 10
                annot_matrix = matrix.shape[0] <= 10 and matrix.shape[1] <= 10
            else:
                annot_matrix = annot  # Use provided annotation setting

            plt.subplot(1, num_matrices, i + 1)
            sns.heatmap(matrix, annot=annot_matrix, fmt=fmt, cmap=cmap, annot_kws=annot_kws)
            plt.title(titles[i] if titles else f'Matrix {i + 1}')
        
        plt.tight_layout()
        plt.show()

    def plot_stiffness_matrices(self, KE=None, KTr=None):
        """
        Visualizes the stiffness matrices using seaborn's heatmap. Only the matrices
        provided as arguments will be visualized.

        Parameters:
        - KE: Element stiffness matrix for stress analysis (default is None).
        - KTr: Element stiffness matrix for topology derivative calculations (default is None).
        """
        plt.figure(figsize=(14, 7))
        plot_number = 1
        
        if KE is not None:
            # Plot for Element Stiffness Matrix (KE)
            plt.subplot(1, 2, plot_number)
            sns.heatmap(KE, annot=True, fmt=".2f", cmap="viridis")
            plt.title('Element Stiffness Matrix (KE)')
            plot_number += 1

        if KTr is not None:
            # Plot for Topology Derivative Stiffness Matrix (KTr)
            plt.subplot(1, 2, plot_number)
            sns.heatmap(KTr, annot=True, fmt=".2f", cmap="viridis")
            plt.title('Topology Derivative Stiffness Matrix (KTr)')

        plt.tight_layout()
        plt.show()

    def plot_displacement(self, U, nelx, nely):
        """
        Visualizes the displacement field from the finite element analysis.

        Parameters:
        - U: Displacement vector obtained from FEA.
        - nelx: Number of elements in the x-direction.
        - nely: Number of elements in the y-direction.
        """
        # Reshape the displacement vector into two fields
        Ux = U[::2].reshape((nely+1, nelx+1))
        Uy = U[1::2].reshape((nely+1, nelx+1))

        # Compute the magnitude of displacement
        U_mag = np.sqrt(Ux**2 + Uy**2)

        # Visualize the magnitude of displacement using matplotlib
        plt.figure(figsize=(10, 4))
        plt.imshow(U_mag, cmap='viridis', extent=[0, nelx, 0, nely], origin='lower')
        plt.colorbar()  # Shows the color scale
        plt.title('Magnitude of Displacement')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.show()

    def plot_level_set_function(self, struc, lsf):
        """
        Visualizes the original structure and the corresponding level set function.

        Parameters:
        - struc: A 2D numpy array representing the structure.
        - lsf: A 2D numpy array representing the level set function.
        """
        plt.figure(figsize=(10, 5))

        # Plot the original structure
        plt.subplot(1, 2, 1)
        plt.imshow(struc, cmap="gray", origin='upper')
        plt.colorbar()
        plt.title('Original Structure')

        # Plot the level set function result
        plt.subplot(1, 2, 2)
        plt.imshow(lsf, cmap="RdBu", origin='upper')
        plt.colorbar()
        plt.title('Level Set Function (LSF)')

        plt.tight_layout()
        plt.show()


    def plot_sensitivity(self, original_sens, smoothed_sens, title_suffix=''):
        """
        Visualizes the original and smoothed sensitivity matrices.

        Parameters:
        - original_sens: A 2D numpy array representing the original sensitivity.
        - smoothed_sens: A 2D numpy array representing the smoothed sensitivity.
        - title_suffix: A string suffix for the plot titles to distinguish different plots if necessary.
        """
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.imshow(original_sens, cmap='viridis', aspect='auto')
        plt.colorbar()
        plt.title(f'Original Sensitivity{title_suffix}')

        plt.subplot(1, 2, 2)
        plt.imshow(smoothed_sens, cmap='viridis', aspect='auto')
        plt.colorbar()
        plt.title(f'Smoothed Sensitivity{title_suffix}')

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    nely = 30
    nelx = 60
    struc = np.ones((nely, nelx))
    toplevelset = TopLevelSet()
    lsf = toplevelset.reinit(struc)
    print(lsf)


# Instantiate the visualizer
    visualizer = Visualizer()

# Plot the stiffness matrices
    visualizer.plot_level_set_function(struc, lsf)
