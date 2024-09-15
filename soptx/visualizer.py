import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


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
