import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from top_levelset import TopLevelSet

class Visualizer:
    def __init__(self):
        pass
    
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

    #def plot_displacement(self, U, nelx, nely):
    #    """
    #    Visualizes the displacement field from the finite element analysis.

    #    Parameters:
    #    - U: Displacement vector obtained from FEA.
    #    - nelx: Number of elements in the x-direction.
    #    - nely: Number of elements in the y-direction.
    #    """
    #    # Reshape the displacement vector into two fields
    #    Ux = U[::2].reshape((nely+1, nelx+1))
    #    Uy = U[1::2].reshape((nely+1, nelx+1))

    #    # Compute the magnitude of displacement
    #    U_mag = np.sqrt(Ux**2 + Uy**2)

    #    # Visualize the magnitude of displacement using matplotlib
    #    plt.figure(figsize=(10, 4))
    #    plt.imshow(U_mag, cmap='viridis', extent=[0, nelx, 0, nely], origin='lower')
    #    plt.colorbar()  # Shows the color scale
    #    plt.title('Magnitude of Displacement')
    #    plt.xlabel('X Coordinate')
    #    plt.ylabel('Y Coordinate')
    #    plt.show()


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
        im = plt.imshow(U_mag, cmap='viridis', aspect='auto', extent=[0, nelx, 0, nely], origin='lower')
        plt.colorbar(im, fraction=0.046, pad=0.04)  # Shows the color scale
        plt.title('Magnitude of Displacement')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.axis('equal')
        plt.axis('tight')
        plt.show()


if __name__ == "__main__":
# 使用示例：
# Ensure reproducibility
    np.random.seed(0)

# Instantiate the class
    material_analysis = TopLevelSet()

# Get the material information
    KE, KTr, lambda_, mu = material_analysis.materialInfo()

# Instantiate the visualizer
    visualizer = Visualizer()

# Plot the stiffness matrices
    visualizer.plot_stiffness_matrices(KE, KTr)
