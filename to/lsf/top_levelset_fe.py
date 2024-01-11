import numpy as np
from fealpy.mesh import QuadrangleMesh

class TopLevelset:

    def __init__(self, nelx: int = 60, nely: int = 30, volReq: float = 0.3, 
                stepLength: int = 3, numReinit: int = 2, topWeight: int = 2):
        '''
        Initialize the topological optimization problem.

        Parameters: 
        - nelx (int): Number of elements in the horizontal direction. Defaults to 60.
        - nely (int): Number of elements in the vertical direction. Defaults to 30.
        - volReq (float): The required volume fraction of the final design. Defaults to 0.3.
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
        self._mesh = QuadrangleMesh.from_box(box = [0, nelx+1, 0, nely+1], nx = nelx, ny = nely)

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
        mesh = self._mesh
        node = mesh.entity('node') # 按列增加
        cell = mesh.entity('cell') # 左下角逆时针
        print("node:", node.shape, "\n", node)
        print("cell:", cell.shape, "\n", cell)

        struc = np.ones((nely, nelx))
        shapeSens = np.zeros((nely, nelx))
        topSens = np.zeros((nely, nelx))

        mesh.celldata['struc'] = struc.flatten('F') # 按列增加
        mesh.celldata['shapeSens'] = shapeSens.flatten('F') # 按列增加
        mesh.celldata['topSens'] = topSens.flatten('F') # 按列增加

        import os
        output = './mesh/'
        if not os.path.exists(output):
            os.makedirs(output)
        fname = os.path.join(output, 'quad_mesh_2.vtu')
        mesh.to_vtk(fname=fname)

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

        KE, KTr, lambda_, mu = self. materialInfo()
        print("KE:", KE.shape, "\n", KE.round(4))

if __name__ == "__main__":

    ts = TopLevelset()
    ts.optimize()

