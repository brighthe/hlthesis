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

if __name__ == "__main__":

    ts = TopLevelset()
    ts.optimize()

