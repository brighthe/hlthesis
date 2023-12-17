import numpy as np
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

    def optimize(self):
        # Initialize optimization parameterse
        nelx, nely, volfrac = self._nelx, self._nely, self._volfrac
        mesh = self._mesh
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
        fname = os.path.join(output, 'quad_mesh_1.vtu')
        mesh.to_vtk(fname=fname)

if __name__ == "__main__":

    ts = TopSimp()
    ts.optimize()

    #import matplotlib.pyplot as plt
    #fig = plt.figure()
    #axes = fig.gca()
    #mesh.add_plot(axes)
    #mesh.find_node(axes, showindex=True, color='r', marker='o', markersize=2, fontsize=8, fontcolor='r')
    #plt.show()
