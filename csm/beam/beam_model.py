import numpy as np

from fealpy.mesh import EdgeMesh

class CantileveredSimplySupportedBeam():
    def __init__(self):
        self.I = 1.186e-4 # Moment of Inertia - m^4
        self.A = 6650e-6 # Cross-sectional area - m^2
        self.E = 200e9 # Elastic Modulus newton - m^2

    #def init_mesh(self):
    #    mesh = EdgeMesh.generate_cantilevered_mesh()
    #    return mesh

    @classmethod
    def generate_cantilevered_mesh(cls):
        # Unit m
        node = np.array([
            [0], [5], [7.5]], dtype=np.float64)
        cell = np.array([
            [0, 1], [1, 2]], dtype=np.int_)
        mesh = cls(node, cell)

        mesh.meshdata['disp_bc'] = (np.array([0, 1], dtype = np.int_), np.zeros(2))
        mesh.meshdata['force_bc'] = (np.array([0, 1, 2], dtype = np.int_), 
                                     np.array([[-62500, -52083], [-93750, 39062], [-31250, 13021]], dtype = np.int_))

        return mesh 

