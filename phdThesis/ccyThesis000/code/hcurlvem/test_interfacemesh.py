import numpy as np
import matplotlib.pyplot as plt

from fealpy.geometry import CircleCurve, FoldCurve
from fealpy.mesh import MeshFactory as mf
from fealpy.mesh import PolygonMesh, HalfEdgeMesh2d 
from fealpy.geometry import CircleCurve, FoldCurve, DoubleCircleCurve, DoubleBandY 
from fealpy.mesh.interface_mesh_generator import interfacemesh2d 
from interface_mesh_2d import HalfEdgeMesh2dWithInterface, InterfaceMesher

interface = CircleCurve([0.5, 0.5], 0.31)
interface = DoubleCircleCurve(0.25, 0.35, np.array([1.5, 0.5]))
#interface = DoubleBandY()
#mesh = HalfEdgeMesh2dWithInterface([-4, 4, -3.1, 3.1], interface, 10, 10)
#mesh.uniform_refine()
mesh = interfacemesh2d([0, 4, 0, 1], interface, 32, 8, meshtype='tri')

#mesher = InterfaceMesher([0, 4, 0, 1], interface, 16, 4, 5)
#mesher = InterfaceMesher([-4, 4, -3, 3], interface, 4, 4, 7)
#mesh = mesher.get_mesh(0)

fig = plt.figure()
axes = fig.gca() 
mesh.add_plot(axes, aspect=1)
#mesh.find_cell(axes, showindex=True)
plt.show()
