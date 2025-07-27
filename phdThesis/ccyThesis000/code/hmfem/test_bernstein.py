import numpy as np
from fealpy.functionspace import BernsteinFESpace
from fealpy.mesh import TetrahedronMesh, TriangleMesh

#mesh = TriangleMesh.from_box([0, 1, 0, 1], 1, 1)
#mesh = TetrahedronMesh.from_box([0, 1, 0, 1, 0, 1], 5, 5, 5)
mesh = TriangleMesh.from_polygon_gmsh([[0, 0], [1, 0], [1, 1], [0, 1]], 0.01)
print(mesh.number_of_cells())
space = BernsteinFESpace(mesh, p=14)

integrator = mesh.integrator(5, 'cell')
bcs, ws = integrator.get_quadrature_points_and_weights()

uh = space.function()

import time
t0 = time.time()
val0 = space.grad_m_basis(bcs, 3)
t1 = time.time()
val1 = space.grad_m_basis_0(bcs, 3)
t2 = time.time()
print("diff : ", np.max(np.abs(val0[..., 0]-val1[..., 0])))
print(t2-t1)
print(t1-t0)
