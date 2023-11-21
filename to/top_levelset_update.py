import os

from fealpy.mesh import QuadrangleMesh


nelx = 60
nely = 30
mesh = QuadrangleMesh.from_box(box = [0, nelx, 0, nely], nx = nelx, ny = nely)
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)
fname = os.path.join(output, 'quad_mesh.vtu')
mesh.to_vtk(fname=fname)

