import os
import numpy as np

from fealpy.mesh import QuadrangleMesh
from fealpy.mesh import TriangleMesh

def cross_mesh(box=[0,2,0,1], nx=202, ny=101):
    qmesh = QuadrangleMesh.from_box(box=box, nx=nx, ny=ny)
    node = qmesh.entity('node')
    cell = qmesh.entity('cell')
    NN = qmesh.number_of_nodes()
    NC = qmesh.number_of_cells()
    bc = qmesh.entity_barycenter('cell') 
    newNode = np.r_['0', node, bc]

    newCell = np.zeros((4*NC, 3), dtype=np.int_) 
    newCell[0:NC, 0] = range(NN, NN+NC)
    newCell[0:NC, 1:3] = cell[:, 0:2]
    
    newCell[NC:2*NC, 0] = range(NN, NN+NC)
    newCell[NC:2*NC, 1:3] = cell[:, 1:3]

    newCell[2*NC:3*NC, 0] = range(NN, NN+NC)
    newCell[2*NC:3*NC, 1:3] = cell[:, 2:4]

    newCell[3*NC:4*NC, 0] = range(NN, NN+NC)
    newCell[3*NC:4*NC, 1:3] = cell[:, [3, 0]] 

    return TriangleMesh(newNode, newCell)

# Check if the directory exists, if not, create it
output = './mesh/'
if not os.path.exists(output):
    os.makedirs(output)

mesh = cross_mesh(box = [0,2,0,1], nx=202, ny=101)
fname = os.path.join(output, 'cross_mesh.vtu')
mesh.to_vtk(fname=fname)

from fealpy.functionspace import LagrangeFESpace

p = 1
vec = LagrangeFESpace(mesh, p=p, spacetype='C', doforder='vdims')
GD = mesh.geo_dimension()

Vvec = GD*(vec, )

def init_cantilever():
    tol = 1E-14  # tolerance for coordinate comparisons
    Lag, Nx, Ny = 40.0, 202, 101  # Specific parameters for cantilever
    lx, ly = 2.0, 1.0  # Dimensions of the domain

    # Mesh and function spaces
    mesh = RectangleMesh(Point(0.0, 0.0), Point(lx, ly), Nx, Ny, 'crossed')
    Vvec = VectorFunctionSpace(mesh, 'CG', 1)

    # Define boundary conditions
    class DirBd(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], .0)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)
    DirBd().mark(boundaries, 1)
    ds = Measure("ds")(subdomain_data=boundaries)
    bcd = [DirichletBC(Vvec, (0.0, 0.0), boundaries, 1)]

    # Load points
    Load = [Point(lx, 0.5)]

    # Initialize level set function
    XX, YY = np.meshgrid(np.linspace(0.0, lx, Nx + 1), np.linspace(0.0, ly, Ny + 1))
    phi_mat = -np.cos(8.0 / lx * pi * XX) * np.cos(4.0 * pi * YY) - 0.4 \
              + np.maximum(200.0 * (0.01 - XX ** 2 - (YY - ly / 2) ** 2), .0) \
              + np.maximum(100.0 * (XX + YY - lx - ly + 0.1), .0) + np.maximum(100.0 * (XX - YY - lx + 0.1), .0)

    return Lag, Nx, Ny, lx, ly, Load, 'cantilever', ds, bcd, mesh, phi_mat, Vvec
