import argparse
import time
import os

from linear_elasticity_model import BoxDomainData3d

from scipy.sparse.linalg import lgmres

from fealpy.functionspace import LagrangeFESpace as Space

from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import VectorSourceIntegrator
from fealpy.fem import BilinearForm
from fealpy.fem import LinearForm
from fealpy.fem import DirichletBC
from fealpy.fem import VectorNeumannBCIntegrator


## 参数解析
parser = argparse.ArgumentParser(description=
        """
        单纯形网格（三角形、四面体）网格上任意次有限元方法
        """)

parser.add_argument('--degree',
        default=1, type=int,
        help='Lagrange 有限元空间的次数, 默认为 1 次.')

parser.add_argument('--GD',
        default=3, type=int,
        help='模型问题的维数, 默认求解 3 维问题.')

parser.add_argument('--scale',
        default=1, type=float,
        help='网格变形系数，默认为 1')

parser.add_argument('--doforder',
        default='vdims', type=str,
        help='自由度排序的约定，默认为 vdims')

parser.add_argument('--output',
        default='./results/', type=str,
        help='output directory for the results. default is ./results/')

args = parser.parse_args()

# Check if the directory exists, if not, create it
if not os.path.exists(args.output):
    os.makedirs(args.output)

p = args.degree
GD = args.GD
scale = args.scale
doforder = args.doforder
output = args.output

pde = BoxDomainData3d()

domain = pde.domain()
mesh = pde.init_mesh()
NN = mesh.number_of_nodes()

start_time_1 = time.time()

# 构建双线性型，表示问题的微分形式
space = Space(mesh, p=p, doforder=doforder)
uh = space.function(dim=GD)
vspace = GD*(space, ) # 把标量空间张成向量空间
bform = BilinearForm(vspace)
bform.add_domain_integrator(LinearElasticityOperatorIntegrator(pde.lam, pde.mu))
bform.assembly()
A = bform.get_matrix()

end_time_1 = time.time()
elapsed_time_1 = end_time_1 - start_time_1
print("Elapsed time_1: {} seconds".format(elapsed_time_1))

start_time_2 = time.time()

# 构建单线性型，表示问题的源项
lform = LinearForm(vspace)
lform.add_domain_integrator(VectorSourceIntegrator(pde.source, q=1))
if hasattr(pde, 'neumann'):
    bi = VectorNeumannBCIntegrator(pde.neumann, threshold=pde.is_neumann_boundary, q=1)
    lform.add_boundary_integrator(bi)
lform.assembly()
F = lform.get_vector()

end_time_2 = time.time()
elapsed_time_2 = end_time_2 - start_time_2
print("Elapsed time_2: {} seconds".format(elapsed_time_2))

start_time_3 = time.time()

if hasattr(pde, 'dirichlet'):
    bc = DirichletBC(vspace, pde.dirichlet, threshold=pde.is_dirichlet_boundary)
    A, F = bc.apply(A, F, uh)

class IterationCounter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.niter = 0

    def __call__(self, rk=None):
        self.niter += 1
        if self._disp:
            print('iter %3i' % (self.niter))

def solve_system(A, F, tol=1e-8):
    counter = IterationCounter(disp = False)
    result, info = lgmres(A, F, tol = tol, callback = counter)
    return result

uh.flat[:] = solve_system(A, F)

end_time_3 = time.time()
elapsed_time_3 = end_time_3 - start_time_3
print("Elapsed time_3: {} seconds".format(elapsed_time_3))

if output != 'None':
    mesh.nodedata['uh'] = uh
    fname = os.path.join(output, f'test.vtu')
    mesh.to_vtk(fname=fname)
