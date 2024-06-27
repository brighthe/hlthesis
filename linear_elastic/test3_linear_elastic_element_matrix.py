from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import BilinearForm
from fealpy.functionspace import LagrangeFESpace as Space
from fealpy.mesh import UniformMesh2d

nelx, nely = 3, 2
domain = [0, 3, 0, 2]
hx = (domain[1] - domain[0]) / nelx
hy = (domain[3] - domain[2]) / nely
mesh = UniformMesh2d(extent=(0, nelx, 0, nely), h=(hx, hy), origin=(domain[0], domain[2]))
p = 1
space = Space(mesh, p=p, doforder='vdims')
GD = 2
uh = space.function(dim=GD)
# 材料参数
E0 = 1.0  # 弹性模量
nu = 0.3  # 泊松比
lambda_ = (E0 * nu) / ((1 + nu) * (1 - 2 * nu))
mu = E0 / (2 * (1 + nu))

integrator = LinearElasticityOperatorIntegrator(lam=lambda_, mu=mu, q=p+3)
vspace = GD * (space,)
bform = BilinearForm(vspace)
bform.add_domain_integrator(integrator)
KK = integrator.assembly_cell_matrix(space=vspace)
print("KK:", KK.shape, "\n", KK[0].round(4))