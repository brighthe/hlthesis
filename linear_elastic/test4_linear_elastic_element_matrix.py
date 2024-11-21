from fealpy.fem import LinearElasticityOperatorIntegrator
from fealpy.fem import BilinearForm
from fealpy.functionspace import LagrangeFESpace
from fealpy.mesh import TriangleMesh

NX = 4
NY = 4
mesh = TriangleMesh.from_box(box=[0, 1, 0, 1], nx=NX, ny=NY)
NC = mesh.number_of_cells()
print("NC:", NC)
NN = mesh.number_of_nodes()
print("NN:", NN)

p = 1
space = LagrangeFESpace(mesh, p=p, ctype='C', doforder='sdofs')


GD = 2
uh = space.function(dim=GD)
print("uh:", uh.shape)

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