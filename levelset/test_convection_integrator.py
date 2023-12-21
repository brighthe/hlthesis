import numpy as np

from fealpy.mesh import TriangleMesh
# from tssim.part import *
# from tssim.part.Assemble import Assemble
from fealpy.fem import ScalarConvectionIntegrator 
from fealpy.decorator import cartesian
from fealpy.functionspace import LagrangeFESpace 
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.fem import BilinearForm
from fealpy.fem import VectorConvectionIntegrator 
from fealpy.fem import LinearForm
from fealpy.fem import VectorSourceIntegrator
from scipy.sparse import csr_matrix
'''
@cartesian
def velocity_field(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = np.sin(y)*np.sin(x) 
    u[..., 1] = -np.sin(x)*np.sin(y)
    return u
@cartesian
def solution(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = np.sin(x)*np.sin(y)**2*np.cos(x) - np.sin(x)**2*np.sin(y)*np.cos(y)
    u[..., 1] = -np.sin(x)*np.sin(y)**2*np.cos(x) + np.sin(x)**2*np.sin(y)*np.cos(y)
    return u
'''
@cartesian
def velocity_field(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = x**2 
    u[..., 1] = y**2
    return u
@cartesian
def solution(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = 2*x**3
    u[..., 1] = 2*y**3
    return u


ns = 32
udegree = 3
GD = 2

mesh = TriangleMesh.from_unit_square(nx=ns, ny=ns)
uspace = LagrangeFESpace(mesh, p=udegree, doforder='sdofs')
Ouspace = LagrangeFiniteElementSpace(mesh, p=udegree)

u0 = uspace.interpolate(velocity_field, dim=2)
Ou0 = Ouspace.interpolation(velocity_field, dim=2)
so = uspace.interpolate(solution, dim=2)
Oso = Ouspace.interpolation(solution, dim=2)

Sbform = BilinearForm(uspace)
Sbform.add_domain_integrator(ScalarConvectionIntegrator(c=u0, q=4))
Sbform.assembly()
C = Sbform.get_matrix()
lform = LinearForm(GD*(uspace,))
lform.add_domain_integrator(VectorSourceIntegrator(solution, q=4))
B = lform.assembly()
print("B:", B.shape)
gdof = uspace.number_of_global_dofs()
print(np.sum(np.abs(C@Ou0[:, 0] - B[:gdof])))

OC = Ouspace.convection_matrix(c=Ou0).T
OB = Ouspace.source_vector(solution, dim=2)
print("OB:", OB.shape)
print(np.sum(np.abs(OC@Ou0[:, 0] - OB[:, 0] )))

qf = mesh.integrator(4) 
bcs, ws = qf.get_quadrature_points_and_weights()
phi = uspace.basis(bcs)
gphi = uspace.grad_basis(bcs)
ps = mesh.bc_to_point(bcs)
cm = mesh.entity_measure('cell')
coef = velocity_field(ps)
CC = np.einsum('q, qck, qcn, qcmn, c -> ckm', ws, phi, coef, gphi, cm) 
cell2dof = uspace.cell_to_dof()
I = np.broadcast_to(cell2dof[:, :, None], shape=CC.shape)
J = np.broadcast_to(cell2dof[:, None, :], shape=CC.shape)
NC = csr_matrix((CC.flat, (I.flat, J.flat)), shape=(gdof, gdof))
val = solution(ps)
BB = np.einsum('q, qcn, qck, c -> cnk', ws, val, phi, cm)
NB = np.zeros((gdof, GD), dtype=np.float64)
for i in range(GD):
    np.add.at(NB[:, i], cell2dof, BB[:, i, :])
print("NB:", NB.shape)
print(np.sum(np.abs(NC@Ou0[:, 0] - NB[:, 0] )))


#C1 = assemble.matrix([udegree, 1], [udegree, 0], Ou0(assemble.bcs)[...,0])
#C2 = assemble.matrix([udegree, 2], [udegree, 0], Ou0(assemble.bcs)[...,1])
#AC = C1+C2
#bcs = assemble.bcs
#ipoint = mesh.bc_to_point(bcs)
#value = solution(ipoint)
#fb1 = assemble.vector([udegree,0],value[...,0])
#fb2 = assemble.vector([udegree,0],value[...,1])

#print(np.sum(np.abs(AC@Ou0[:,0] - fb1)))
#print(np.sum(np.abs(AC-C)))
#print(np.sum(np.abs(AC-OC)))
print("C-OC:", np.sum(np.abs(C-OC)))
print("C-NC:", np.sum(np.abs(C-NC)))
print("OC-NC:", np.sum(np.abs(OC-NC)))
print("B-OB:", np.sum(np.abs(B-OB.flatten('F'))))
print("B-NB:", np.sum(np.abs(B-NB.flatten('F'))))
print("OB-NB:", np.sum(np.abs(OB.flatten('F')-NB.flatten('F'))))
