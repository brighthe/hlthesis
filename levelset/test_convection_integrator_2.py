import numpy as np

from fealpy.mesh import TriangleMesh
from fealpy.fem import ScalarConvectionIntegrator 
from fealpy.decorator import cartesian
from fealpy.functionspace import LagrangeFESpace 
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.fem import BilinearForm
from fealpy.fem import VectorConvectionIntegrator 
from fealpy.fem import LinearForm
from fealpy.fem import VectorSourceIntegrator
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
udegree = 1
qq = 1
mesh = TriangleMesh.from_unit_square(nx=ns, ny=ns)
uspace = LagrangeFESpace(mesh,p=udegree,doforder='sdofs')

Ouspace = LagrangeFiniteElementSpace(mesh,p=udegree)

u0=uspace.interpolate(velocity_field,dim=2)
Ou0=Ouspace.interpolation(velocity_field,dim=2)
so = uspace.interpolate(solution,dim=2)
Oso = Ouspace.interpolation(solution, dim=2)

Sbform = BilinearForm(uspace)
Sbform.add_domain_integrator(ScalarConvectionIntegrator(u0, qq))
Sbform.assembly()
C = Sbform.get_matrix()
lform = LinearForm((uspace,)*2)
lform.add_domain_integrator(VectorSourceIntegrator(solution, qq))
B = lform.assembly()
gdof = uspace.number_of_global_dofs()
print(np.sum(np.abs(C@Ou0[:,0] - B[:gdof])))

OC = Ouspace.convection_matrix(Ou0,qq).T
OB = Ouspace.source_vector(solution,dim=2,q=qq)
print(np.sum(np.abs(OC@Ou0[:,0] - OB[:,0] )))



print(np.sum(np.abs(C-OC)))
print(np.sum(np.abs(B-OB.flatten('F'))))

    
