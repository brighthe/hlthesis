import numpy as np
import gmsh

from fealpy.mesh import TetrahedronMesh
from fealpy.decorator import cartesian

class BoxDomainData3d():
    def __init__(self):
        self.L = 1
        self.W = 0.2

        self.mu = 1
        self.rho = 1

        delta = self.W/self.L
        gamma = 0.4*delta**2
        beta = 1.25

        self.lam = beta
        self.g = gamma
        self.d = np.array([0.0, 0.0, -1.0])

    def domain(self):
        return [0.0, self.L, 0.0, self.W, 0.0, self.W]

    def to_TetrahedronMesh(self):
        ntags, vxyz, _ = gmsh.model.mesh.getNodes()
        node = vxyz.reshape((-1,3))
        vmap = dict({j:i for i,j in enumerate(ntags)})
        tets_tags,evtags = gmsh.model.mesh.getElementsByType(4)
        evid = np.array([vmap[j] for j in evtags])
        cell = evid.reshape((tets_tags.shape[-1],-1))
        return TetrahedronMesh(node,cell)

    def init_mesh(self):
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.013)
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.011)


        lc = 1.0
        gmsh.model.geo.addPoint(0,0,0,lc,1)
        gmsh.model.geo.addPoint(1,0,0,lc,2)
        gmsh.model.geo.addPoint(1,0,1,lc,3)
        gmsh.model.geo.addPoint(0,0,1,lc,4)
        gmsh.model.geo.addPoint(0,0,0.51,lc,5)
        gmsh.model.geo.addPoint(0.25,0,0.5,lc,6)
        gmsh.model.geo.addPoint(0,0,0.49,lc,7)

        gmsh.model.geo.addLine(1, 2, 1)
        gmsh.model.geo.addLine(2, 3, 2)
        gmsh.model.geo.addLine(3, 4, 3)
        gmsh.model.geo.addLine(4, 5, 4)
        gmsh.model.geo.addLine(5, 6, 5)
        gmsh.model.geo.addLine(6, 7, 6)
        gmsh.model.geo.addLine(7, 1, 7)

        gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, 6, 7], 1)
        gmsh.model.geo.addPlaneSurface([1], 1)
        gmsh.model.geo.extrude([(2,1)],0,1,0);

        gmsh.model.geo.synchronize()

        gmsh.model.mesh.generate(3)
        mesh = self.to_TetrahedronMesh()
        mesh.to_vtk(fname='model.vtu')

        gmsh.finalize()
        return mesh

    @cartesian
    def displacement(self, p):
        pass

    @cartesian
    def jacobian(self, p):
        pass

    @cartesian
    def strain(self, p):
        pass

    @cartesian
    def stress(self, p):
        pass

    @cartesian
    def source(self, p):
#        shape = len(p.shape[:-1])*(1,) + (-1, )
        val = self.d*self.g*self.rho
        return val 

    @cartesian
    def dirichlet(self, p):
        val = np.array([0.0, 0.0, 0.0])
        return val

    @cartesian
    def is_dirichlet_boundary(self, p):
        return np.abs(p[..., 0]) < 1e-12
    
    @cartesian
    def neumann(self, p, n):
        val = np.array([0.0, -50, 0.0], dtype=np.float64)
        return val

    @cartesian
    def is_neumann_boundary(self, p):
        x = p[..., 0]
        y = p[..., 1]
        z = p[..., 2]
        flag = np.abs(z - 1) < 1e-13
        return flag

