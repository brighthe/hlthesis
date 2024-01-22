import numpy as np

from scipy.sparse import csr_matrix

from fealpy.functionspace.femdof import multi_index_matrix2d
from fealpy.functionspace import LagrangeFESpace

from fealpy.quadrature import FEMeshIntegralAlg

class HuZhangDof2d():
    def __init__(self, mesh, p):
        self.p = p
        self.mesh = mesh
        self.multiindex = multi_index_matrix2d(p)

    def node_to_dof(self, index=np.s_[:]):
        """
        @brief 节点上自由度
        """
        NN = self.mesh.number_of_nodes()
        nldof = self.number_of_local_dofs("node")
        n2d = np.arange(nldof*NN).reshape(NN, nldof)[index]
        return n2d 

    def edge_to_internal_dof(self, index=np.s_[:]):
        """
        @brief 边内部的自由度
        """
        NN = self.mesh.number_of_nodes()
        NE = self.mesh.number_of_edges()

        nldof = self.number_of_local_dofs("node")
        eldof = self.number_of_local_dofs("edge")

        e2d = np.arange(nldof*NN, nldof*NN+eldof*NE).reshape(NE, eldof)[index]
        return e2d

    def number_of_local_dofs(self, doftype='all'):
        p = self.p
        if doftype == 'all': # number of all dofs on a cell 
            return 3*(p+1)*(p+2)//2 
        elif doftype in {'cell', 2}: # number of dofs inside the cell 
            return 3*(p+1)*(p+2)//2 - 9 - 6*(p-1) 
        elif doftype in {'face', 'edge', 1}: # number of dofs on each edge 
            return 2*(p-1)
        elif doftype in {'node', 0}: # number of dofs on each node
            return 3

    def number_of_global_dofs(self):
        NC = self.mesh.number_of_cells()
        NE = self.mesh.number_of_edges()
        NN = self.mesh.number_of_nodes()
        nldof = self.number_of_local_dofs("node")
        eldof = self.number_of_local_dofs(doftype='edge')
        cldof = self.number_of_local_dofs(doftype='cell')
        return NE*eldof + NC*cldof + NN*nldof

    def cell_to_dof(self):
        p = self.p
        NC = self.mesh.number_of_cells()
        NE = self.mesh.number_of_edges()
        NN = self.mesh.number_of_nodes()

        cldof = self.number_of_local_dofs('cell')
        eldof = self.number_of_local_dofs('edge')
        nldof = self.number_of_local_dofs('node')
        ldof = self.number_of_local_dofs()
        gdof = self.number_of_global_dofs()

        n2dof = self.node_to_dof()
        e2dof = self.edge_to_internal_dof()

        cell = self.mesh.entity("cell")
        edge = self.mesh.entity("edge")
        c2e = self.mesh.ds.cell_to_edge()

        c2d = np.zeros((NC, ldof), dtype=np.int_)

        # 顶点自由度
        c2d[:, 0:3] = n2dof[cell[:, 0]]
        c2d[:, 3:6] = n2dof[cell[:, 1]]
        c2d[:, 6:9] = n2dof[cell[:, 2]]

        # 边自由度
        locEdge = np.array([[1, 2], [0, 2], [0, 1]])
        if eldof>0:
            c2d[:, 9:eldof+9]           = e2dof[c2e[:, 0]]
            c2d[:, 9+eldof:eldof*2+9]   = e2dof[c2e[:, 1]]
            c2d[:, 9+eldof*2:eldof*3+9] = e2dof[c2e[:, 2]]

            flag = cell[:, locEdge[0, 0]] == edge[c2e[:, 0], 0]
            c2d[flag, 9:eldof//2+9]                 = e2dof[c2e[flag, 0]][:, :eldof//2]
            c2d[flag, 9+eldof//2:eldof+9]           = e2dof[c2e[flag, 0]][:, eldof//2:]
            flag = cell[:, locEdge[1, 0]] == edge[c2e[:, 1], 0]
            c2d[flag, 9+eldof:eldof+eldof//2+9]     = e2dof[c2e[flag, 1]][:, :eldof//2]
            c2d[flag, eldof+eldof//2+9:eldof*2+9]   = e2dof[c2e[flag, 1]][:, eldof//2:]
            flag = cell[:, locEdge[2, 0]] == edge[c2e[:, 2], 0]
            c2d[flag, 9+eldof*2:eldof*2+eldof//2+9] = e2dof[c2e[flag, 2]][:, :eldof//2]
            c2d[flag, eldof*2+eldof//2+9:eldof*3+9] = e2dof[c2e[flag, 2]][:, eldof//2:]

        # 单元自由度
        if cldof>0:
            c2d[:, eldof*3+9:] = np.arange(NE*eldof+NN*nldof, gdof).reshape(NC, -1)
        return c2d

    @property
    def cell2dof(self):
        return self.cell_to_dof()

    def boundary_dof(self):
        eidx = self.mesh.ds.boundary_edge_index()
        e2d = self.edge_to_dof(index=eidx)
        return e2d.reshape(-1)

    def is_boundary_dof(self):
        bddof = self.boundary_dof()

        gdof = self.number_of_global_dofs()
        flag = np.zeros(gdof, dtype=np.bool_)

        flag[bddof] = True
        return flag

class HuZhangMassOperatorIntegrator:
    def __init__(self, mesh, p, is_corner_node=None, space=None):
        self.p = p
        self.mesh = mesh
        self.dof = HuZhangDof2d(mesh, p)

        ## 角点
        NN = mesh.number_of_nodes()
        if is_corner_node is None:
            self.is_corner_node = np.zeros(NN, dtype=np.bool_)
        else:
            self.is_corner_node = is_corner_node

        ## 标量空间
        if space is None:
            self.lspace = LagrangeFESpace(mesh, p)
        else:
            self.lspace = space
        self.cellmeasure = mesh.entity_measure('cell')
        self.integralalg = FEMeshIntegralAlg(mesh, p+3, cellmeasure=self.cellmeasure)
        self.integrator = self.integralalg.integrator

    def mass_matrix(self, space):
        mesh = space.mesh
        NC = mesh.number_of_cells()
        ldof = space.dof.number_of_local_dofs()
        gdof = space.dof.number_of_global_dofs()
        cm = space.cellmeasure
        c2d = space.dof.cell_to_dof() #(NC, ldof)
        bcs, ws = space.integrator.get_quadrature_points_and_weights()
        phi = space.basis(bcs) #(NQ, NC, ldof, 3)
        num = np.array([1, 2, 1])
        mass = np.einsum("qclg, qcdg, c, q, g->cld", phi, phi, cm, ws, num)

        I = np.broadcast_to(c2d[:, :, None], shape=mass.shape)
        J = np.broadcast_to(c2d[:, None, :], shape=mass.shape)
        M = csr_matrix((mass.flat, (I.flat, J.flat)), shape=(gdof, gdof))
        return mass, M

