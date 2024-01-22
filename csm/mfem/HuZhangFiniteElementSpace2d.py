import numpy as np

from fealpy.functionspace.Function import Function
from fealpy.quadrature import FEMeshIntegralAlg

from fealpy.functionspace.femdof import multi_index_matrix2d

from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

from fealpy.decorator import cartesian, barycentric
from fealpy.functionspace import LagrangeFiniteElementSpace

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

class HuZhangFiniteElementSpace2d():
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
            self.lspace = LagrangeFiniteElementSpace(mesh, p)
        else:
            self.lspace = space
        self.cellmeasure = mesh.entity_measure('cell')
        self.integralalg = FEMeshIntegralAlg(mesh, p+3, cellmeasure=self.cellmeasure)
        self.integrator = self.integralalg.integrator

    @barycentric
    def basis(self, bc):
        """
        @brief 基函数
        """
        p    = self.p
        mesh = self.mesh
        node = mesh.entity("node")
        cell = mesh.entity("cell")
        NC   = mesh.number_of_cells()
        midx = mesh.multi_index_matrix(p, 2)

        ldof  = self.dof.number_of_local_dofs()
        gdof  = self.dof.number_of_global_dofs()
        eldof = self.dof.number_of_local_dofs("edge")
        cldof = self.dof.number_of_local_dofs("cell")

        c2frame = self.cell_frame() #(NC, 7, 3, 3)
        sphi = self.lspace.basis(bc) #(NQ, NC, ldof//3)
        phi = np.zeros(bc.shape[:-1]+(NC, ldof, 3), dtype=np.float_) #(NQ, NC, ldof, 3)

        ## 节点 
        sndofidx = [0, ldof//3-p-1, ldof//3-1] # 三个顶点自由度在标量元中的编号
        for i in range(3):
            frame_ni   = c2frame[:, i]                        # (NC, 3, 3) 
            sphi_ni    = sphi[..., sndofidx[i], None, None]   # (NQ, NC, 1, 1) 
            phi[..., 3*i:3*i+3, :] = sphi_ni*frame_ni         # (NQ, NC, 3, 3) 

        ## 边
        if eldof>0:
            eldof_2 = eldof//2
            Ndof = eldof*3+9
            for i in range(3):
                N = i*eldof+9
                edofidx  = np.where(midx[:, i]==0)[0][1:-1] 
                sphi_ei  = sphi[..., edofidx, None]    # (NQ, NC, eldof_2, 1)     
                frame_ei = c2frame[:, i+3, None]        # (NC, 1, 3, 3)
                phi[..., N:N+eldof_2, :] = sphi_ei*frame_ei[..., 0, :] #(NQ, NC, eldof_2, 3) 
                phi[..., N+eldof_2:N+eldof, :] = sphi_ei*frame_ei[..., 1, :] #(NQ, NC, eldof_2, 3) 
                phi[..., Ndof+i*eldof_2:Ndof+(i+1)*eldof_2, :] = sphi_ei*frame_ei[..., 2, :] #(NQ, NC, eldof_2, 3) 

        cldof -= 3*(eldof//2)
        if cldof>0:
            cldof_3 = cldof//3
            cdofidx = np.where(np.all(midx>0, axis=1))[0]
            frame_c = c2frame[:, 6, None]              # (NC, 1, 3, 3)
            sphi_c  = sphi[..., cdofidx, None]        # (NQ, NC, cldof//3, 1)     
            phi[..., -cldof:-2*cldof_3, :]    = sphi_c*frame_c[..., 0, :] # (NQ, NC, cldof//3, 3)
            phi[..., -2*cldof_3:-cldof_3, :] = sphi_c*frame_c[..., 1, :] # (NQ, NC, cldof//3, 3)
            phi[..., -cldof_3:, :]           = sphi_c*frame_c[..., 2, :] # (NQ, NC, cldof//3, 3)
        return phi

    def cell_frame(self): #TODO 更高效的实现: 只需要边界点标架和边上的标架 
        """
        @breif 每个单元有 7 个标架，三个顶点，三条边，内部 
        """
        c    = np.sqrt(2)/2.0
        mesh = self.mesh

        NC = mesh.number_of_cells()
        NN = mesh.number_of_nodes()

        cell = mesh.entity('cell')
        edge = mesh.entity('edge')

        c2e = mesh.ds.cell_to_edge()
        e2n = mesh.edge_unit_normal()
        e2t = mesh.edge_unit_tangent()

        c2en = e2n[c2e] #(NC, 3, 2)
        c2et = e2t[c2e] #(NC, 3, 2)

        c2frame = np.zeros((NC, 7, 3, 3), dtype=np.float_)
        ## 顶点和单元内部的标架
        c2frame[:] = np.identity(3)

        ## 边上标架
        c2frame[:, 3:6, 0, 0] =  c2en[..., 0]*c2en[..., 0]
        c2frame[:, 3:6, 0, 1] =  c2en[..., 0]*c2en[..., 1]
        c2frame[:, 3:6, 0, 2] =  c2en[..., 1]*c2en[..., 1]

        c2frame[:, 3:6, 1, 0] =  2*c*c2en[..., 0]*c2et[..., 0]
        c2frame[:, 3:6, 1, 1] =  c*(c2en[..., 0]*c2et[..., 1]+c2en[..., 1]*c2et[..., 0])
        c2frame[:, 3:6, 1, 2] =  2*c*c2en[..., 1]*c2et[..., 1]

        c2frame[:, 3:6, 2, 0] =  c2et[..., 0]*c2et[..., 0]
        c2frame[:, 3:6, 2, 1] =  c2et[..., 0]*c2et[..., 1]
        c2frame[:, 3:6, 2, 2] =  c2et[..., 1]*c2et[..., 1]

        ## 边界点的标架
        is_corner_node = self.is_corner_node
        is_bd_edge = self.mesh.ds.boundary_edge_flag()
        is_bd_node = self.mesh.ds.boundary_node_flag()

        NBN = np.sum(is_bd_node)
        bdnidxmap = np.zeros(NN, dtype=np.int_)
        bdnidxmap[is_bd_node] = np.arange(NBN)
        n2t = np.zeros((NBN, 2), dtype=np.float_) # 边界点的切向
        n2n = np.zeros((NBN, 2), dtype=np.float_) # 边界点的法向

        n2n[bdnidxmap[edge[is_bd_edge, 0]]]  = e2n[is_bd_edge]
        n2n[bdnidxmap[edge[is_bd_edge, 1]]] += e2n[is_bd_edge]
        n2n *= 0.5
        n2t[:, 0] = -n2n[:, 1]
        n2t[:, 1] =  n2n[:, 0]

        for i in range(3):
            flag = is_bd_node[cell[:, i]] & (~is_corner_node[cell[:, i]])
            c2nn = n2n[bdnidxmap[cell[flag, i]]]
            c2nt = n2n[bdnidxmap[cell[flag, i]]]
            c2frame[flag, i, 0, 0] =  c2nn[:, 0]*c2nn[:, 0]
            c2frame[flag, i, 0, 1] =  c2nn[:, 0]*c2nn[:, 1]
            c2frame[flag, i, 0, 2] =  c2nn[:, 1]*c2nn[:, 1]

            c2frame[flag, i, 1, 0] =  2*c*c2nn[:, 0]*c2nt[:, 0]
            c2frame[flag, i, 1, 1] =  c*(c2nn[:, 0]*c2nt[:, 1]+c2nn[:, 1]*c2nt[:, 0])
            c2frame[flag, i, 1, 2] =  2*c*c2nn[:, 1]*c2nt[:, 1]

            c2frame[flag, i, 2, 0] =  c2nt[:, 0]*c2nt[:, 0]
            c2frame[flag, i, 2, 1] =  c2nt[:, 0]*c2nt[:, 1]
            c2frame[flag, i, 2, 2] =  c2nt[:, 1]*c2nt[:, 1]

        return c2frame

    @barycentric
    def div_basis(self, bc):
        """
        @brief 基函数
        """
        p    = self.p
        mesh = self.mesh
        node = mesh.entity("node")
        cell = mesh.entity("cell")
        NC   = mesh.number_of_cells()
        midx = mesh.multi_index_matrix(p, 2)

        ldof  = self.dof.number_of_local_dofs()
        gdof  = self.dof.number_of_global_dofs()
        eldof = self.dof.number_of_local_dofs("edge")
        cldof = self.dof.number_of_local_dofs("cell")

        c2frame = self.cell_frame() #(NC, 7, 3, 3)
        sgphi = self.space.grad_basis(bc) #(NQ, NC, ldof//3, 2)
        dphi = np.zeros(sphi.shape[:-1]+(ldof, 2), dtype=np.float_) #(NQ, NC, ldof, 2)

        ## 节点 
        sndofidx = [0, ldof//3-p-1, ldof//3-1] # 三个顶点自由度在标量元中的编号
        for i in range(3):
            frame_ni   = c2frame[:, i]                          # (NC, 3, 3) 
            sgphi_ni    = sgphi[..., sndofidx[i], None]         # (NQ, NC, 1, 2) 
            gphi[..., 3*i:3*i+3, 0] = np.sum(phi_ni*frame_ni[..., :2], axis=-1) # (NQ, NC, 3) 
            gphi[..., 3*i:3*i+3, 1] = np.sum(phi_ni*frame_ni[..., 1:], axis=-1) # (NQ, NC, 3) 

        ## 边
        eldof_2 = eldof//2
        Ndof = eldof*3+9
        for i in range(3):
            N = i*eldof+9
            edofidx  = np.where(midx[:, i]==0)[0][1:-1] 
            sgphi_ei  = sgphi[..., edofidx, None]    # (NQ, NC, eldof_2, 1, 2)     
            frame_ei = c2frame[: i+3, None]        # (NC, 1, 3, 3)
            gphi[..., N:N+eldof_2, 0] = np.sum(sgphi_ei*frame_ei[..., 0, :2], axis=-1) #(NQ, NC, eldof_2) 
            gphi[..., N:N+eldof_2, 1] = np.sum(sgphi_ei*frame_ei[..., 0, 1:], axis=-1) #(NQ, NC, eldof_2) 

            gphi[..., N+eldof_2:N+eldof, 0] = np.sum(sgphi_ei*frame_ei[..., 1, :2], axis=-1) #(NQ, NC, eldof_2) 
            gphi[..., N+eldof_2:N+eldof, 1] = np.sum(sgphi_ei*frame_ei[..., 1, 1:], axis=-1) #(NQ, NC, eldof_2) 

            gphi[..., Ndof+i*eldof_2:Ndof+(i+1)*eldof_2, 0] = np.sum(sgphi_ei*frame_ei[..., 2, :2], axis=-1) #(NQ, NC, eldof_2) 
            gphi[..., Ndof+i*eldof_2:Ndof+(i+1)*eldof_2, 1] = np.sum(sgphi_ei*frame_ei[..., 2, 1:], axis=-1) #(NQ, NC, eldof_2) 

        if cldof>0:
            cldof_3 = cldof//3
            cdofidx = np.where(np.all(midx>0, axis=1))[0]
            frame_c = c2frame[: 6, None]              # (NC, 1, 3, 3)
            sgphi_c  = sgphi[..., cdofidx, None]        # (NQ, NC, cldof//3, 1)     
            gphi[..., -cldof:2*cldof_3, 0]    = np.sum(sgphi_c*frame_c[..., 0, :2], axis=-1)# (NQ, NC, cldof//3)
            gphi[..., -cldof:2*cldof_3, 1]    = np.sum(sgphi_c*frame_c[..., 0, 1:], axis=-1)# (NQ, NC, cldof//3)
            gphi[..., -2*cldof_3:-cldof_3, 0] = np.sum(sgphi_c*frame_c[..., 1, :2], axis=-1) # (NQ, NC, cldof//3)
            gphi[..., -2*cldof_3:-cldof_3, 1] = np.sum(sgphi_c*frame_c[..., 1, 1:], axis=-1) # (NQ, NC, cldof//3)
            gphi[..., -cldof_3:, 0]           = np.sum(sgphi_c*frame_c[..., 2, :2], axis=-1) # (NQ, NC, cldof//3)
            gphi[..., -cldof_3:, 1]           = np.sum(sgphi_c*frame_c[..., 2, 1:], axis=-1) # (NQ, NC, cldof//3)
        return val

    def cell_to_dof(self):
        return self.dof.cell2dof

    def number_of_global_dofs(self):
        return self.dof.number_of_global_dofs()

    def number_of_local_dofs(self, doftype='all'):
        return self.dof.number_of_local_dofs(doftype)

    @barycentric
    def value(self, uh, bc, index=np.s_[:]):
        '''@
        @brief 计算一个有限元函数在每个单元的 bc 处的值
        @param bc : (..., GD+1)
        @return val : (..., NC, GD)
        '''
        phi = self.basis(bc)
        c2d = self.dof.cell_to_dof()
        # uh[c2d].shape = (NC, ldof); phi.shape = (..., NC, ldof, GD)
        val = np.einsum("cl, ...clk->...ck", uh[c2d], phi, optimize=True)
        return val

    @barycentric
    def div_value(self, uh, bc, index=np.s_[:]):
        pass

    @barycentric
    def grad_value(self, uh, bc, index=np.s_[:]):
        pass

    @barycentric
    def edge_value(self, uh, bc, index=np.s_[:]):
        pass

    @barycentric
    def face_value(self, uh, bc, index=np.s_[:]):
        pass

    def mass_matrix(self):
        mesh = self.mesh
        NC = mesh.number_of_cells()
        ldof = self.dof.number_of_local_dofs()
        gdof = self.dof.number_of_global_dofs()
        cm = self.cellmeasure
        c2d = self.dof.cell_to_dof() #(NC, ldof)

        bcs, ws = self.integrator.get_quadrature_points_and_weights()
        phi = self.basis(bcs) #(NQ, NC, ldof, 3)
        num = np.array([1, 2, 1])
        mass = np.einsum("qclg, qcdg, c, q, g->cld", phi, phi, cm, ws, num)

        I = np.broadcast_to(c2d[:, :, None], shape=mass.shape)
        J = np.broadcast_to(c2d[:, None, :], shape=mass.shape)
        M = csr_matrix((mass.flat, (I.flat, J.flat)), shape=(gdof, gdof))
        return M 

    def tr_mass_matrix(self):
        mesh = self.mesh
        NC = mesh.number_of_cells()
        ldof = self.dof.number_of_local_dofs()
        gdof = self.dof.number_of_global_dofs()
        cm = self.cellmeasure
        c2d = self.dof.cell_to_dof() #(NC, ldof)

        bcs, ws = self.integrator.get_quadrature_points_and_weights()
        phi = self.basis(bcs) #(NQ, NC, ldof, 3)
        trphi = phi[..., 0]+phi[..., -1]
        mass = np.einsum("qcl, qcd, c, q->cld", trphi, trphi, cm, ws)

        I = np.broadcast_to(c2d[:, :, None], shape=mass.shape)
        J = np.broadcast_to(c2d[:, None, :], shape=mass.shape)
        M = csr_matrix((mass.flat, (I.flat, J.flat)), shape=(gdof, gdof))
        return M 

    def div_matrix(self, spaces):
        """
        @brief 计算 (div \sigma, u) u 属于 space, space 是一个元祖: sdof
          形式的空间
        """
        space = spaces[0]
        mesh = self.mesh
        NC = mesh.number_of_cells()
        ldof0 = self.dof.number_of_local_dofs()
        ldof1 = space.dof.number_of_local_dofs()
        gdof0 = self.dof.number_of_global_dofs()
        gdof1 = space.dof.number_of_global_dofs()
        cm = self.cellmeasure

        c2d = self.dof.cell_to_dof() #(NC, ldof)
        c2d_space = space.dof.cell_to_dof()
        c2d_space = np.c_[c2d_space, c2d_space+gdof1]

        bcs, ws = self.integrator.get_quadrature_points_and_weights()

        if space.basis.coordtype == 'barycentric':
            uphi = space.basis(bcs) #(NQ, NC, ldof1, 2)
        else:
            points = self.mesh.bc_to_point(bcs)
            uphi = space.basis(points)

        dphi = self.div_basis(bcs) #(NQ, NC, ldof, 2)
        A = np.zeros((NC, ldof0, ldof1*2), dtype=np.float_) 
        A[..., :ldof1] = np.einsum("qcl, qcm, c, q->clm", dphi[..., 0], uphi, cm, ws, optimize=True)
        A[..., ldof1:] = np.einsum("qcl, qcm, c, q->clm", dphi[..., 1], uphi, cm, ws, optimize=True)

        I = np.broadcast_to(c2d[:, :, None], shape=A.shape)
        J = np.broadcast_to(c2d_space[:, None, :], shape=A.shape)
        B = csr_matrix((A.flat, (I.flat, J.flat)), shape=(gdof0, gdof1*2))
        return B

    def source_vector(self, f):
        mesh = self.mesh
        cm = self.cellmeasure
        ldof = self.dof.number_of_local_dofs()
        gdof = self.dof.number_of_global_dofs()
        bcs, ws = self.integrator.get_quadrature_points_and_weights()
        c2d = self.dof.cell_to_dof() #(NC, ldof)

        p = mesh.bc_to_point(bcs) #(NQ, NC, GD)
        fval = f(p) #(NQ, NC, GD)

        phi = self.basis(bcs) #(NQ, NC, ldof, GD)
        val = np.einsum("qcg, qclg, q, c->cl", fval, phi, ws, cm)# (NC, ldof)
        vec = np.zeros(gdof, dtype=np.float_)
        np.add.at(vec, c2d, val)
        return vec

    def projection(self, f, method="L2"):
        M = self.mass_matrix()
        b = self.source_vector(f)
        x = spsolve(M, b)
        return self.function(array=x)

    def function(self, dim=None, array=None, dtype=np.float_):
        if array is None:
            gdof = self.dof.number_of_global_dofs()
            array = np.zeros(gdof, dtype=np.float_)
        return Function(self, dim=dim, array=array, coordtype='barycentric', dtype=dtype)

    def interplation(self, f):
        pass

    def L2_error(self, u, uh):
        '''@
        @brief 计算 ||u - uh||_{L_2}
        '''
        mesh = self.mesh
        cm = self.cellmeasure
        bcs, ws = self.integrator.get_quadrature_points_and_weights()
        p = mesh.bc_to_point(bcs) #(NQ, NC, GD)
        uval = u(p) #(NQ, NC, GD)
        uhval = uh(bcs) #(NQ, NC, GD)
        errval = np.sum((uval-uhval)*(uval-uhval), axis=-1)#(NQ, NC)
        val = np.einsum("qc, q, c->", errval, ws, cm)
        return np.sqrt(val)

    def set_neumann_bc(self, g):
        bcs, ws = self.integralalg.faceintegrator.get_quadrature_points_and_weights()

        edof = self.dof.number_of_local_dofs('edge')
        eidx = self.mesh.ds.boundary_edge_index()
        phi = self.edge_basis(bcs, index=eidx) #(NQ, NE0, edof, GD)
        e2n = self.mesh.edge_unit_normal(index=eidx)
        phi = np.einsum("qelg, eg->qel", phi, e2n) #(NQ, NE0, edof)

        point = self.mesh.edge_bc_to_point(bcs, index=eidx)
        gval = g(point) #(NQ, NE0)

        em = self.mesh.entity_measure("edge")[eidx]
        integ = np.einsum("qel, qe, e, q->el", phi, gval, em, ws)

        e2d = np.ones((len(eidx), edof), dtype=np.int_)
        e2d[:, 0] = edof*eidx
        e2d = np.cumsum(e2d, axis=-1)

        gdof = self.dof.number_of_global_dofs()
        val = np.zeros(gdof, dtype=np.float_)
        np.add.at(val, e2d, integ)
        return val

    def set_dirichlete_bc(self, uh, gD):
        """
        @brief gD : [[f0, eidx0], [f1, eidx1], [f2, eidx2], ...]
        """
        is_bd_dof = np.zeros(gdof, dtype=np.bool_)
        n2d = self.node_to_dof()
        e2d = self.edge_to_internal_dof()

        is_bd_node = self.mesh.ds.boundary_node_flag()
        is_bd_edge = self.mesh.ds.boundary_edge_flag()
        
        is_bd_dof[n2d[is_bd_node, 0]] = True
        is_bd_dof[e2d[is_bd_edge, :2]] = True

        is_bd_dof[n2d[self.is_corner_node]] = True

        uh[is_bd_dof] = 0
        return is_bd_dof



