import numpy as np
from scipy.sparse import csr_matrix
from fealpy.functionspace.Function import Function

from fealpy.quadrature import PolygonMeshIntegralAlg

class VEDof():
    """!
    @brief 最低阶 H(curl) 协调的虚单元空间
    """
    def __init__(self, mesh):
        self.mesh = mesh

    def number_of_local_dofs(self):
        """
        @brief 返回每个单元的自由度个数
        """
        return self.mesh.ds.number_of_vertices_of_cells()

    def number_of_global_dofs(self):
        return int(self.mesh.number_of_edges())

    def cell_to_dof(self):
        """
        @brief 返回每个单元的自由度编号 cell2dof, cell2dofLoc, 
               其中第 i 个单元的自由度为 cell2dof[cell2dof[i] : cell2dof[i+1]] 
        """
        return self.mesh.ds.cell_to_edge()

class HCurlVirtualElementSpace2d():
    """!
    @brief 最低阶 H(curl) 协调的虚单元空间
    """
    def __init__(self, mesh):
        self.mesh = mesh
        self.dof = VEDof(mesh)

        self.cellbarycenter = mesh.entity_barycenter('cell')
        self.cellmeasure = mesh.entity_measure('cell')
        self.integralalg = PolygonMeshIntegralAlg(self.mesh, 2, 
                    cellmeasure=self.cellmeasure,
                    cellbarycenter=self.cellbarycenter)

        def fx(p, index) : return p[..., 0]
        def fy(p, index) : return p[..., 1]
        self.cellbarycenter[:, 0] = self.integralalg.integral(fx,
                celltype=True)/self.cellmeasure
        self.cellbarycenter[:, 1] = self.integralalg.integral(fy,
                celltype=True)/self.cellmeasure

        self.PI = self.matrix_PI()

    def curl_val(self, cell2dofLoc):
        """
        @brief 返回每个单元上每个基函数的旋度值
        """
        mesh = self.mesh
        cellmea = self.cellmeasure
        val = (mesh.ds.cell_to_edge_sign().astype(np.float_)-0.5)*2
        val = np.split(val, cell2dofLoc[1:-1])

        f = lambda x : x[0]/x[1] 
        val = list(map(f, zip(val, cellmea)))
        return val

    def curl_matrix(self, beta=None, dtype=np.float_):
        """!
        @brief 返回每个单元的旋度矩阵 C : 
                        C_{Kij} = (beta curl phi_i, curl phi_j)_K
        """
        dof  = self.dof
        mesh = self.mesh
        ldof = dof.number_of_local_dofs()
        gdof = dof.number_of_global_dofs()
        cellbar = self.cellbarycenter 
        cellmea = self.cellmeasure
        betaK = 1 if beta is None else beta(cellbar).reshape(-1)

        cell2edgeSign = (mesh.ds.cell_to_edge_sign().astype(np.float_)-0.5)*2
        cell2edge, cell2edgeLoc = mesh.ds.cell_to_edge()
        c2d       = np.hsplit(cell2edge, cell2edgeLoc[1:-1])
        c2eS      = np.split(cell2edgeSign, cell2edgeLoc[1:-1])

        f1 = lambda x: x[1]*(x[0][:, None]@x[0][None, :]).flat
        f2 = lambda x: np.repeat(x, x.shape[0])
        f3 = lambda x: np.tile(x, x.shape[0])

        C = np.concatenate(list(map(f1, zip(c2eS, betaK/cellmea))))
        I = np.concatenate(list(map(f2, c2d)))
        J = np.concatenate(list(map(f3, c2d)))
        C = csr_matrix((C, (I, J)), shape=(gdof, gdof), dtype=dtype)
        return C

    def matrix_PI(self):
        """!
        @brief 
        """
        mesh = self.mesh
        NC   = mesh.number_of_cells()
        NE   = mesh.number_of_edges()
        cellmea = self.cellmeasure
        cellbar = self.cellbarycenter
        edgebar = mesh.entity_barycenter("edge")
        edgemea = mesh.entity_measure("edge")
        edge2cell = mesh.ds.edge_to_cell()

        # edgeint[:, 0] 每条边上左边单元的函数 (x1 - barx1) 的积分
        # edgeint[:, 1] 每条边上右边单元的函数 (x1 - barx1) 的积分
        # edgeint[:, 2] 每条边上左边单元的函数 (barx0 - x0) 的积分
        # edgeint[:, 3] 每条边上右边单元的函数 (barx0 - x0) 的积分
        edgeint = np.zeros((NE, 4), dtype=np.float_)
        edgeint[:, 0] = edgebar[:, 1]-cellbar[edge2cell[:, 0], 1]
        edgeint[:, 1] = cellbar[edge2cell[:, 1], 1]-edgebar[:, 1]
        edgeint[:, 2] = cellbar[edge2cell[:, 0], 0]-edgebar[:, 0]
        edgeint[:, 3] = edgebar[:, 0]-cellbar[edge2cell[:, 1], 0]

        cell2edge, cell2edgeLoc = mesh.ds.cell_to_edge()

        # c2eint[:, 0] 每个单元每条边函数 (x1 - barx1) 的积分
        # c2eint[:, 1] 每个单元每条边函数 (barx0 - x0) 的积分
        c2eint = np.zeros((cell2edgeLoc[-1], 2), dtype=np.float_)
        c2eint[cell2edgeLoc[edge2cell[:, 1]]+edge2cell[:, 3], 0] = edgeint[:, 1]
        c2eint[cell2edgeLoc[edge2cell[:, 0]]+edge2cell[:, 2], 0] = edgeint[:, 0]
        c2eint[cell2edgeLoc[edge2cell[:, 1]]+edge2cell[:, 3], 1] = edgeint[:, 3]
        c2eint[cell2edgeLoc[edge2cell[:, 0]]+edge2cell[:, 2], 1] = edgeint[:, 2]
        c2eint = np.split(c2eint, cell2edgeLoc[1:-1])

        f = lambda x : (x[0][:, None]@x[3][None, :]/x[1] - x[2].T)/x[1]
        f = lambda x : (-x[1].T)/x[0]
        PI = list(map(f, zip(cellmea, c2eint)))
        return PI

    def mass_matrix(self, alpha=None, dtype=np.float_):
        dof  = self.dof
        mesh = self.mesh
        gdof = dof.number_of_global_dofs()
        cellbar = self.cellbarycenter 
        cellmea = self.cellmeasure
        edgetan = mesh.edge_tangent()
        edgemea = mesh.entity_measure('edge')
        alphaK = 1 if alpha is None else alpha(cellbar).reshape(-1)

        cell2edge, cell2edgeLoc = mesh.ds.cell_to_edge()
        c2d = np.hsplit(cell2edge, cell2edgeLoc[1:-1])
        tN  = np.split(edgetan[cell2edge], cell2edgeLoc[1:-1]) # \tilde{N}^K 

        f0 = lambda x : np.eye(x[0].shape[0]) - x[0]@x[1] 
        N = list(map(f0, zip(tN, self.PI)))

        f1 = lambda x: x[0]*(x[1].T@x[1]).flat + np.sqrt(x[3])*((x[2].T)@x[2]/edgemea[x[4], None]).flat
        f2 = lambda x: np.repeat(x, x.shape[0])
        f3 = lambda x: np.tile(x, x.shape[0])

        M = np.concatenate(list(map(f1, zip(alphaK*cellmea, self.PI, N, cellmea,
            c2d))))
        I = np.concatenate(list(map(f2, c2d)))
        J = np.concatenate(list(map(f3, c2d)))
        M = csr_matrix((M, (I, J)), shape=(gdof, gdof), dtype=dtype)
        return M

    def projection_to_smspace(self, uh):
        cell2dof, cell2dofLoc = self.dof.cell_to_dof()
        uhK = np.split(uh[cell2dof], cell2dofLoc[1:-1])

        f = lambda x : x[0]@x[1]
        uhK = np.concatenate(list(map(f, zip(self.PI, uhK)))).reshape(-1, 2)
        return uhK

    def L2_error(self, u, uh, celltype=False):
        cell2dof, cell2dofLoc = self.dof.cell_to_dof()
        uhK = np.split(uh[cell2dof], cell2dofLoc[1:-1])

        f = lambda x : x[0]@x[1]
        uhK = np.concatenate(list(map(f, zip(self.PI, uhK)))).reshape(-1, 2)

        def e(pp, cidx) : return np.sum((u(pp)-uhK[cidx][None, :])**2, axis=-1)
        if celltype:
            return self.integralalg.integral(e, celltype=True)
        else:
            return np.sqrt(self.integralalg.integral(e))

    def curl_error(self, cu, uh, celltype=False):
        cellmea = self.cellmeasure
        cell2dof, cell2dofLoc = self.dof.cell_to_dof()
        uhK = np.split(uh[cell2dof], cell2dofLoc[1:-1])
        val = self.curl_val(cell2dofLoc)

        f = lambda x : x[0]@x[1]
        cuh = np.array(list(map(f, zip(uhK, val))))

        def e(pp, cidx) : return (cu(pp)-cuh[cidx][None, :])**2
        if celltype:
            return self.integralalg.integral(e, celltype=True)
        else:
            return np.sqrt(self.integralalg.integral(e))

    def source_vector(self, f, dtype=np.float_):
        """
        @brief 右端向量
        """
        NC = self.mesh.number_of_cells()
        gdof = self.dof.number_of_global_dofs()
        cell2dof, _ = self.dof.cell_to_dof()

        def f0(pp, cidx) : return f(pp)[..., 0]
        def f1(pp, cidx) : return f(pp)[..., 1]

        F_ = np.zeros((NC, 2), dtype=dtype)
        F_[:, 0] = self.integralalg.integral(f0, celltype=True)
        F_[:, 1] = self.integralalg.integral(f1, celltype=True)

        f2 = lambda x : (x[0][None, :]@x[1]).flat
        FK = np.concatenate(list(map(f2, zip(F_, self.PI))))
        F = np.zeros(gdof, dtype=dtype)
        np.add.at(F, cell2dof, FK)
        return F

    def set_dirichlet_bc(self, gD, uh, threshold=None, q=None):
        """
        """
        p = 1
        mesh = self.mesh

        qf = self.integralalg.edgeintegrator if q is None else mesh.integrator(q, 'edge')
        bcs, ws = qf.get_quadrature_points_and_weights()

        index = self.mesh.ds.boundary_edge_index()

        t = mesh.edge_tangent(index=index)
        ps = mesh.bc_to_point(bcs, index=index)
        val = gD(ps, t) # (NQ, NBE)
        if len(val.shape)==3:
            val = np.einsum('...ed, ed->...e', val, t)

        gdof = self.number_of_global_dofs()
        uh[index, 0] = np.einsum('q, qe->e', ws, val, optimize=True)
        isDDof = np.zeros(gdof, dtype=np.bool_)
        isDDof[index] = True
        return isDDof

    def interpolation(self, u):
        eb = self.mesh.entity_barycenter('edge')
        uval = u(eb)
        t = self.mesh.edge_tangent()
        uI = self.function()
        uI[:] = np.sum(uval*t, axis=1)
        return uI

    def number_of_global_dofs(self):
        return self.dof.number_of_global_dofs()

    def function(self, dim=None, array=None, dtype=np.float64):
        return Function(self, dim=dim, array=array, coordtype='cartesian',
                        dtype=dtype)

    def get_curl_value_on_cell(self, uh):
        cellmea = self.cellmeasure
        cell2dof, cell2dofLoc = self.dof.cell_to_dof()
        uhK = np.split(uh[cell2dof], cell2dofLoc[1:-1])
        val = self.curl_val(cell2dofLoc)

        f = lambda x : x[0]@x[1]
        cuh = np.array(list(map(f, zip(uhK, val, cellmea))))
        return cuh

    def array(self, dim=None, dtype=np.float64):
        gdof = self.number_of_global_dofs()
        if dim is None:
            shape = gdof
        elif type(dim) is int:
            shape = (gdof, dim)
        elif type(dim) is tuple:
            shape = (gdof,) + dim
        return np.zeros(shape, dtype=dtype)








