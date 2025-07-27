import time
import string
import itertools
import numpy as np

from typing import Union, Callable

from fealpy.functionspace import BernsteinFESpace
from scipy.special import factorial, comb
from fealpy.common.tensor import *
from fealpy.functionspace.Function import Function 
from fealpy.decorator import barycentric
from fealpy.quadrature.TriangleQuadrature1 import TriangleQuadrature

from scipy.sparse import csc_matrix, csr_matrix
from scipy.linalg import solve_triangular

class CmConformingDof():
    def __init__(self, mesh, p, m):
        self.mesh = mesh
        self.m = m
        self.p = p

    def number_of_local_dofs(self, etype):
        p = self.p
        m = self.m
        if etype=='cell':
            return (p+1)*(p+2)//2
        if etype=='edge':
            ndof = (2*m+1)*(2*m+2)//2
            eidof = (m+1)*(2*p-7*m-2)//2
            return ndof + eidof
        if etype=='node':
            return (2*m+1)*(2*m+2)//2

    def number_of_internal_dofs(self, etype):
        p = self.p
        m = self.m
        if etype=='cell':
            ndof = (2*m+1)*(2*m+2)//2
            eidof = (m+1)*(2*p-7*m-2)//2
            return (p+1)*(p+2)//2 - (ndof+eidof)*3
        if etype=="edge":
            return (m+1)*(2*p-7*m-2)//2
        if etype=='node':
            return (2*m+1)*(2*m+2)//2

    def number_of_global_dofs(self):
        mesh = self.mesh
        NN = mesh.number_of_nodes()
        NE = mesh.number_of_edges()
        NC = mesh.number_of_cells()
        ndof = self.number_of_internal_dofs('node')
        eidof = self.number_of_internal_dofs('edge')
        cidof = self.number_of_internal_dofs('cell')
        return NN*ndof + eidof*NE + cidof*NC

    def node_to_dof(self):
        m = self.m
        mesh = self.mesh
        ndof = self.number_of_internal_dofs('node') 

        NN = mesh.number_of_nodes()
        n2d = np.arange(NN*ndof, dtype=np.float_).reshape(NN, ndof)
        return n2d.astype(np.int_)

    def edge_to_internal_dof(self):
        p = self.p
        m = self.m
        ndof = self.number_of_internal_dofs('node') 
        eidof = self.number_of_internal_dofs('edge') 

        mesh = self.mesh
        NN = mesh.number_of_nodes()
        NE = mesh.number_of_edges()
        e2id = np.arange(NN*ndof, NN*ndof + NE*eidof).reshape(NE, eidof)
        return e2id.astype(np.int_)

    def edge_to_dof(self):
        n2d = self.node_to_dof()
        e2id = self.edge_to_internal_dof()

        edge = self.mesh.entity('edge')
        e2d = np.c_[n2d[edge[:, 0]], n2d[edge[:, 1]], e2id]
        return e2d.astype(np.int_)

    def cell_to_dof(self):
        p, m = self.p, self.m
        mesh = self.mesh

        ldof = self.number_of_local_dofs('cell')
        cidof = self.number_of_internal_dofs('cell')
        ndof = self.number_of_internal_dofs('node') 
        eidof = self.number_of_internal_dofs('edge') 
        n2d = self.node_to_dof()
        e2id = self.edge_to_internal_dof()

        NN = mesh.number_of_nodes()
        NE = mesh.number_of_edges()
        NC = mesh.number_of_cells()
        cell = mesh.entity('cell')
        c2e = mesh.ds.cell_to_edge()
        c2eSign = mesh.ds.cell_to_edge_sign()
        c2dof = np.zeros((NC, ldof), dtype=np.int_)

        ## 内部自由度
        c2dof[:, ldof-cidof:] = np.arange(NN*ndof + NE*eidof, NN*ndof + NE*eidof
                + NC*cidof).reshape(NC, cidof)

        ## 顶点自由度
        for v in range(3):
            c2dof[:, ndof*v:ndof*(v+1)] = n2d[cell[:, v]]
        ## 边上的自由度
        for e in range(3):
            N = ndof*3+eidof*e
            flag = ~c2eSign[:, e]
            c2dof[:, N:N+eidof] = e2id[c2e[:, e]] 
            n0, n1 = 0, p-4*m-1
            for l in range(m+1):
                c2dof[flag, N+n0:N+n0+n1] = c2dof[flag, N+n0:N+n0+n1][:, ::-1]
                n0 += n1
                n1 += 1
        return c2dof

    def is_boundary_dof(self, r):
        p = self.p
        m = self.m
        nn = (r+1)*(r+2)//2
        ne = 0
        for i in range(r+1): ne += p-4*m-1+i

        gdof = self.number_of_global_dofs()
        isBdNode = self.mesh.ds.boundary_node_flag()
        isBdEdge = self.mesh.ds.boundary_edge_flag()
        bndof = self.node_to_dof()[isBdNode][:, :nn]
        bedof = self.edge_to_internal_dof()[isBdEdge][:, :ne]
        isBdDof = np.zeros(gdof, dtype=np.bool_)
        isBdDof[bndof] = True
        isBdDof[bedof] = True
        return isBdDof

class CmConformingFiniteElementSpace2d(object):
    def __init__(self, mesh, p=1, m=0):
        assert(p>4*m)
        self.mesh = mesh
        self.m = m
        self.p = p
        self.bspace = BernsteinFESpace(mesh, p)
        self.dof = CmConformingDof(mesh, p, m)
        self.coeff = self.coefficient_matrix()
        self.integrator = TriangleQuadrature(p+4)

    def coefficient_matrix(self):
        p = self.p
        m = self.m
        mesh = self.mesh

        NC = mesh.number_of_cells()
        ldof = self.dof.number_of_local_dofs('cell')

        tem = np.eye(ldof, dtype=np.float_)
        coeff = np.tile(tem, (NC, 1, 1)) 

        ## 重新对多重指标进行编号
        multiIndex = mesh.multi_index_matrix(p, 2)
        midx2num = lambda a : (a[:, 1]+a[:, 2])*(1+a[:, 1]+a[:, 2])//2 + a[:, 2]

        S02m0 = multiIndex[multiIndex[:, 0]>=p-2*m] # S_0^{2m}(delta_0)
        S02m1 = S02m0[:, [1, 0, 2]]
        S02m2 = S02m0[:, [1, 2, 0]]
        S1m0 = multiIndex[(multiIndex[:, 0]<=m)&np.all(multiIndex[:, 1:]<p-2*m, axis=1)]
        S1m0 = S1m0[:, [0, 2, 1]][::-1]
        S1m1 = S1m0[:, [2, 0, 1]]
        S1m2 = S1m0[:, [1, 2, 0]]
        S2 = multiIndex[np.all(multiIndex>m, axis=1)]

        dof2midx = np.r_[S02m0, S02m1, S02m2, S1m0, S1m1, S1m2, S2] #自由度对应的多重指标
        dof2num = midx2num(dof2midx) # 自由度对应的多重指标的编号
        dof2num = np.argsort(dof2num)

        S02m = [S02m0, S02m1, S02m2]
        ## 节点处的基函数
        for v in range(3):
            flag = np.ones(3, dtype=np.bool_)
            flag[v] = False
            for gamma in S02m[v]:
                i = midx2num(gamma[None, :]) # 泛函编号
                i = dof2num[i]
                gamma = gamma[flag]
                r = np.sum(gamma)
                betas = mesh.multi_index_matrix(r, 2)
                for beta in betas:
                    if np.any(gamma-beta[flag]<0):
                        continue
                    sign = (-1)**beta[v]
                    beta[v] += p-r
                    j = midx2num(beta[None, :]) # bernstein 多项式编号
                    j = dof2num[j]
                    beta = beta[flag]
                    coeff[:, i, j] = sign*factorial(r)*comb(p, r)*np.prod(comb(gamma, beta)) 

        glambda = mesh.grad_lambda() #(NC, 3, 2)
        c2e = mesh.ds.cell_to_edge() #(NC, 3)
        n = mesh.edge_normal()[c2e] #(NC, 3, 2)
        N = np.einsum('cfd, ced->cef', glambda, n, optimize=True) #(NC, 3, 3)

        S1m = [S1m0, S1m1, S1m2]
        ## 边上的基函数
        for e in range(3):
            flag = np.ones(3, dtype=np.bool_)
            flag[e] = False
            for gamma in S1m[e]:
                i = midx2num(gamma[None, :]) # 泛函编号
                i = dof2num[i]
                r = gamma[e]
                betas = mesh.multi_index_matrix(r, 2)
                for beta in betas:
                    val = factorial(r)**2*comb(p, r)/np.prod(factorial(beta))
                    val = val*np.prod(N[:, e]**beta, axis=1)
                    beta += gamma
                    beta[e] -= r
                    j = midx2num(beta[None, :]) # bernstein 多项式编号
                    j = dof2num[j]
                    coeff[:, i, j] = val[:, None] 
        coeff = np.linalg.inv(coeff)
        coeff = np.transpose(coeff, (0, 2, 1))

        ## 自由度转换系数
        node = mesh.entity('node')
        cell = mesh.entity('cell')

        t0 = node[cell[:, 2]]-node[cell[:, 1]]
        t1 = node[cell[:, 0]]-node[cell[:, 2]]
        t2 = node[cell[:, 1]]-node[cell[:, 0]]

        tdelta = np.zeros((NC, 3, 2, 2), dtype=np.float_)
        tdelta[:, 0, 0] = t2
        tdelta[:, 0, 1] = -t1
        tdelta[:, 1, 0] = -t2
        tdelta[:, 1, 1] = t0
        tdelta[:, 2, 0] = t1
        tdelta[:, 2, 1] = -t0

        ndof = self.dof.number_of_local_dofs('node') 
        coeff1 = np.zeros((NC, 3*ndof, 3*ndof), dtype=np.float_)
        for v in range(3):
            tdeltav = tdelta[:, v]
            k, kk = 1, 1 # k 是当前顶点的第 k 个自由度，kk 是当前顶点的第r行，第一个自由度的编号
            coeff1[:, ndof*v, ndof*v] = 1 
            for r in range(2*m+1)[1:]:
                ff = lambda x: np.array([int(ss) for ss in np.binary_repr(x, r)], dtype=np.int_)
                idx = np.array(list(map(ff, np.arange(2**r))))
                isSymteryValue = np.ones(len(idx), dtype=np.bool_)
                for i in range(r)[1:]:
                    isSymteryValue = isSymteryValue & (idx[:, i] <= idx[:, i-1])
                multiidx = mesh.multi_index_matrix(r, 1)
                num = factorial(r)/np.prod(factorial(multiidx), axis=1)
                for midx in multiidx:
                    tdeltav_sym = symmetry_span_array(tdeltav, midx).reshape(-1, 2**r)
                    coeff1[:, ndof*v+k, ndof*v+kk:ndof*v+kk+isSymteryValue.sum()] = num[None, :]*tdeltav_sym[:, isSymteryValue]
                    k += 1
                kk += r+1

        coeff[:, :3*ndof] = np.einsum("cji, cjk->cik", coeff1, coeff[:,
            :3*ndof], optimize=True)
        return coeff[:, :, dof2num];

    def basis(self, bcs):
        coeff = self.coeff
        bphi = self.bspace.basis(bcs)
        return np.einsum('cil, qcl->qci', coeff, bphi, optimize=True)

    def grad_m_basis(self, bcs, m):
        coeff = self.coeff
        bgmphi = self.bspace.grad_m_basis(bcs, m)
        return np.einsum('cil, qclg->qcig', coeff, bgmphi, optimize=True)

    def grad_m_matrix(self, m):
        p = self.p
        gdof = self.dof.number_of_global_dofs()
        cell2dof = self.dof.cell_to_dof()
        integrator = self.integrator
        bcs, ws = integrator.get_all_gauss_point_and_weight()

        GD = self.mesh.geo_dimension()
        idx = self.mesh.multi_index_matrix(m, GD-1)
        num = factorial(m)/np.prod(factorial(idx), axis=1)

        cm = self.mesh.entity_measure("cell")
        gmphi = self.grad_m_basis(bcs, m)
        M = np.einsum('qclg, qcmg, g, q, c->clm', gmphi, gmphi, num, ws, cm, optimize=True)

        I = np.broadcast_to(cell2dof[:, :, None], M.shape) 
        J = np.broadcast_to(cell2dof[:, None, :], M.shape) 
        M = csr_matrix((M.flat, (I.flat, J.flat)), shape=(gdof, gdof), dtype=np.float_)
        return M

    def mass_matrix(self):
        # TODO 优化效率, Bernstein 基的质量矩阵每个单元是相同的
        p = self.p
        gdof = self.dof.number_of_global_dofs()
        cell2dof = self.dof.cell_to_dof()
        integrator = self.integrator
        bcs, ws = integrator.get_all_gauss_point_and_weight()

        cm = self.mesh.entity_measure("cell")
        phi = self.basis(bcs)
        M = np.einsum('qcl, qcm, q, c->clm', phi, phi, ws, cm, optimize=True)

        I = np.broadcast_to(cell2dof[:, :, None], M.shape) 
        J = np.broadcast_to(cell2dof[:, None, :], M.shape) 
        M = csr_matrix((M.flat, (I.flat, J.flat)), shape=(gdof, gdof), dtype=np.float_)
        return M

    def source_vector(self, f):
        ## TODO 优化效率, 这里的 phi 冗余了
        p = self.p
        gdof = self.dof.number_of_global_dofs()
        cell2dof = self.dof.cell_to_dof()
        integrator = self.integrator
        bcs, ws = integrator.get_all_gauss_point_and_weight()

        cm = self.mesh.entity_measure("cell")
        phi = self.basis(bcs)

        point = self.mesh.bc_to_point(bcs)
        fval = f(point)
        inte = np.einsum('qcl, qc, q, c->cl', phi, fval, ws, cm, optimize=True)

        F = np.zeros(gdof, dtype=np.float_)
        np.add.at(F, cell2dof, inte)
        return F

    def function(self, dim=None, array=None, dtype=np.float64):
        return Function(self, dim=dim, array=array, 
                coordtype='barycentric', dtype=dtype)

    def array(self, dim=None, dtype=np.float64):
        gdof = self.dof.number_of_global_dofs()
        return np.zeros(gdof, dtype=dtype)

    @barycentric
    def value(self, 
            uh: np.ndarray, 
            bc: np.ndarray, 
            index: Union[np.ndarray, slice]=np.s_[:]
            ) -> np.ndarray:
        """
        @brief Computes the value of a finite element function `uh` at a set of
        barycentric coordinates `bc` for each mesh cell.

        @param uh: numpy.ndarray, the dof coefficients of the basis functions.
        @param bc: numpy.ndarray, the barycentric coordinates with shape (NQ, TD+1).
        @param index: Union[numpy.ndarray, slice], index of the entities (default: np.s_[:]).
        @return numpy.ndarray, the computed function values.

        This function takes the dof coefficients of the finite element function `uh` and a set of barycentric
        coordinates `bc` for each mesh cell. It computes the function values at these coordinates
        and returns the results as a numpy.ndarray.
        """
        phi = self.basis(bc) # (NQ, 1, ldof)
        cell2dof = self.dof.cell_to_dof()
        val = np.einsum('qcl, cl->qc', phi, uh[cell2dof], optimize=True)
        return val

    @barycentric
    def grad_m_value(self, uh, bcs, m):
        gmphi = self.grad_m_basis(bcs, m) # (NQ, 1, ldof)
        cell2dof = self.dof.cell_to_dof()
        val = np.einsum('qclg, cl->qcg', gmphi, uh[cell2dof])
        return val

    @barycentric
    def grad_value(self, uh, bcs):
        return self.grad_m_value(uh, bcs, 1)

    def boundary_interpolate(self, gD: list, uh: np.ndarray, threshold) -> np.ndarray:
        #TODO 只处理的边界为 0 的情况
        r = len(gD)
        isDDof = self.dof.is_boundary_dof(r-1)
        uh[isDDof] = 0
        return isDDof

    set_dirichlet_bc = boundary_interpolate

    def number_of_global_dofs(self):
        return self.dof.number_of_global_dofs()

