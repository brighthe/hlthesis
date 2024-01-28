import numpy as np

from scipy.sparse import diags, lil_matrix,csr_matrix
from scipy.sparse import spdiags

class Mimetic():
    def __init__(self, mesh):
        self.mesh = mesh
    
    def M_f(self, ac=None):
        mesh = self.mesh
        NC = mesh.number_of_cells()
        NE = mesh.number_of_edges()
        cell2edge = mesh.ds.cell_to_edge()
        norm = mesh.edge_unit_normal()
        edge_centers = mesh.entity_barycenter(etype=1)
        cell_centers = mesh.entity_barycenter(etype=2) # (NC, GD)
        flag = np.where(mesh.ds.cell_to_edge_sign().toarray(), 1, -1)
        edge_measure = mesh.entity_measure(etype=1)
        cell_measure = mesh.entity_measure(etype=2)
        result = np.zeros((NE, NE))
        if ac is None:
            ac = cell_measure
        
        for i in range(NC):
            LNE = len(cell2edge[i])
            cell_out_flag = flag[i][cell2edge[i]] # (LNE, )
            cell_edge_measure = edge_measure[cell2edge[i]] # (LNE, )
            cell_edge_centers = edge_centers[cell2edge[i]] # (LNE, GD)

            R = np.einsum('l, lg, l -> lg', cell_out_flag, cell_edge_centers-cell_centers[i, :],
                          cell_edge_measure) # (LNE, GD)
            N = norm[cell2edge[i], :] # (LNE, GD)

            M_consistency = R @ np.linalg.inv(R.T @ N) @ R.T
            M_stability = np.trace(R@R.T) /cell_measure[i] * (np.eye(LNE) - N @ np.linalg.inv(N.T @ N) @ N.T)
            M = M_consistency + M_stability # (LNE, LNE)
            indexi, indexj = np.meshgrid(cell2edge[i], cell2edge[i])
            result[indexi, indexj] += M 

        return result
    
    def M_c(self):
        mesh = self.mesh
        cell_measure = mesh.entity_measure("cell")
        result = np.diag(cell_measure)

        return result
    
    def div_operator(self):
        mesh = self.mesh
        NC = mesh.number_of_cells()
        NE = mesh.number_of_edges()
        cell_measure = mesh.entity_measure(etype=2)
        edge_measure = mesh.entity_measure(etype=1)
        cell2edge = mesh.ds.cell_to_edge()
        flag = np.where(mesh.ds.cell_to_edge_sign().toarray(), 1, -1)
        result = np.zeros((NC, NE))
        for i in range(NC):
            cell_out_flag = flag[i][cell2edge[i]] # (LNE, )
            cell_edge_measure = edge_measure[cell2edge[i]] # (LNE, )
            result[i, cell2edge[i]] = cell_out_flag * cell_edge_measure/cell_measure[i] # (NC, NE)

        return result

    def gard_operator(self):
        mesh = self.mesh
        edge = mesh.entity('edge')
        NE = mesh.number_of_edges()
        NN = mesh.number_of_nodes()
        edge_measure = mesh.entity_measure(etype=1)
        
        gradh = np.zeros((NE, NN))

        for i in range(NE):
            gradh[i, edge[i, 0]] = -1 / edge_measure[i]
            gradh[i, edge[i, 1]] = 1 / edge_measure[i]

        return gradh
    
    def source(self, fun, gddof, D):
        mesh = self.mesh
        NE = mesh.number_of_edges()
        cell_centers = mesh.entity_barycenter(etype=2)
        cell_measure = mesh.entity_measure('cell') 
        
        edge_measure = mesh.entity_measure('edge') 
        edge_centers = mesh.entity_barycenter(etype=1)

        f = fun(cell_centers)
        b = cell_measure * f

        gamma = np.zeros(NE)
        gamma[gddof] = edge_measure[gddof] * D(edge_centers[gddof])
        result = np.hstack((-gamma, -b))

        return result
