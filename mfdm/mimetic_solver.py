import numpy as np

from scipy.sparse import diags, lil_matrix,csr_matrix
from scipy.sparse import spdiags

class Mimetic():
    def __init__(self, mesh):
        self.mesh = mesh
    
    def M_f(self, ac=None):
        """

        Results
        - result (ndarray, (NE, NE) )
        """
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

    def u_M_f(self, velocity, ac=None):
        """

        Results
        - result (ndarray, (NE, NE) )
        """
        mesh = self.mesh
        node = mesh.entity('node')
        edge = mesh.entity('edge')
        cell = mesh.entity('cell')
        print("cell:", cell)
        print("edge:", edge)
        NC = mesh.number_of_cells()
        NE = mesh.number_of_edges()
        cell2edge = mesh.ds.cell_to_edge()
        print("cell2edge:", cell2edge)
        cell2node = mesh.ds.cell_to_node()
        print("cell2node:", cell2node)
        norm = mesh.edge_unit_normal()
        edge_centers = mesh.entity_barycenter(etype=1)
        cell_centers = mesh.entity_barycenter(etype=2) # (NC, GD)
        flag = np.where(mesh.ds.cell_to_edge_sign().toarray(), 1, -1)
        print("flag:", flag.shape, "\n", flag)
        edge_measure = mesh.entity_measure(etype=1)
        cell_measure = mesh.entity_measure(etype=2)
        edge_norm = mesh.edge_unit_normal()
        print("edge_norm:", edge_norm.shape, "\n", edge_norm)
        result = np.zeros((NE, NE))
        if ac is None:
            ac = cell_measure

        u_nodes = velocity(node)
        print("u_nodes:", u_nodes.shape, "\n", u_nodes)
        u_edges = (u_nodes[edge[:, 0]] + u_nodes[edge[:, 1]]) / 2
        print("u_edges:", u_edges.shape, "\n", u_edges)
        #u0 = np.einsum('ij, ij -> i', u_edges, edge_norm)
        #print("u0:", u0.shape, "\n", u0)

        result = []
        for i in range(NC):
            #node_c = node[cell2node[i]] # (LNN, )
            #print("node_c:", node_c.shape, "\n", node_c)
            #u_nodes = velocity(node_c) # (LNN, 2)
            #print("u_nodes:", u_nodes)
            edge_c = edge[cell[i]] # (LNE, 2)
            #print("edge_c:", edge_c.shape, "\n", edge_c)
            #print(edge_c[:, 0])
            u_edges_c = (u_nodes[edge_c[:, 0]] + u_nodes[edge_c[:, 1]]) / 2 # (LNE, 2)
            #print("u_edges_c:", u_edges_c)
            tmp1 = edge_norm[cell2edge[i]]
            print("tmp1:", tmp1)
            tmp2 = flag[i, cell2edge[i]].reshape(-1, 1)
            print("tmp2:", tmp2)
            edge_norm_c = tmp1 * tmp2
            print("edge_norm_c:", edge_norm_c) # (LNE, 2)
            LNE = len(cell2edge[i])
            cell_out_flag = flag[i][cell2edge[i]] # (LNE, )
            print("cell_out_flag:", cell_out_flag)
            cell_edge_measure = edge_measure[cell2edge[i]] # (LNE, )
            cell_edge_centers = edge_centers[cell2edge[i]] # (LNE, GD)

            R = np.einsum('l, lg, l -> lg', cell_out_flag, cell_edge_centers-cell_centers[i, :],
                          cell_edge_measure) # (LNE, GD)
            N = norm[cell2edge[i], :] # (LNE, GD)

            M_consistency = R @ np.linalg.inv(R.T @ N) @ R.T
            M_stability = np.trace(R@R.T) /cell_measure[i] * (np.eye(LNE) - N @ np.linalg.inv(N.T @ N) @ N.T)
            M = M_consistency + M_stability # (LNE, LNE)

            #uhc = uh[cell2edge[i]] # (LNE, )
            uh_c = np.einsum('ij, ij -> i', u_edges_c, edge_norm_c) # (LNE, )
            #print("uh_c:", uh_c) # (LNE, )

            uhc_M = M @ uh_c # (LNE, )

            result.append(uhc_M)

        #print("result:", result)
        UM = np.zeros((NC, NE))
        for i, (cell_edges, uhc_M) in enumerate(zip(cell2edge, result)):
            UM[i, cell_edges] = uhc_M

        return UM

    
    def M_c(self):
        """

        Results
        - result (ndarray, (NC, NC) )
        """
        mesh = self.mesh
        cell_measure = mesh.entity_measure("cell")
        result = np.diag(cell_measure)

        return result
    
    def div_operator(self):
        """
        离散散度算子

        Results
        - result (ndarray, (NC, NE) )
        """
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
            result[i, cell2edge[i]] = cell_out_flag * cell_edge_measure/cell_measure[i]

        return result

    def gard_operator(self):
        """
        离散梯度算子

        Results
        - result (ndarray, (NE, NN) )
        """
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
        """

        Results
        - result (ndarray, (NC+NE, ) )
        """
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

    def source_neumann(self, fun):
        """

        Results
        - result (ndarray, (NC+NE, ) )
        """
        mesh = self.mesh
        NE = mesh.number_of_edges()
        cell_centers = mesh.entity_barycenter(etype=2)
        cell_measure = mesh.entity_measure('cell') 

        f = fun(cell_centers)
        b = cell_measure * f

        gamma = np.zeros(NE)
        result = np.hstack((-gamma, -b))

        return result

