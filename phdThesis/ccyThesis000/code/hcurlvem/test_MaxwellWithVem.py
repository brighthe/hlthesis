#!/usr/bin/env python3
# 

import time
import sys
import copy
import numpy as np
import matplotlib.pyplot as plt

from fealpy.mesh import TriangleMesh, HalfEdgeMesh2d, PolygonMesh
from fealpy.decorator import cartesian, barycentric
from fealpy.tools.show import showmultirate, show_error_table
from fealpy.boundarycondition import DirichletBC #导入边界条件包
from fealpy.pde.MaxwellPDE2d import SinData as PDE

from scipy.sparse.linalg import spsolve, cg
from HCurlVirtualElementSpace2d import HCurlVirtualElementSpace2d

def delete_boundary_edge(mesh):
    NN = mesh.number_of_nodes()
    NE = mesh.number_of_edges()
    NC = mesh.number_of_cells()
    cell, loc = mesh.entity('cell')

    cstart = mesh.ds.cellstart
    halfedge = mesh.entity('halfedge')
    node = mesh.entity('node')
    ver, cel, nex, per, opp = halfedge[:, 0], halfedge[:, 1], halfedge[:, 2], halfedge[:, 3], halfedge[:, 4]
    idx, = np.where(cel[opp] < cstart)

    v0 = node[ver[idx]]-node[ver[per[idx]]]
    v1 = node[ver[nex[idx]]]-node[ver[idx]]

    hidx = idx[np.abs(np.cross(v0, v1))<0.000001]

    nidx = ver[hidx]
    isMarkedNode = np.zeros(NN, dtype=np.bool_)
    isMarkedNode[nidx] = True
    node = node[~isMarkedNode]

    nidxmap = np.zeros(NN, dtype=np.int_)-1
    nidxmap[~isMarkedNode] = np.arange(np.sum(~isMarkedNode))

    mark = isMarkedNode[cell]
    mark = np.split(mark, loc[1:-1])

    f = lambda x : np.sum(x)
    num = np.array((list(map(f, mark))))
    num = np.cumsum(num)
    loc[1:] = loc[1:]-num
    cell = cell[~isMarkedNode[cell]]
    cell = nidxmap[cell]

    mesh = PolygonMesh(node, cell, loc)
    mesh = HalfEdgeMesh2d.from_mesh(mesh)
    return mesh

def get_mesh(i):
    mesh = HalfEdgeMesh2d.from_mesh(pde.init_mesh(nx = 2**i, ny = 2**i,
        meshtype='noconvex'))
    return mesh

def compute(mesh, pde):
    space = HCurlVirtualElementSpace2d(mesh)
    gdof = space.dof.number_of_global_dofs()
    NDof[i] = gdof

    bc = DirichletBC(space, pde.dirichlet) 

    M = space.mass_matrix()
    C = space.curl_matrix()
    b = space.source_vector(pde.source)
    B = C-M 

    B, b = bc.apply(B, b)
    uh = spsolve(B, b)
    # 计算误差
    err0 = space.L2_error(pde.solution, uh)
    err1 = space.curl_error(pde.curl_solution, uh)
    return err0, err1

pde = PDE()
maxit = 6
errorType = ['$|| u - u_h||_0$',
             '$||rot u - rot u_h||_0$']
errorMatrix = np.zeros((len(errorType), maxit), dtype=np.float_)
NDof = np.zeros(maxit, dtype=np.int_)

for i in range(maxit):
    mesh = get_mesh(i)
    errorMatrix[0, i], errorMatrix[1, i] = compute(mesh, pde)

showmultirate(plt, 2, NDof, errorMatrix,  errorType, propsize=20)
show_error_table(NDof, errorType, errorMatrix)
plt.show()

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.savefig("aaa.png", dpi=400)

