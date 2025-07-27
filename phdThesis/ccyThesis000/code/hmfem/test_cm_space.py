import sys
import numpy as np
import matplotlib.pyplot as plt

from fealpy.mesh import TriangleMesh
from cm_conforming_fespace import CmConformingFiniteElementSpace2d

def test(m=0, p=1, kk=1):
    #mesh = TriangleMesh.from_one_triangle(meshtype='equ')
    mesh = TriangleMesh.from_box([0, 1, 0, 1], 1, 1)
    mesh.node[1, 1] = 1.2
    mesh.node[3, 0] = 0.8
    space = CmConformingFiniteElementSpace2d(mesh, m=m, p=p)
    
    uh = space.function()
    #uh[:] = np.random.rand(len(uh))
    uh[:] = 0
    uh[kk] = 10

    cell = mesh.entity('cell')
    edge = mesh.entity('edge')
    print(edge)
    print(cell)
    print(space.dof.cell_to_dof())

    bcs = np.array([[0, 1/2, 1/2]], dtype=np.float_)

    val = space.grad_m_value(uh, bcs, 1)
    print(val)


    fig = plt.figure()
    axes = fig.gca()
    mesh.add_plot(axes)
    mesh.find_node(axes, showindex=True)
    #plt.show()

if __name__ == "__main__" :
    m = int(sys.argv[1])
    p = int(sys.argv[2])
    kk = int(sys.argv[3])
    test(m, p, kk)

