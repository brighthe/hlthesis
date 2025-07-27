import gmsh

import matplotlib.pyplot as plt

from fealpy.mesh import TriangleMesh
from fealpy.mesh.gmsh_interface import gmsh_to_fealpy

def generate_mesh():
    gmsh.initialize()
    gmsh.model.add("mesh_generation")

    # 使用occ模块创建盒子和圆形
    gmsh.model.occ.addRectangle(0, 0, 0, 4, 1, 1)
    gmsh.model.occ.addDisk(1.25, 0.5, 0, 0.35, 0.35, 2)
    gmsh.model.occ.addDisk(1.75, 0.5, 0, 0.35, 0.35, 3)

    gmsh.model.occ.fuse([(2, 2)], [(2, 3)])
    gmsh.model.occ.fragment([(2, 1)], [(2, 2)])
    gmsh.model.occ.synchronize()

    # 设置区域的网格尺寸
    tags = gmsh.model.getEntities(0)  # 获取所有三维实体的标签
    gmsh.model.mesh.setSize(tags, 0.05)

    gmsh.option.setNumber("Mesh.Optimize", 0)


    # 生成网格
    gmsh.model.mesh.generate(2)
    mesh = gmsh_to_fealpy(gmsh.model, TriangleMesh, 2, 3)


    fig = plt.figure()
    axes = fig.gca()
    #axes.set_aspect(1)
    mesh.add_plot(axes, cellcolor='g')
    plt.savefig('aaa.svg')

    plt.show()

    gmsh.finalize()

generate_mesh()
