import numpy as np
import matplotlib.pyplot as plt

from fealpy.mesh import TriangleMesh

mesh = TriangleMesh.from_box([0, 1, 0, 1], 1, 1)
print(mesh.entity('edge'))
print(mesh.entity('cell'))

ip = mesh.interpolation_points(3)
print(mesh.cell_to_ipoint(3))

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes, cellcolor='w')
mesh.find_node(axes, node=ip, showindex=True)
#mesh.find_edge(axes, showindex=True)
#mesh.find_cell(axes, showindex=True)
plt.savefig('out1.svg')
plt.show()





