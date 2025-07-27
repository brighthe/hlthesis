import gmsh
import numpy as np
from fealpy.mesh import TriangleMesh
import matplotlib.pyplot as plt
# Initialize gmsh
gmsh.initialize()

# Create a new model
model = gmsh.model

# Set the model's dimension to 2D
model.add("2D")

lc = 0.1

# Define the parameters of the rectangle
L1 = 1  # Length of the first rectangle
H1 = 0.5  # Height of the first rectangle
L2 = 0.5  # Length of the second rectangle
H2 = 1  # Height of the second rectangle
gap = 0.1  # Gap between the rectangles

x = np.array([0.0, 4.0])
y = np.array([0.11, 0.13, 0.24, 0.26, 0.49, 0.51, 0.74, 0.76, 0.87, 0.89])

box = model.occ.addRectangle(0.0, 0.0, 0.0, 4.0, 1.0)
lines = []
line = []
for i in range(len(y)):
    p0 = model.occ.addPoint(0.0, y[i], 0.0)
    p1 = model.occ.addPoint(4.0, y[i], 0.0)
    l = model.occ.addLine(p0, p1)
    lines.append(l)
    line.append(l)
y = np.array([37/200, 0.375, 0.625, (76+87)/200])
for i in range(len(y)):
    p0 = model.occ.addPoint(0.0, y[i], 0.0)
    p1 = model.occ.addPoint(4.0, y[i], 0.0)
    l = model.occ.addLine(p0, p1)
    lines.append(l)
model.occ.fragment([(2, box)], [(1, l) for l in lines])

# Get the common edge of the two rectangles
model.occ.synchronize()

model.mesh.setSize(gmsh.model.getEntities(0), 0.05)
#
model.mesh.field.add("Distance", 1)
model.mesh.field.setNumbers(1, "CurvesList", line)
model.mesh.field.setNumber(1, "Sampling", 1000)

model.mesh.field.add("Threshold", 2)
model.mesh.field.setNumber(2, "InField", 1)
model.mesh.field.setNumber(2, "DistMin", 0.01)
model.mesh.field.setNumber(2, "DistMax", 0.06)
model.mesh.field.setNumber(2, "SizeMin", 0.015)
model.mesh.field.setNumber(2, "SizeMax", 0.05)
model.mesh.field.setNumber(2, "StopAtDistMax", 0)
model.mesh.field.setAsBackgroundMesh(2)

#model.mesh.filed.add("Constant", 3)
#model.mesh.field.setNumber(3, "VIn", 1/16)
#model.mesh.field.setNumber(3, "VOut", 1/8)

model.mesh.field.setAsBackgroundMesh(2)
model.mesh.generate(2)

# Visualize the mesh using Gmsh's built-in viewer
#gmsh.fltk.initialize()
#gmsh.fltk.run()
node = gmsh.model.mesh.get_nodes()[1].reshape(-1, 3)[:, :2]
for yy in y:
    flag = np.abs(node[:, 1] - yy)<1e-14
    node[flag, 1] = node[flag, 1] + (np.random.rand(flag.sum())-0.5)*1.5e-2
NN = node.shape[0]

nid2tag = gmsh.model.mesh.get_nodes()[0]
tag2nid = np.zeros(NN*2, dtype = np.int_)
tag2nid[nid2tag] = np.arange(NN)

cell = gmsh.model.mesh.get_elements(2, -1)[2][0].reshape(-1, 3)
cell = tag2nid[cell]
mesh = TriangleMesh(node, cell)

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes, aspect=1, linewidths=0.15, cellcolor=[128/256, 230/256, 115/256])
plt.savefig("bbb.png", dpi=400)
plt.show()

gmsh.finalize()

