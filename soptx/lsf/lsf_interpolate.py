
import numpy as np

from shape_gradient import TopLsfShapeGrad

# Cantilever 的默认参数
domain_width = 3
domain_hight = 2
nelx = 30
nely = 20
ew = domain_width / nelx
eh = domain_hight / nely
volReq = 0.5
stepLength = 2;
numReinit = 3
ts = TopLsfShapeGrad(domain_width = domain_width, domain_hight = domain_hight, \
                     nelx = nelx, nely = nely, \
                     volReq=volReq, stepLength=stepLength, numReinit=3)

fe_domain = [0, domain_width, 0, domain_hight]
fe_mesh = ts.generate_mesh(domain=fe_domain, nelx=nelx, nely=nely)
fe_NC = fe_mesh.number_of_cells()
fe_node = fe_mesh.entity('node') # 左下角按列增加
fe_cell = fe_mesh.entity('cell') # 左下角逆时针
fe_center_node = fe_mesh.entity_barycenter('cell')
fe_center_x = fe_center_node[:, 0]
fe_center_y = fe_center_node[:, 1]
fe_x = fe_node[:, 0]
fe_y = fe_node[:, 1]
print("fe_x:", fe_x.shape, "\n", fe_x.round(4))
print("fe_y:", fe_y.shape, "\n", fe_y.round(4))

ls_domain = [-ew/2, domain_width+ew/2, -eh/2, domain_hight+eh/2]
ls_mesh = ts.generate_mesh(domain=ls_domain, nelx=nelx+1, nely=nely+1)
ls_NC = ls_mesh.number_of_cells()
ls_NN = ls_mesh.number_of_nodes()
ls_node = ls_mesh.entity('node') # 左下角按列增加
ls_cell = ls_mesh.entity('cell') # 左下角逆时针
ls_center_node = ls_mesh.entity_barycenter('cell')
ls_x = ls_node[:, 0]
ls_y = ls_node[:, 1]
print("ls_x:", ls_x.shape, "\n", ls_x.round(4))
print("ls_y:", ls_y.shape, "\n", ls_y.round(4))

from scipy.interpolate import griddata
def f(x, y):
    return np.sin(x) + np.cos(y)
ls_values = f(ls_x, ls_y)
# 插值到有限元网格上
points = np.vstack((ls_x, ls_y)).T
values = ls_values
xi = np.vstack((fe_x, fe_y)).T
fe_values_interpolated = griddata(points, values, xi, method='cubic')
print("fe_values_interpolated:", fe_values_interpolated)
# 有限元网格上的精确值
fe_values_exact = f(fe_x, fe_y)
print("fe_values_exact:", fe_values_exact)
# 计算误差
error = np.sum(np.abs(fe_values_interpolated - fe_values_exact))
print("error:", error)


ls_Phi = ts.init_lsf(mesh = ls_mesh)
print("ls_Phi0:", ls_Phi.shape, "\n", ls_Phi.round(4))

# 水平集函数值 Phi 从水平集节点投影到有限元节点
from scipy.interpolate import griddata
fe_Phi = griddata(points, ls_Phi, xi, method='cubic')
print("fe_Phi:", fe_Phi.shape, "\n", fe_Phi.round(4))

from scipy.io import loadmat
mat = loadmat('fe_phi_max.mat')
data = mat['FENd'][0, 0]
phi_data = data['Phi'].astype(np.float64)
fe_Phi_mat = phi_data[:, 0]
print("fe_Phi_mat:", fe_Phi_mat.shape, "\n", fe_Phi_mat.round(4))

error_2 = np.sum(np.abs(fe_Phi - fe_Phi_mat))
print("error_2:", error_2)


