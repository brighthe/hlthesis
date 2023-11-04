import numpy as np
import matplotlib.pyplot as plt
from top_simp import TopSimp

# 设定参数
nelx, nely = 60, 20
volfrac = 0.5
x = np.full((nely, nelx), volfrac)

# 创建一个随机灵敏度分布
np.random.seed(0)
dc = -np.random.rand(nely, nelx) # 在实际的拓扑优化问题中，dc（灵敏度）应该是负数，因为我们通常希望最小化目标函数

# 创建 TopSimp 对象
tsp = TopSimp()

# 调用 OC 函数
xnew = tsp.OC(nelx, nely, volfrac, x, dc)

# 绘制原始和更新后的设计变量分布
fig, ax = plt.subplots(1, 2, figsize=(10, 5))
ax[0].imshow(x, cmap='viridis')
ax[0].set_title('Original Design Variables')
ax[1].imshow(xnew, cmap='viridis')
ax[1].set_title('Updated Design Variables')

# 添加颜色条以便于比较
cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
fig.colorbar(ax[1].get_images()[0], cax=cbar_ax)

plt.show()
