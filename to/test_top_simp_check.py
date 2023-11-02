import numpy as np
import matplotlib.pyplot as plt

from top_simp import TopSimp


nelx, nely = 60, 20
rmin = 1.5
volfrac = 0.5
x = np.full((nely, nelx), volfrac)

# 创建一个随机灵敏度分布
np.random.seed(0)
dc = np.random.rand(nely, nelx)

# 创建 TopSimp 对象
tsp = TopSimp()

# 调用 check 函数
dcn = tsp.check(nelx, nely, rmin, x, dc)

# 绘制原始和过滤后的灵敏度分布
fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# 原始灵敏度
im1 = ax[0].imshow(dc, cmap='viridis')
ax[0].set_title('Original Sensitivities')
fig.colorbar(im1, ax=ax[0], orientation='vertical', fraction=0.046, pad=0.04)

# 过滤后的灵敏度
im2 = ax[1].imshow(dcn, cmap='viridis')
ax[1].set_title('Filtered Sensitivities')
fig.colorbar(im2, ax=ax[1], orientation='vertical', fraction=0.046, pad=0.04)

plt.show()
