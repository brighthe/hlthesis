import numpy as np
import matplotlib.pyplot as plt
import matplotlib
def  velocity_field(p):
    x = p[..., 0]
    y = p[..., 1]
    u = np.zeros(p.shape)
    u[..., 0] = np.sin((np.pi*x))**2 * np.sin(2*np.pi*y)
    u[..., 1] = -np.sin((np.pi*y))**2 * np.sin(2*np.pi*x)
    return u

# 生成网格点
x = np.linspace(0, 1, 60)
y = np.linspace(0, 1, 60)
xx, yy = np.meshgrid(x, y)
# h1 = x[0] - y[0]

# 对速度场进行采样
points = np.stack([xx, yy], axis=-1)
u = velocity_field(points)
plt.quiver(x,y,u[...,0],u[...,1])
plt.show()
# 为了得到大小，我们需要计算速度矢量的模
v = np.linalg.norm(u, axis=-1)

def circle(p):
    x = p[...,0]
    y = p[...,1]
    val = np.sqrt((x-0.5)**2+(y-0.75)**2)-0.15
    return val
    
# 生成网格点
x = np.linspace(0, 1, 62)
y = np.linspace(0, 1, 62)
xx, yy = np.meshgrid(x, y)
h2 = x[1] - x[0]

# 对circle函数进行采样
points = np.stack([xx, yy], axis=-1)
lsf = circle(points)


def evolve(v, lsf):
    vFull = np.pad(v, ((1,1),(1,1)), mode='constant', constant_values=0)
    
    dt = 0.001

    plt.figure(figsize=(1, 1))
    
    for i in range(1000):
        dpx = ( np.roll(lsf, shift=(0, -1), axis=(0, 1)) - lsf ) / h2 # 向前差分
        dmx = ( lsf - np.roll(lsf, shift=(0, 1), axis=(0, 1)) ) / h2  # 向后差分
        dpy = ( np.roll(lsf, shift=(-1, 0), axis=(0, 1)) - lsf ) / h2
        dmy = ( lsf - np.roll(lsf, shift=(1, 0), axis=(0, 1)) ) / h2
        
        lsf = lsf - dt * np.minimum(vFull, 0) * np.sqrt( np.minimum(dmx, 0)**2 + np.maximum(dpx, 0)**2 + np.minimum(dmy, 0)**2 + np.maximum(dpy, 0)**2 ) \
                  - dt * np.maximum(vFull, 0) * np.sqrt( np.maximum(dmx, 0)**2 + np.minimum(dpx, 0)**2 + np.maximum(dmy, 0)**2 + np.minimum(dpy, 0)**2 )
        
        # 每隔若干迭代绘制一次结果
        if i % 1 == 0:  
            plt.imshow(lsf, cmap='hot', origin='lower')
            plt.colorbar()
            plt.title(f"Iteration {i}")
            plt.pause(0.01)  # 暂停图像0.1秒，
            # plt.show()
            plt.clf()  # 清除当前图像

    plt.show()  

evolve(v=v, lsf=lsf)
