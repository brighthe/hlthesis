import numpy as np

from simp_top_cantilever_fixed_hole import TopSimp

# Cantilever beam with a fixed hole
nelx = 45
nely = 30
volfrac = 0.5
penal = 3.0
rmin = 1.5
ts = TopSimp(nelx=nelx, nely=nely, volfrac=volfrac, penal=penal, rmin=rmin)

# 初始化优化参数
nelx, nely, volfrac, penal, rmin = ts._nelx, ts._nely, ts._volfrac, ts._penal, ts._rmin
mesh = ts._mesh_top2

node = mesh.entity('node') # 按列增加
cell = mesh.entity('cell') # 左下角逆时针

# 根据体积分数 volfrac 初始化设计变量场
x = np.full((nely, nelx), volfrac)

# 利用 passive 单元给结构打孔洞
passive = np.zeros((nely, nelx))
for ely in range(nely):
    for elx in range(nelx):
        if np.sqrt((ely+1 - nely/2.)**2 + (elx+1 - nelx/3.)**2) < nely / 3.:
            passive[ely, elx] = 1
            x[ely, elx] = 0.001

# 初始化每个单元的目标函数的灵敏度为 0
dc = np.zeros((nely, nelx))

from mbb_beam_operator_integrator import MbbBeamOperatorIntegrator
E0 = 1.0
nu = 0.3
integrator = MbbBeamOperatorIntegrator(nu=nu, E0=E0, nelx=nelx, nely=nely, penal=penal, x=x)
KE = integrator.stiff_matrix()

loop = 0 # 迭代次数
change = 1.0 # 迭代中设计变量的最大变化

# 优化循环，直到 change < 0.01 迭代终止
while change > 0.01:
    loop += 1
    xold = np.copy(x)

    # 对当前设计执行有限元分析计算位移
    U, Ue = ts.FE(mesh=mesh, x=x, penal=penal)

    # 计算目标函数
    c = 0 # 初始化目标函数为 0
    temp1 = x ** penal
    temp2 = np.einsum('ij, jk, ki -> i', Ue, KE, Ue.T).reshape(nelx, nely).T
    c = c + np.einsum('ij, ij -> ', temp1, temp2)

    # 计算每个单元的目标函数的灵敏度
    temp3 = -penal * x ** (penal-1)
    dc[:] = np.einsum('ij, ij -> ij', temp3, temp2)

    # 灵敏度过滤
    dc = ts.check(rmin=rmin, x=x, dc=dc)

    # 使用 Optimality Criteria Method 更新设计
    x = ts.OC(volfrac=volfrac, x=x, dc=dc, passive=passive)

    # 计算当前的 volume fraction
    volfrac = np.sum(x) / (nelx*nely)

    # 打印当前迭代的结果
    change = np.max(np.abs(x - xold))
    print(f' Iter.: {loop:4d} Objective.: {c:10.4f} Volfrac.: {volfrac:6.3f} change.: {change:6.3f}')
    
    # 可视化材料密度分布
    import matplotlib.pyplot as plt
    plt.imshow(-x, cmap='gray')
    plt.axis('off')
    plt.axis('equal')
    plt.draw()
    plt.pause(1e-5)
        
plt.ioff()
plt.show()
