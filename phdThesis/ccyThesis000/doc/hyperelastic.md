# 超弹性模型及其有限元方法模拟
超弹性材料是一类具有特殊力学行为的材料，
它们在受力后能够表现出非常大的弹性变形而不发生永久性变形或损坏。
这种性质使得超弹性材料在工程应用中具有广泛的应用前景，
例如在弹性垫片、胶体材料、橡胶制品、生物医学器械等领域。

## 一、超弹性模型
考虑一块材料的变形过程, 其在 $t$ 时刻的位置为 $\Omega(t)$, 将 $\Omega(0)$ 称为参考构型，
当前 $t$ 时刻 $\Omega(t)$ 称为当前构型。参考构型上的每个物质点 $\boldsymbol{X}$
对应于 $\Omega(t)$ 上的一个物质点 $\boldsymbol{x}$，那么存在一个变形函数:
$$
\boldsymbol{x} = \boldsymbol{\varphi}(\boldsymbol{X}, t)
$$
假设变形过程中没有发生错位断裂，那么该映射是一个参考构型到当前构型的微分同胚。

### （1）变形几何中的概念
记参考构型上的标架为 $\{\boldsymbol{e}_i\}_{i=0}^2$, 当前构型标架为 $\{\boldsymbol{e}_i'\}_{i=0}^2$. 定义材料位移为 
$$
\boldsymbol{u} = \boldsymbol{x-X}
$$
变形梯度为
$$
\boldsymbol{F} = \frac{\mathrm{d}\boldsymbol{x}}{\mathrm{d}\boldsymbol{X}} =
\frac{\partial x_i}{\partial X_j}\boldsymbol{e}_i'\otimes \boldsymbol{e}_j =
\nabla_{\boldsymbol{X}}\boldsymbol{u} - \boldsymbol{I}
$$
体积比 
$$
J = \mathrm{det}(\boldsymbol{F})
$$
Cauchy-Green 变形张量:
$$
\boldsymbol{C} = \boldsymbol{F}^T \boldsymbol{F} = 
F_{ik}F_{kj} \boldsymbol{e}_i\otimes \boldsymbol{e}_j
$$
以上概念仅与变形有关，在连续介质力学中还有一个非常重要的概念: 应力. 
应力不仅与变形有关还和和材料本身的性质有关。连续介质力学中有很多种应力概念，
本文使用的应力为第一类 PK 应力 
$$
\boldsymbol{P} = P_{ij}\boldsymbol{e}_i'\otimes \boldsymbol{e}_j
$$ 

<!-- 
和第二类 PK 应力 
$$
\boldsymbol{S} = S_{ij}\boldsymbol{e}_i\otimes \boldsymbol{e}_j
$$ 
他们之间的关系为:
$$
\boldsymbol{S} = \boldsymbol{F}^{-1}\boldsymbol{P} 
$$
-->

### (2) 平衡方程
考虑参考构型上的平衡方程:
$$
-\mathbf{div}(\boldsymbol{P}) = \boldsymbol{f}
$$
其中:$\mathbf{div}(\boldsymbol{P}) = \partial_{x_i}P_{ij} \boldsymbol{e}_{j}$,
$\boldsymbol{f}$ 为体力. 边界条件为:
$$
\begin{aligned}
\boldsymbol{u} & = \boldsymbol{g}_D\quad  \quad \mathrm{on}\ \ \Gamma_D\\
\boldsymbol{P}\cdot \boldsymbol{n} & = \boldsymbol{g}_N\quad  \quad \mathrm{on}\ \ \Gamma_N
\end{aligned}
$$
其中 $\Gamma_D \cup \Gamma_N = \partial \Omega$, $\boldsymbol{n}$ 为外法向。

平衡方程就是本文要求解的方程，本文以位移 $\bm u$
为未知量，需要将方程修改为 $\bm u$ 的方程，即找到
$\bm P$ 和 $\bm u$ 之间的关系，这个在连续介质力学中就是本构关系。

### (3) 本构方程
本构关系描述的是一个材料发生形变(位移)以后产生的应力有多大。
超弹性材料的一个特点就是本构关系可以由一个应变能泛函 $W$ 描述，且满足:
$$
\boldsymbol{P} = \frac{\partial W}{\partial \boldsymbol{F}}
$$
应变能泛函有很多种，适用于不同的情况，本文使用常用的 Neo-Hookean 应变能:
$$
W = \frac{\mu}{2}(\mathrm{tr}(\boldsymbol{C})-3) - \mu\ln(J) + \frac{\lambda}{2}(\ln(J))^2
$$
其中 $\mathrm{tr}(\boldsymbol{C}) = \sum_iC_{ii}$, $\mu, \lambda$ 是材料的参数，
不同材料有不同的值。那么第一类 PK 应力就可以被计算:
$$
\boldsymbol{P} = \frac{\partial W}{\partial \boldsymbol{F}} = 
\mu(\boldsymbol{F}-\boldsymbol{F}^{-T})+\lambda\ln(J)\boldsymbol{F}^{-T}
$$
详细过程见附录。注意上式中 $\bm F = \nabla\bm u - \bm I$, 所以 $\bm P$ 和 $\bm u$ 之间的关系被找到了。

## 二、有限元求解
本节将会考虑 $\boldsymbol{u} \in (H^1(\Omega))^d$,对问题进行变分得到可用有限元离散的格式。

根据前面的讨论，$-\mathbf{div}(\boldsymbol{P})$ 可以认为是 $(H^1(\Omega))^d$
到$(H^{-1}(\Omega))^d$ 上的一个非线性二阶微分算子, 其中
$-\mathbf{div}$是一个线性算子，非线性的来源是 $\boldsymbol{P}$，$\boldsymbol{P}$
可以看做是从 $(H^1(\Omega))^d$ 到$(L^2(\Omega))^{d\times d}$ 的非线性算子.
那么问题的求解需要将 $\boldsymbol{P}$ 线性化, 使用的工具就是变分。

### （1）变分的数学解释
对于两个 Banaha 空间 $A, B$ 之间的一个算子 $\mathcal L$, 记 $\mathcal L_{u}$
为$\mathcal L$ 在 $u$ 处的变分, 那么上一节中出现的变量 $\boldsymbol{F}$, $J,
\boldsymbol{F^{-T}}$, $\boldsymbol{P}$ 都可以认为是 $(H^1(\Omega))^d$
到其他空间的算子。令 $\delta \boldsymbol{u}$ 为位移的一个增量，对于
$(H^1(\Omega))^d$ 上的一个算子 $\mathcal L$, 记: 
$$
\delta \mathcal L := \mathcal{L}(\boldsymbol{u}+\delta \boldsymbol u) -
\mathcal L(\boldsymbol{u}) = \mathcal L_u(\delta \boldsymbol{u}) + o\|\delta \bm
u\|^2
$$
那么有:
$$
\delta \boldsymbol{F} = \nabla\delta \boldsymbol{u}, \quad
\delta J = 
J \boldsymbol{F}^{-T}:\delta \boldsymbol{F}, \quad 
\delta \boldsymbol{F}^{-T} = \boldsymbol{F}^{-T}\delta 
\boldsymbol{F} \boldsymbol{F}^{-T}
$$
$$
\delta \boldsymbol{P} = \mu(\delta \boldsymbol{F} - 
\boldsymbol{F}^{-T}\delta \boldsymbol{F} \boldsymbol{F}^{-T}) +
\lambda(
\ln J \boldsymbol{F}^{-T}\delta \boldsymbol{F} \boldsymbol{F}^{-T} + 
(\boldsymbol{F}^{-T}: \delta\boldsymbol{F})\boldsymbol{F}^{-T}
)+o\|\delta \bm u\|^2
$$
详细计算过程见附录。
### （2）弱形式
定义 $H_D^1(\Omega) = \{v\in H^1(\Omega) : v|_{\Gamma_D} = 0\}$
将平衡方程两边同时与 $\boldsymbol{v} \in (H_D^1(\Omega))^d$ 做内积:
$$(-\mathbf{div}(\boldsymbol{P}), \boldsymbol{v}) = (\boldsymbol{f}, \boldsymbol{v})$$
分部积分后得到:
$$(\boldsymbol{P}, \nabla \boldsymbol{v}) = 
(\boldsymbol{f}, \boldsymbol{v})+
(\boldsymbol{g}_N, \boldsymbol{v})_{\Gamma_N}$$
定义 $(H^1(\Omega))^d\times(H_D^1(\Omega))^d$ 上的非线性泛函:
$$
a(\boldsymbol{u}, \boldsymbol{v}) = (\boldsymbol{P}, \nabla \boldsymbol{v})
$$
定义 $(H^1_D(\Omega))^d$ 上的线性型:
$$b(\boldsymbol{v}) = 
(\boldsymbol{f}, \boldsymbol{v})+
(\boldsymbol{g}_N, \boldsymbol{v})_{\Gamma_N}$$
非线性变分问题为: 找到 $\boldsymbol{u} \in (H^1(\Omega))^d$ 满足:
$$
a(\boldsymbol{u}, \boldsymbol{v}) = b(\boldsymbol{v}), \quad\quad \forall 
\boldsymbol{v} \in (H^1_D(\Omega))^d \tag{1}
$$
### （3）迭代算法
方程 (1) 是一个非线性问题，需要线性化后迭代求解。
任取 $\boldsymbol{u}_0\in H^1(\Omega)$，设 
$\delta \boldsymbol{u}_1 = \boldsymbol{u-u}_0$，那么有:
$$
a(\boldsymbol{u}, \boldsymbol{v}) = a(\boldsymbol{u}_0, \boldsymbol{v}) +
(\delta \bm P_1, \nabla \bm v)
$$
写成另一种形式:
$$
(\delta \bm P_1, \nabla \bm v) = 
a(\boldsymbol{u}, \boldsymbol{v})-a(\boldsymbol{u}_0, \boldsymbol{v})
 = b(\bm v) - a(\bm u_0, \bm v) \tag{2}
$$
删去 $\delta\bm P_1$ 中的高阶项, 并定义:
$$
a_{\bm u_0}(\delta \bm u_1, \bm v) := (\bm P_{\bm u_0}(\delta \bm u_1), \nabla \bm v) \approx (\delta \bm P_1, \nabla \bm v)
$$
那么 $a_{\bm u_0}$ 是一个双线性泛函，求解与方程 (2) 相似的线性问题: 找到 $\delta\boldsymbol{u}_1
\in (H^1_D(\Omega))^d$ 满足:
$$
a_{\boldsymbol{u}_0}(\delta\boldsymbol{u}_1, \boldsymbol{v}) = 
b(\boldsymbol{v}) - a(\boldsymbol{u}_0, \boldsymbol{v}), 
\quad\quad \forall \boldsymbol{v} \in (H^1_D(\Omega))^d
$$
得到 $\delta \boldsymbol{u}_1$ 后，令 $\boldsymbol{u}_1 =\boldsymbol{u}_0+\delta
\boldsymbol{u}_0$, 这是一个 $\bm u$ 的近似。重复这个过程得到一个函数列 
$\{\boldsymbol{u}_k\}$, 当 $a$ 满足某些性质时函数列收敛于真解。
### （4）有限元离散
使用单纯形网格 $\mathcal T_h$ 离散区域 $\Omega$, 在 $\mathcal T_h$ 上建立
Lagrange 有限元空间 $V_h$ 离散
$(H^1(\Omega))^d$，记 
$$
V_{h, D} = (H^1_D(\Omega))^d \cap V_h
$$ 
使用上述迭代算法求解超弹性模型，选定初值 $\bm u_h^0 = 0$，
迭代过程中假设我们已经知道 $\bm u_h^k$，我们需要求解的有限元问题为: 找到
$\delta\boldsymbol{u}_h^k \in V_h$ 满足:
$$
a_{\boldsymbol{u}_h^{k-1}}(\delta\boldsymbol{u}_h^k, \boldsymbol{v}) = 
b(\boldsymbol{v}) - a(\boldsymbol{u}_h^{k-1}, \boldsymbol{v}),  
\quad\quad \forall \boldsymbol{v} \in V_{h, D}
$$
现在假设

</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
