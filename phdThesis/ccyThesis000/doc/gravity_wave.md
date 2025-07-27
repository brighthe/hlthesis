# 引力波学习文档
广义相对论是天文学中的一个重要理论，它将宇宙描述为一个四维的空间时间流形，
具有弯曲的度量，其曲率描述了引力场的本质。
## 宇宙的几何描述
宇宙是一个四维的空间时间流形 $M$，配备了一个 Lorentz 度量 $g_{ab}$，
$$
g_{ab} = g_{\mu\nu}(dx^{\mu})_a\otimes(dx^{\nu})_b
$$
其中 $g_{\mu\nu}$ 是一个对称非退化矩阵，$\{(dx^{\mu})_a\}_{\mu=0}^3$ 
是 $M$ 上的余切向量基。最简单的 Lorentz 度量是 Minkowski 度量 $\eta_{ab}$:
$$
(\eta_{\mu\nu}) = \mathrm{diag}(-1, 1, 1, 1)
$$
对于向量 $v^a \in T_pM$, 若 
- $g_{ab}v^av^b = 0$，则称 $v^a$ 是类光向量
- $g_{ab}v^av^b > 0$，则称 $v^a$ 是类空向量
- $g_{ab}v^av^b < 0$，则称 $v^a$ 是类时向量

如果一个超曲面 $\Sigma$ 的法向量 $n^a$ 是类时向量，那么称 $\Sigma$ 是类空超曲面。
此时度量 $g_{ab}$ 在 $\Sigma$ 上的限制是一个 Riemann 度量。

## Einstein 方程
Einstein 方程是广义相对论的核心, 
其描述了空间中质量和能量的分布如何影响时空的几何形状:
$$
G_{ab} = \frac{8\pi G}{c^4}T_{ab} 
$$
其中 $G_{ab} = R_{ab} - \frac{1}{2}g_{ab}R$ 是 Einstein 张量，
$T_{ab}$ 是能量动量张量，$R_{ab}$ 是 Ricci 张量，$R$ 是标量曲率。
在真空中 $T_{ab} = 0$，Einstein 方程变为:
$$
G_{ab} = 0
$$
显然 $\mathrm{tr}(G) = 0$，即 
$g^{ab}G_{ab} = g^{ab}R_{ab} - \frac{1}{2}g^{ab}g_{ab}R = 0$，
所以 $R = 0, R_{ab} = 0$。反过来，如果 Ricci
张量为零，那么 Einstein 张量也为零。所以真空中的 Einstein 方程等价于:
$$
R_{ab} = 0
$$
## 线性化 Einstein 方程
假设 $g_{ab}$ 是 Minkowski 度量的微小扰动:
$$
g_{ab} = \eta_{ab} + \epsilon h_{ab} + O(\epsilon^2)
$$
其中 $\eta_{ab}$ 是 Minkowski 度量，$0<\epsilon\ll 1$，未知量变成了 $h_{ab}$。
将上面的度量代入 Einstein 方程，令 $\nabla$ 是 $\eta_{ab}$ 的联络，
那么与 $g_{ab}$ 适配的联络与 $\nabla$ 之间的克氏符为:
$$
\Gamma_{abc} = \frac{1}{2}(\nabla_ag_{bc} + \nabla_bg_{ac} - \nabla_cg_{ab})
= \epsilon\frac{1}{2}(\nabla_ah_{bc} + \nabla_bh_{ac} - \nabla_ch_{ab}) + O(\epsilon^2)
$$
$\Gamma_{abc} = O(\epsilon)$，所以 Riemann 张量为:
$$
\begin{aligned}
R_{abcd} & = \nabla_c\Gamma_{bda} - \nabla_d\Gamma_{bca} + O(\epsilon^2)\\
& = \frac{1}{2} \epsilon(\nabla_c\nabla_bh_{da} - \nabla_c\nabla_ah_{bd} -
\nabla_d\nabla_bh_{ca} + \nabla_d\nabla_ah_{bc}) + O(\epsilon^2)\\
\end{aligned}
$$
所以 Ricci 张量为:
$$
\begin{aligned}
R_{bd} = R^{c}_{\ bcd} & = \frac{1}{2}\epsilon(\nabla_c\nabla_bh^{c}_{\ d} 
- \nabla_c\nabla^ch_{bd} - \nabla_d\nabla_bh^{c}_{\ c} + \nabla_d\nabla^ch_{bc}) + O(\epsilon^2)\\
& = \epsilon\nabla_c\nabla_{(b}h^{\ c}_{d)} - 
\frac{1}{2}\epsilon(\nabla_c\nabla^ch_{bd} + \nabla_d\nabla_b h^c_{\ c}) + 
O(\epsilon^2)
\end{aligned}
$$
$$
R = R^{a}_{\ a} = \epsilon(\nabla_c\nabla^ah^{\ c}_{a} - \nabla_c\nabla^ch^{a}_{\ a}) + O(\epsilon^2)
$$
$R = O(\epsilon)$ 那么 Einstein 张量为:
$$
\begin{aligned}
G_{ab} & = R_{ab} - \frac{1}{2}g_{ab}R\\
& = \epsilon\nabla_c\nabla_{(a}h^{\ c}_{b)} 
- \frac{1}{2}\epsilon(\nabla_c\nabla^ch_{ab} 
+ \nabla_b\nabla_a h^c_{\ c}) - \frac{1}{2}\epsilon\eta_{ab}( \nabla_c\nabla^dh^{c}_{\ d} - \nabla_c\nabla^ch^{d}_{\ d})
+ O(\epsilon^2)\\ 
\end{aligned}
$$
令 $h = h^a_a, \bar{h}_{ab} = h_{ab} - \eta_{ab} h/2$ 那么 $\bar{h} = -h, 
h_{ab} = \bar{h}_{ab} + \eta_{ab}\bar{h}/2$，所以:
$$
\nabla_{a}\nabla_{b} h_{cd} = \nabla_{a}\nabla_{b} \bar{h}_{cd} + \frac{1}{2}\eta_{ab}\nabla_{c}\nabla_{d}\bar{h}
$$
所以有:
$$
\frac{\mathrm{d}}{\mathrm{d}\epsilon}G_{ab} = \nabla^c\nabla_{(a}\bar{h}_{b)c} - 
\frac{1}{2}(\eta_{ab}\nabla^c\nabla^d\bar{h}_{cd} + \nabla^c\nabla_c\bar{h}_{ab})
=: H(\bar{h}_{ab})
$$
显然 $H$ 是一个线性算子，所以 Einstein 方程的线性化形式为:
$$
H(\bar{h}_{ab}) = \frac{8\pi G}{c^4}\mathcal{T}_{ab}
$$
其中 $\mathcal{T}_{ab} = \frac{\mathrm{d}}{\mathrm{d}\epsilon}T_{ab}$。
需要注意的是，上式保证了 $\nabla_a\mathcal{T}^{ab} = 0$，所以线性化的 Einstein 方程
可以保证能量，动量，角动量守恒。

## 规范自由性
假设 $\phi$ 是两个微分流形 $M, N$ 之间的微分同胚，
那么 $M$ 上的度量 $g$ 可以推前到 $N$，$(N, \phi_*g)$ 与 $(M, g)$ 
等距同构。那么如果 $g$ 是 Einstein
方程的解，那么由于:
$$
R(\phi_*g) = \phi_*R(g)
$$
所以 $\phi_*g$ 也是 Einstein 方程的解。也是 Einstein 方程的解，他们互相等价，称为规范相关。
对于线性化的 Einstein 方程，任取向量场 $v^a$，定义 $h'_{ab} = \mathcal{L}_v\beta_{ab} = \nabla_{(a}v_{b)}$，
那么有:
$$
H(h'_{ab}) = 0 
$$
即 $h'_{ab}\in \ker(H)$，所以 Einstein 方程的解不是唯一的，
这个问题不能由边界条件来解决。我们可以施加一个规范条件，在等价的解中选择一个。
类似的情况在电磁学中也有，当我们求解磁矢势 $A$ 时：
$$
\boldsymbol{B} = \nabla \times \boldsymbol{A}
$$
这个方程的解不是唯一的，当 $\boldsymbol{A}'\in\ker(\nabla\times)$, 
$\boldsymbol{A+A}'$ 也是方程的解。
这时可以加上一个条件：
$$
\nabla\cdot\boldsymbol{A} = 0
$$
这个条件就是规范条件。在线性 Einstein 方程中，常用的规范条件是：
$$
\nabla^b\bar{h}_{ab} = 0
$$
在这个规范条件下，线性化的 Einstein 方程的形式更加简单:
$$
H(\bar{h}_{ab}) = \square\bar{h}_{ab} = \frac{8\pi G}{c^4}\mathcal{T}_{ab} 
$$
其中 $\square = \nabla^a\nabla_a$ 是 d'Alembert 算子。
> **注1:** 对于任意 $h_{ab}$，选择 $v^a$ 满足 $\nabla^b\nabla_b v^a = -\nabla^b h_{ab}$，
> 令 $h'_{ab} = h_{ab}+\nabla_{(a}v_{b)}$ 那么有:
> $$
> \begin{aligned}
> \bar{h'}_{ab} 
> & = h'_{ab} - \frac{1}{2}\eta_{ab}h'\\
> & = h_{ab} + \nabla_{(a}v_{b)} - \frac{1}{2}\eta_{ab}(h_{c}^{\ c} + 
> \eta^{cd}\nabla_{c}v_{d})\\
> & = \bar{h}_{ab} + \nabla_{(a}v_{b)} - \frac{1}{2}\eta_{ab}\nabla^{c}v_{c}\\
> \end{aligned}
> $$
> 所以:
> $$
> \begin{aligned}
> \nabla^b \bar{h'}_{ab} & = \nabla^b\bar{h}_{ab} + \nabla^b\nabla_{(a}v_{b)} - \frac{1}{2}\nabla_a\nabla^bv_b\\
> & = \nabla^b\bar{h}_{ab} + \frac{1}{2}\nabla_b\nabla^bv_a\\
> & = 0
> \end{aligned}
> $$
> 所以，对于任意 $h_{ab}$，总可以找到一个等价的 $h'_{ab}$ 满足规范条件。

> **注2:** 关于规范相关有一个等价的解释：选择一个合适的坐标系，
>  Einstein 方程可以写成分量形式:
> $$
> R_{\mu\nu}(x) - \frac{1}{2}g_{\mu\nu}(x)R(x) = \frac{8\pi G}{c^4}T_{\mu\nu}(x)
> $$
> 其中每个分量都是 $x$ 的函数。这个方程的未知量是 $g_{\mu\nu}(x)$，因为 $g$ 对称,
> 所以有 $10$ 个未知量。选择适当的边界条件，上式同时也是 10 个偏微分方程，
> 原本非常合理，10 个方程求解 10 个未知量，但是因为有 Bianchi 恒等式:
> $$
> \nabla_{[a}R_{bc]d}^{\quad\ e} = 0
> $$
> 这可以推出 $\nabla_aG^{ab} = 0$, 其分量形式是四个方程:
> $$
> \nabla_{\mu}G^{\mu\nu} = 0
> $$
> 这表明 10 个方程中只有 6 个是独立的，所以实际上只有 6 个方程求解 10 个未知量。
> 实际上，如果 $g_{\mu\nu}(x)$ 是 Einstein 方程的解，那么 $g_{\mu\nu}(x)$ 的任意坐标变换
> $g'_{\mu\nu}(x)$ 也是 Einstein 方程的解，那么 $g_{\mu\nu}(x)$ 和 $g'_{\mu\nu}(x)$ 是
> 规范相关。
>
> 最后添加的规范条件就是为了消除这种规范相关性，使得解是唯一的。规范条件的分量形式是:
> $$
> \nabla^{\mu}\bar{h}_{\mu\nu} = 0
> $$
> 这是四个方程，所以加上规范条件后，10 个方程求解 10 个未知量。

## 线性 Einstein 方程的平面波解
令 $k^a$ 是一个类光矢量，即 $k^ak_a = 0$，$A_{ab}$ 是一个常对称张量，即 
$\nabla_cA_{ab} = 0$。
对于任意的 $C^2$ 函数 $f(x)$，定义:
$$
\bar{h}_{ab} = A_{ab}f(k_cx^c)
$$
那么:
$$
\square \bar{h}_{ab} = A_{ab}\square f(k_cx^c) = A_{ab}k^ck_c f''(k_cx^c) = 0
$$
而且：
$$
\nabla^b\bar{h}_{ab} = A_{ab}\nabla^bf(k_cx^c) = A_{ab}k^bf'(k_cx^c) = 0
$$
所以 $\bar{h}_{ab}$ 是线性化 Einstein 方程的解。这种解称为平面波解。

在平面波解中，我们可以转化原始的规范条件为：
1. $A^a_{\ a} = 0$
2. 对于某个类时矢量 $u^b$, $A^{ab}u^b = 0$.

这个条件被称为 TT 规范条件。
选择 $u^b = (1, 0, 0, 0)$，那么前面定义的平面波解 $\bar{h}_{ab}$ 可以找打一个
规范相关的解 
$\bar{h}'_{ab} = A'_{ab}f(k_cx^c)$ 满足 TT 规范条件。

<!--
## 李导数
流形上一个单参微分同胚群 $\phi : M\times\mathbb{R}\to M$ 满足:
- $\phi_t$ 是 $M$ 上的微分同胚
- $\phi_{t+s} = \phi_t\circ\phi_s$

**注:**$\phi$ 是两个参数的映射，固定 $t$，上面第一个性质表明 
$\phi_t$ 是 $M$ 上的微分同胚，固定 $p$，上面第二个性质表明
$\phi_p$ 是 $\mathbb{R}$ 到 $M$ 的映射, 其像是一个曲线，而且 
$\forall q \in \mathcal{R}(\phi_p)$, $\mathcal{R}(\phi_q) = \mathcal{R}(\phi_{p})$。
一个向量场 $\boldsymbol{X}$ 可以定义一个单参微分同胚群 $\phi$: $\phi_t(p) = \exp(t\boldsymbol{X})p$，

对于向量场 $\boldsymbol{X}$ 定义的单参微分同胚群，一个张量 $\boldsymbol{T}$ 的李导数 $\mathcal{L}_{\boldsymbol{X}}\boldsymbol{T}$ 定义为:
$$
\mathcal{L}_{\boldsymbol{X}}\boldsymbol{T} = \lim_{t\to 0}\frac{1}{t}((\phi_{-t})^*\boldsymbol{T} - \boldsymbol{T}) = \lim_{t\to 0}\frac{1}{t}((\phi_{t})_*\boldsymbol{T} - \boldsymbol{T})
$$
当 $\boldsymbol{T}$ 是一个向量场时:
$$
\mathcal{L}_{\boldsymbol{X}}\boldsymbol{Y} = [\boldsymbol{X}, \boldsymbol{Y}]
= \nabla_{\boldsymbol{X}}\boldsymbol{Y} - \nabla_{\boldsymbol{Y}}\boldsymbol{X}
$$
对于一个度量张量 $g$:
$$
\mathcal{L}_{\boldsymbol{X}}g = \nabla_{\boldsymbol{X}}g + \nabla_{g}\boldsymbol{X}
$$
-->


















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
</br>
</br>
</br>




