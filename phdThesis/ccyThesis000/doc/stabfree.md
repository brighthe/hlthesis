在本文中，我们首先对如何去除虚单元方法中稳定化项进行了讨论，其关键是找到一个
多面体 $K$ 上的有限维空间 $\mathbb{V}_{k-1}(K)$，
一个虚单元空间 $V_k$ 的梯度空间到 $\mathbb{V}_{k-1}(K)$ 的投影算子 $Q_K$ 满足:
1. $\|Q_K\nabla v\|_{0,K}=\|\nabla v\|_{0,K}\quad\forall v\in V_k(K)$
2. $Q_K \nabla v$ 可以根据 $v$ 在 $V_k$ 中的自由度进行计算。

根据这个要求，我们将多边形 $K$ 分成单纯形网格 $\mathcal{T}_K$，
在 $\mathcal{T}_K$ 上构造了宏元空间 $\mathbb{V}_{k-1}(K)$:
$$
\mathbb{V}_{k-1}(K):=\{\phi\in\boldsymbol{H}(\mathrm{div},K):
\phi|_T\in\mathbb{P}_{k-1}(T;\mathbb{R}^d)
\text{ for each }
T\in\mathcal{T}_K, \mathrm{div}\phi \in \mathbb{P}_{k-2}(K)\}
$$
我们证明了：
$$
\mathbb{V}_{k-1}(K)=\mathrm{div}\ \mathrm{skw}
V_{k}^{d-2}(K)\oplus(\boldsymbol{x}-\boldsymbol{x}_{K})\mathbb{P}_{\max\{k-2,0\}}(K).
$$
其中 $\mathrm{div}\ \mathrm{skw} V_{k}^{d-2}(K)$ 当 $d=2$ 时， 是 Lagrange 
元空间的旋度空间，当 $d=3$ 时，是 Nedelec 元空间的旋度空间。
对于 $\phi\in \mathbb{V}_{k-1}(K)$ 我们证明了如下范数等价：
$$
\|\phi\|_{0,K}\simeq
h_{K}\|\operatorname{div}\phi\|_{0,K}+\sup_{\psi\in\operatorname{div}\hat{\boldsymbol{V}}_{k}^{d-2}(K)}\frac{(\phi,\boldsymbol{\psi})_{K}}{\|\boldsymbol{\psi}\|_{0,K}}+\sum_{F\in\mathcal{F}^{\partial}(\mathcal{T}_{K})}h_{F}^{1/2}\|\phi\cdot\boldsymbol{n}\|_{0,F}.
$$
我们提出了二阶椭圆方程的
$H^1$协调、非协调虚单元方法，其与传统虚单元方法相比，区别是将梯度投影到了 $\mathbb{V}_{k-1}(K)$ 上，去掉了稳定化项。
根据上面的范数等价，我们可以证明了提出的方法的稳定性和误差估计。
最后我们给出了数值算例，验证了算法的收敛性，稳定性。



