# $H(\mathrm{div}, \mathbb{S})$ 的对称张量虚单元
## 一、de Rham 复形的虚单元
$$
\mathbb{R} \xrightarrow{}H^1\xrightarrow{\nabla}H(\mathrm{curl})\xrightarrow{\mathrm{curl}}
H(\mathrm{div})\xrightarrow{\mathrm{div}} L^2 \xrightarrow{}0
$$

1. $H^1$ 协调虚单元空间
$$
V_k^{vert}(K) := \{v \in H^1(K) : v|_e \in \mathbb{P}_k(e)\ \forall e \in \partial K, \Delta v \in \mathbb{P}_{k-2}(K)\}
$$
2. $H(\mathrm{curl})$ 协调的虚单元空间
$$
\begin{aligned}
V_k^{edge}(K) := \{\bm v \in H({\mathrm{div}}, K)\cap H({\mathrm{curl}}, K)
: \mathrm{div}\bm v \in \mathbb{P}_{k-1}(K), \mathrm{curl}\bm v \in
\mathbb{P}_{k-1}(K) \\
\bm v \times \bm n \in \mathbb{P}_k(e) \ \forall e \in \partial K
\} 
\end{aligned}
$$

3. $H(\mathrm{div})$ 协调的虚单元空间
$$
\begin{aligned}
V_k^{face}(K) := \{\bm v \in H({\mathrm{div}}, K)\cap H({\mathrm{curl}}, K)
: \mathrm{div}\bm v \in \mathbb{P}_{k-1}(K), \mathrm{curl}\bm v \in
\mathbb{P}_{k-1}(K)\\
\bm v \cdot \bm n \in \mathbb{P}_k(e) \ \forall e \in \partial K
\} 
\end{aligned}
$$

4. $L^2$ 空间 : $\mathbb{P}_k(K)$

其中 $V_k^{face}$ 的存在性关键在于两个问题：
$$
\left\{
\begin{aligned}
\mathrm{rot}\bm v &= f \quad \Omega\\
\mathrm{div}\bm v &= g \quad \Omega\\
\bm v \cdot \bm n &= h \quad \partial\Omega
\end{aligned}
\right.
$$
对于$\Omega$上的任意光滑函数 $f, g$ 以及 $\partial\Omega$ 上的光滑函数 $h(\bm h)$ 
都存在唯一的 $v$.
## 二、有限元弹性复形
$$
\bm{RM}\xrightarrow{\subset}H^1(\Omega, \mathbb{R}^3)\xrightarrow{\mathrm{def}}
H(\mathrm{inc}, \Omega; \mathbb{S})\xrightarrow{\mathrm{inc}}H(\mathrm{div},
\Omega; \mathbb{S}) \xrightarrow{\mathrm{div}} L^2(\Omega, \mathbb{R}^3)\to 0
$$
其中 
$$
\mathrm{def}(\bm v) = \nabla \bm v + \nabla^{T} \bm v
$$
$$
\mathrm{inc}(\bm v) = \nabla\times\bm v \times\nabla
$$

在这个复形下有类似的问题：
$$
\left\{
\begin{aligned}
\mathrm{inc}\bm v &= f \quad \Omega\\
\mathrm{div}\bm v &= g \quad \Omega\\
\bm v \cdot \bm n &= h \quad \partial\Omega
\end{aligned}
\right.
$$
对于$\Omega$上的任意光滑函数 $f \in \mathrm{ker}(div), g$ 以及 $\partial\Omega$ 上的光滑函数 $h(\bm h)$ 
都存在唯一的 $v$ 满足上面三个等式。

**证明**: 将证明分两步：
1. $\mathrm{inc}\bm v = f$  这个存在性显然
2. 
    $$
    \left\{
    \begin{aligned}
    \mathrm{inc}\bm v &= f \quad \Omega\\
    \mathrm{div}\bm v &= g \quad \Omega\\
    \bm v \cdot \bm n &= h \quad \partial\Omega
    \end{aligned}
    \right.
    $$
    首先我们可以找到一个光滑函数 $\bm v_0$ 满足 $\mathrm{inc}\bm v_0 = f$, 令
    $\tilde{\bm g} = \bm g-\mathrm{div}\bm v_0，$$\tilde{\bm
    h} = \bm{h}-\bm{v}_0\cdot\bm{n}$那么问题可以修改为找到 $\bm v$ 满足:
    $$
    \left\{
    \begin{aligned}
    \mathrm{inc}(\bm v-\bm{v}_0)&= 0\quad \Omega\\
    \mathrm{div}(\bm v-\bm{v}_0)&= \tilde{g} \quad \Omega\\
    (\bm v-\bm{v}_0) \cdot \bm n &= \tilde{h} \quad \partial\Omega
    \end{aligned}
    \right.
    $$
    再把$(\bm v-\bm{v}_0)$换成 $\bm v$, 问题变为对于任意 $g$ 找到 
    $\bm v \in \mathrm{ker}(\mathrm{inc})$ 满足
    $\mathrm{div}\bm v = g, \bm v\cdot \bm n =
    \tilde{h}$。根据复形的恰当性存在 $\bm u \in H^1(\Omega, \mathbb{R}^3)$
    满足 $\mathrm{def}\bm u =\bm v$ 那么问题可以修改为:
    $$
    \left\{\begin{aligned}
    \mathrm{div(def}\bm{u}) &= \tilde{g}\quad \Omega\\
    \mathrm{def}\bm{u}\cdot \bm{n} &= \tilde{h} \quad\partial\Omega
    \end{aligned}
    \right.
    $$
    $a(u, v) = (\mathrm{def}\bm{u}, \mathrm{def}\bm{v})$在 $H_0^1(\Omega, \mathbb{R}^3)$上的是一个椭圆算子，所以以上问题解存在唯一，问题得证。
