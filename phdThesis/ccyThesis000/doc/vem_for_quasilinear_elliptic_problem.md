# 二阶半线性椭圆方程的后验误差估计

## 问题描述
$$
\begin{aligned}
-\nabla\cdot(\mu(\boldsymbol{x}), |\nabla u|)\nabla u) &= 
f \quad \text{in } \Omega, \\
u &= 0 \quad \text{on } \partial\Omega,
\end{aligned}
$$

问题假设:
- $\mu \in C^0(\bar\Omega\times[0, \infty))$
- 存在常数 $0 < m_\mu \le M_\mu < \infty$ 使得
  $$
  m_\mu(t-s) \le \mu(\boldsymbol{x}, t)t-\mu(\boldsymbol{x}, s)t \le
  M_\mu(t-s) \quad \forall \boldsymbol{x} \in \bar\Omega,\ \  \forall t \ge s
  \ge 0.
  $$
> 第二条假设可以证明存在 $C_1, C_2$ 使得 $\forall v, w \in \mathbb{R}^2$, 
> $x \in \bar\Omega$, 
> $$
> |\mu(\boldsymbol{x}, |v|)v - \mu(\boldsymbol{x}, |w|)w| \le C_1|v-w|,
> $$
> $$
> C_2|v-w|^2 \le (\mu(\boldsymbol{x}, |v|)v - \mu(\boldsymbol{x},
> |w|)w)\cdot(v-w).
> $$
> 另外令 $s = 0$ 可以得到
> $$
> m_\mu \le \mu(\boldsymbol{x}, t) \le M_\mu.
> $$
> 显然这个性质会与问题适定性有关。
对应的变分问题为：找 $u \in H_0^1(\Omega)$ 使得
$$
a(u; u, v) = f(v) \quad \forall v \in H_0^1(\Omega),
$$
其中 $a(u; v, w) = (\mu(\boldsymbol{x}, |\nabla u|)\nabla u, \nabla v)$，
$f(v) = (f, v)$。

## 网格假设
对于多边形网格 $\mathcal{T}_h$，假设存在常数 $\rho > 0$ 使得
- 每个单元 $K \in \mathcal{T}_h$ 都关于一个半径为 $\rho h_K$ 的球是 star shaped 的。
- $h_e \ge \rho h_K$, 对于任一单元 $K$ 的边 $e$ 都成立。
 
在上面的网格假设下：对于任意 $w \in H^s(K)$, $s \ge 0$, 都有如下估计成立：
$$
\|w-Q_k^K w\|_{0, K} + h_K|w-Q_k^K w|_{1, K} \le C h_K^l|w|_{l, K}.
$$
其中 $Q_k^K : L^2(K) \to \mathbb{P}_k(K)$ 是一个 $L^2$ 正交投影，
$l = \min(k, s)$。

## 虚单元方法
使用常见的虚单元空间，定义单元双线性型 $a_K(u, v)$:
$$
a_h(z_h; u_h, v_h) = \sum_{K \in \mathcal{T}_h} a_h^K(z_h; u_h, v_h),
$$
$$
a_h^K(z_h;v_h,w_h):=(\mu(\boldsymbol{x},|\Pi_1^Kz_h|)\Pi_1^Kv_h,\Pi_1^Kw_h)_K+S^K(z_h;(I-\Pi_0^K)v_h,(I-\Pi_0^K)w_h)
$$
其中 $\Pi_0^K v_h$ 是 $v_h$ 到 $\mathbb{P}_k(K)$ 的 $L^2$ 投影，$\Pi_1^K v_h$ 是
$\nabla v_h$ 到 $(\mathbb{P}_{k-1}(K))^2$ 的 $L^2$ 投影。
稳定化项有两种选取方式：
$$
S^K(z_h;v_h,w_h):=M_\mu m_\mu\sum_{i=1}^{N}dof_i^K(v_h)dof_i^K(w_h)
$$
$$
S^K(z_h;v_h,w_h):=\mu_K(\boldsymbol{x}, |\Pi_1^{K, 0}z_h|)\sum_{i=1}^{N}dof_i^K(v_h)dof_i^K(w_h)
$$
这两个稳定化项都有:
$$
\beta_*\int\limits_K\nabla v_h\cdot\nabla v_h \mathrm{d}\boldsymbol{x}\leq S^K(z_h;v_h,v_h)\leq\beta^*\int\limits_K\nabla v_h\cdot\nabla v_h \mathrm{d}\boldsymbol{x}.
$$
所以离散双线性型有如下性质：
$$
\alpha_*a^K(z_h;v_h,v_h)\leq a_h^K(z_h;v_h,v_h)\leq\alpha^*a^K(z_h;v_h,v_h).
$$

离散问题为：找 $u_h \in V_h$ 使得
$$
a_h(u_h; u_h, v_h) = f_h(v_h) \quad \forall v_h \in V_h.
$$
其中 $f_h(v_h) = \sum_{K \in \mathcal{T}_h} (Q_kf, v_h)_K$。


## 问题
1. Below (11), why $C_4 := \min((C_*M_\mu)^{1/2}, (C_*m_\mu)^{1/2})$? 



























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



