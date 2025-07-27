# VEM for Heat Equatuion with Interface

## 1. Model
Consider a heat euqation: 
$$
\left\{
\begin{aligned}
\frac{\partial u}{\partial t}-\Delta u & = f\quad & (\boldsymbol{x}, t) \in \Omega\times(0, 1)\\
\quad u & = u(x, 0) & t = 0\\
u & = g & \boldsymbol{x} \in \partial\Omega
\end{aligned}
\right.
$$
对时间项进行向后欧拉离散，方程变为:
$$
u^{n+1} - \Delta t\Delta u^{n+1} = \Delta t f^{n+1} + u^n
$$
假设 $u^n$ 已知. 对应的变分问题为: find $u^{n+1} \in H^1(\Omega)$ statify 
$\forall v\in H^1(\Omega)$:
$$
\begin{aligned}
  (u^{n+1}, v) + \Delta t(\nabla u^{n+1}, \nabla v) = 
  \Delta t(f^{n+1}, v) + (u^n, v)
\end{aligned}
$$
## 2. VEM
将 $\Omega$ 离散为多边形网格 $\mathcal T_h$, 定义一个多边形单元 $K$ 上的虚单元空间：
$$
V_k = \{\}
$$
其中 ...

Define a bilinear form on $V_k(K)$:
$$
a_h(u_h, v_h) = (\Pi^0u, \Pi^0 v) + 
\Delta t(\nabla \Pi^1u, \nabla \Pi^1 v) + S((I-\Pi^1)u_h, (I-\Pi^1)v_h)
$$
where :
$$
S(u, v) = \mathrm{dof}(u)\mathrm{dof}(v)
$$
A linear form:
$$
b_h^{n+1}(v_h) = (f^{n+1}, \Pi^0v_h) + (\Pi^0 u_h^{n}, \Pi^0v_h)
$$
Find $u \in V_k, \forall v \in V_{k, 0}$ statify :
$$
a_h(u_h^{n+1}, v_h) = b_h^{n+1}(v_h)
$$




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



