# Agris有限元基函数构造文档
## 一、符号定义
考虑三角形网格 $\mathcal{T}_h$中的一个三角形单元 $T$，三个顶点分别是
$\{\bm{x}_0, \bm{x}_1, \bm{x}_2\}$，记不包含顶点 $i$的边为
$e_i$。定义 $\bm{t}_{ij} = \bm{x}_j - \bm{x}_i$，$\bm{n}_i， \bm{m}_i$ 
分别表示 $e_i$ 
的全局法向和中点（如下图$\bm{n}_i$表示边 $[v, j]$的法向）。

<center>
<img src="tri.png" alt="图片描述" style="width:30%;">
</center>

令$\{\lambda_0, \lambda_1, \lambda_2\}$ 表示 $T$上的重心坐标函数，即 $\lambda_i$
是线性函数满足$\lambda_i(\bm{x}_j) = \delta_{ij}$。记 $\bm{\Lambda}$为一个
$3\times3$的矩阵，$\Lambda_{ij} = \nabla\lambda_i\cdot \bm{n}_j$。对于对称矩阵
$\bm{S}$，定义 $\mathrm{vec}(\bm{S}) = (S_{00}, 2S_{01}, S_{11})^T$ 为其对称部分，
对于矩阵 $\bm M$定义 $\rm{sym}(\bm{M}):= \frac{1}{2}(\bm{M} + \bm{M}^T)$为其对称化。

## 二、Agris 元基函数
Agris 元的形函数空间为 $\mathbb{P}_5(T)$，自由度如下：
$$
\left\{ 
\begin{aligned}
v(\bm{x}_i), \quad i = 0, 1, 2\\
\nabla v(\bm{x}_i), \quad i = 0, 1, 2\\
\nabla^2 v(\bm{x}_i), \quad i = 0, 1, 2\\
\frac{\partial v}{
\partial{\bm{n}_i}} (\bm{m}_i), \quad i = 0, 1, 2
\end{aligned}
\right.
$$
现在构造其对应的基函数，因为三条边，三个顶点上的基函数形式相似，因此我们仅考虑边第
$v$ 条边和第$v$个点上的基函数，$0\leq i < j \leq 2$ 与 $v$不相等。定义边$e_v =
[i, j]$上的函数：
$$
\phi^{e_v} := \frac{16}{\Lambda_{vv}}\lambda_{i}^2\lambda_j^2\lambda_v
$$
顶点 $v$ 上的 6 个间接函数：
$$
\begin{aligned}
\tilde{\phi}_0^v &= \frac{1}{2} \lambda_v^3 \lambda_i^2 - \left( \frac{3}{32} \Lambda_{vj} + \frac{1}{16} \Lambda_{ij} \right) \phi_{e_j} \\
\tilde{\phi}_1^v &= \lambda_v^3 \lambda_i \lambda_j - \frac{1}{16} \Lambda_{jj} \phi_{e_j} - \frac{1}{16} \Lambda_{ii} \phi_{e_i}  \\
\tilde{\phi}_2^v &= \frac{1}{2} \lambda_v^3 \lambda_j^2 - \left( \frac{3}{32} \Lambda_{vi} + \frac{1}{16} \Lambda_{ji} \right) \phi_{e_i} \\
\tilde{\phi}_3^v &= \lambda_v^4 \lambda_i + 8 \tilde{\phi}_0^v + 4 \tilde{\phi}_1^v - \frac{1}{16} \Lambda_{ii} \phi_{e_i} - \left( \frac{1}{4} \Lambda_{vj} + \frac{1}{16} \Lambda_{ij} \right) \phi_{e_j}\\
\tilde{\phi}_4^v &= \lambda_v^4 \lambda_j + 8 \tilde{\phi}_2^v + 4 \tilde{\phi}_1^v - \frac{1}{16} \Lambda_{jj} \phi_{e_j} - \left( \frac{1}{4} \Lambda_{vi} + \frac{1}{16} \Lambda_{ji} \right) \phi_{e_i} \\
\tilde{\phi}_5^v &= \lambda_v^5 - 20 \tilde{\phi}_0^v - 20 \tilde{\phi}_1^v - 20 \tilde{\phi}_2^v + 5 \tilde{\phi}_3^v + 5 \tilde{\phi}_4^v - \frac{5}{16} \Lambda_{vi} \phi_{e_i} - \frac{5}{16} \Lambda_{vj} \phi_{e_j}
\end{aligned}
$$
定义系数矩阵：
$$
\bm{C}^{v,2} = \begin{pmatrix}
\mathrm{vec}(\bm{t}_{vi}^T\bm{t}_{vi}) &
\mathrm{vec}(\mathrm{sym}(\bm{t}_{vi}^T\bm{t}_{vj})) &
\mathrm{vec}(\bm{t}_{vj}^T\bm{t}_{vj})
\end{pmatrix}_{3\times 3}^T, 
\quad 
\bm{C}^{v,1} = 
\begin{pmatrix}
\bm{t}_{vi} & \bm{t}_{vj}
\end{pmatrix}_{2\times 2}^T
$$
那么顶点上的基函数可以表示为：
$$
(\phi^v_0, \phi^v_1, \phi^v_2)^T = 
(\tilde{\phi}^v_0, \tilde{\phi}^v_1, \tilde{\phi}^v_2)
\bm{C}^{v,2},
\quad
(\phi^v_3, \phi^v_4) = 
(\tilde{\phi}^v_3, \tilde{\phi}^v_4)\bm{C}^{v,1},
\quad
\phi^v_5 = \tilde{\phi}^v_5
$$
可以验证 $\phi_0^v, \phi_1^v, \phi_2^v, \phi_3^v, \phi_4^v,\phi_5^v$ 是顶点自由度的对偶基函数。
因此 Agris 元基函数为：
$$
\{\phi^0_i\}_{0\leq i \leq 5}, \cup \{\phi^1_i\}_{0\leq i \leq 5}, \cup
\{\phi^2_i\}_{0\leq i \leq 5}, \cup \{\phi^{[1, 2]}, \phi^{[0, 2]}, \phi^{[0,
1]}\}
$$
