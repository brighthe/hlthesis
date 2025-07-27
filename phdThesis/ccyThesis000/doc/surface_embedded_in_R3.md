# 嵌入在三维空间中的曲面

## 一、经典曲面论
在经典曲面论中，三维空间中的一些曲面 $S$ 可以用参数方程表示为
$$
\boldsymbol{r}(u_0, u_1) = 
x_0(u_0, u_1) \boldsymbol{i} + x_1(u_0, u_1) \boldsymbol{j} + x_2(u_0, u_1) \boldsymbol{k}
$$
其中 $x_i(u_0, u_1)$ 是 $\mathbb{R}^2$ 的一个开集 $D$ 上的连续函数。

根据参数方程我们可以得到曲面上的两个切向量: $\boldsymbol{r}_{u_0},
\boldsymbol{r}_{u_1}$，以及曲面上的法向量 
$\boldsymbol{n} = \boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1}$。映射的 Jacobi 矩阵
为:
$$
J = \frac{\partial \boldsymbol{r}}{\partial \boldsymbol{u}} =
(\frac{\partial x_i}{\partial u_j})
$$
第一基本形式为：
$$
I = \boldsymbol{r}_{u_0}\cdot \boldsymbol{r}_{u_0}(\mathrm{d}u_0)^2 +
\boldsymbol{r}_{u_0}\cdot \boldsymbol{r}_{u_1}\,\mathrm{d}u_0\,\mathrm{d}u_1 +
\boldsymbol{r}_{u_1}\cdot \boldsymbol{r}_{u_1}(\mathrm{d}u_1)^2
$$

## 二、参数曲面的微分几何
重写 $\boldsymbol{r}(u_0, u_1)$ 为 $D$ 到 $S$ 上的映射 $\phi$，现在定义 $D$ 
上的基矢为 $\partial_{u_0}, \partial_{u_1}$，对偶基矢为 $\mathrm{d}u_0, \mathrm{d}u_1$，
那么 $\{\phi_{*}\partial_{u_0}, \phi_{*}\partial_{u_1}\}$，
$\{\phi_{*} \mathrm{d}u_0, \phi_{*} \mathrm{d}u_1\}$ 是 $S$
上的基矢和对偶基矢，为了方便我们仍然记为 $\{\partial_{u_0}, \partial_{u_1}\}$，
$\{\mathrm{d}u_0, \mathrm{d}u_1\}$。其中 $\partial_{u_0}, \partial_{u_1}$ 是 $S$ 上的切向量，
他们分别对应与三维空间中的向量 $\boldsymbol{r}_{u_0}, \boldsymbol{r}_{u_1}$。

记 $R^3$ 中的度量为 $\delta_{ab}$，那么 $S$ 上的度量是 $\delta_{ab}$ 在 $S$ 上的限制。
记 $S$ 上的度量为 $g_{ab}$，那么
$$
g_{ab}(\partial_{u_i})^a(\partial_{u_j})^b = \boldsymbol{r}_{u_i}\cdot \boldsymbol{r}_{u_j} 
$$
所以：
$$
g_{ab} = \boldsymbol{r}_{u_i}\cdot \boldsymbol{r}_{u_j}(\mathrm{d}u_i)_a\otimes
(\mathrm{d}u_j)_b
$$
> 可见第一基本形式就是黎曼度量。

即：$g_{ij} := \boldsymbol{r}_{u_i}\cdot \boldsymbol{r}_{u_j}$,
定义 $g^{ij}$ 为 $g_{ij}$ 的逆矩阵, $g^{ij}g_{jk} = \delta^i_k$，
令 $\{\tilde{\boldsymbol{r}}_{u_i}\}$ 是 $\{\boldsymbol{r}_{u_i}\}$ 的对偶基，
因为
$$
\boldsymbol{v} = (\boldsymbol{v}\cdot \tilde{\boldsymbol{r}}_{u_i})\boldsymbol{r}_{u_i}
$$
所以
$$
(\boldsymbol{r}_{u_i}\cdot \boldsymbol{r}_{u_j})
(\boldsymbol{\tilde{r}}_{u_j}\cdot \tilde{\boldsymbol{r}}_{u_k})
=
\boldsymbol{r}_{u_i}\cdot
(\boldsymbol{r}_{u_j}(\tilde{\boldsymbol{r}}_{u_j}\cdot
\tilde{\boldsymbol{r}}_{u_k}))
= \boldsymbol{r}_{u_i}\cdot \tilde{\boldsymbol{r}}_{u_k} = \delta_{ik}
$$
所以可以写出 $g^{ab}$ 的分量形式 $g^{ij} = \tilde{\boldsymbol{r}}_{u_i}\cdot
\tilde{\boldsymbol{r}}_{u_j}$。
$$
|g^{ij}| = |\tilde{\boldsymbol{r}}_{u_0}\times\tilde{\boldsymbol{r}}_{u_1}|^2 =
|\boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1}|^{-2}
$$
> 证明：对于任意两个向量 $\boldsymbol{a}, \boldsymbol{b}$，
> 令 $\boldsymbol{c} = \frac{\boldsymbol{a}\times \boldsymbol{b}}{|\boldsymbol{a}\times \boldsymbol{b}|}$，
> 那么矩阵 $A = \begin{pmatrix} \boldsymbol{a} \\ \boldsymbol{b} \\ \boldsymbol{c} \end{pmatrix}$ 的行列式为 $(\boldsymbol{a}\times \boldsymbol{b})\cdot \boldsymbol{c = 
> |\boldsymbol{a}\times \boldsymbol{b}|}$。
>  令 $B = \begin{pmatrix} \boldsymbol{a} \\ \boldsymbol{b} \end{pmatrix}$， 
>  那么 $|B^TB| = |A^TA| = |A|^2$, 所以上面第一个等式成立。
>  令 $\tilde{\boldsymbol{a}}, \tilde{\boldsymbol{b}}$ 是 $\boldsymbol{a}, \boldsymbol{b}$ 
>  在所处平面上的对偶向量，令 $C = \begin{pmatrix} \tilde{\boldsymbol{a}} \\ 
>  \tilde{\boldsymbol{b}} \\ \boldsymbol{c} \end{pmatrix}$，那么显然有
>  $|C| = |\tilde{\boldsymbol{a}}\times \tilde{\boldsymbol{b}}|$, 
>  且 $AC = I$，所以 $|\tilde{\boldsymbol{a}}\times \tilde{\boldsymbol{b}}| = |C| = 
>  |A|^{-1} = |\boldsymbol{a}\times \boldsymbol{b}|^{-1}$。
>
>
$D$ 上的体元为:
$$
\tau_{ab} = (\mathrm{d}u_1)_a\wedge (\mathrm{d}u_2)_b = (\mathrm{d}u_1)_a\otimes
(\mathrm{d}u_2)_b - (\mathrm{d}u_2)_a\otimes (\mathrm{d}u_1)_b
$$
那么
$$
\begin{aligned}
\tau^{ab} = g^{ac}g^{bd}\tau_{cd}
& = g^{ac}g^{bd}(\mathrm{d}u_1)_c\otimes
    (\mathrm{d}u_2)_d - g^{ac}g^{bd}(\mathrm{d}u_2)_c\otimes (\mathrm{d}u_1)_d\\
& = g^{i1}g^{j2} (\boldsymbol{r}_{u_i})^a\otimes (\boldsymbol{r}_{u_j})^b
    - g^{i2}g^{j1} (\boldsymbol{r}_{u_i})^a\otimes (\boldsymbol{r}_{u_j})^b\\
& = (g^{i1}g^{j2} - g^{i2}g^{j1}) (\boldsymbol{r}_{u_i})^a\otimes (\boldsymbol{r}_{u_j})^b\\
\end{aligned}
$$
那么
$$
\tau^{ab}\tau_{ab} = (g^{11}g^{22} - g^{12}g^{21}) -
(g^{21}g^{12} - g^{22}g^{11}) = g^{11}g^{22} - g^{12}g^{21} = |g^{ij}|
$$
所以 $S$ 上的体元为 $|\boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1}| \mathrm{d}u_0\wedge \mathrm{d}u_1$。
> **注:** 从另一个角度计算，假设 $\varepsilon_{abc}$ 是 $\mathbb{R}^3$ 中的体元，那么
> $S$ 上的面元为 $\varepsilon_{abc}$ 的诱导面元: $\varepsilon_{abc}n^a$ ，其中 $n^a$ 是 $S$ 在 $\mathbb{R}$ 
> 上的单位法向量，因为 $\varepsilon_{abc}n^a \partial_{u_0}^b \partial_{u_1}^c = 
> (\boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1})\cdot \boldsymbol{n} = 
> |\boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1}|$，所以 $\varepsilon_{abc}n^a =
> |\boldsymbol{r}_{u_0}\times \boldsymbol{r}_{u_1}|\mathrm{d}u_0\wedge \mathrm{d}u_1$。

### 标量场推前的切向梯度
在 $\mathbb{R}^3$ 上定义一组标准正交标架场 $\bm n, \bm t_0, \bm t_1$，那么
$\mathbb{R}^3$ 上的一个函数 $f$ 的梯度可以写作：
$$
\partial_n f \bm n + \partial_{t_0} f \bm t_0 + \partial_{t_1} f \bm t_1
$$
当 $\bm n, \bm t_0, \bm t_1$ 不是标准正交标架，但 $\bm n$ 是单位向量与 $\bm t_0, \bm t_1$ 正交，
那么 $f$ 的梯度可以写作：
$$
\partial_n f \bm n + \partial_{t_0} f \tilde{\bm t}_0 + \partial_{t_1} f \tilde{\bm t}_1
$$
其中 $\tilde{\bm t}_0, \tilde{\bm t}_1$ 是 $\bm t_0, \bm t_1$ 所在平面的对偶基。

现在考虑 $D$ 为一个二维区域，$S$ 为 $\mathbb{R}^3$ 中的一个曲面，$\phi$ 是 $D$ 到 $S$ 的一个微分同胚，
$$
\bm x = \phi(\bm u)
$$
在$D$ 上定义标量场 $\hat{f}$，那么 $f$ 在 $S$ 上的推前为
$$
f(\bm x) = \hat{f}(\phi^{-1}(\bm x))
$$
那么 $f$ 可以看做是一个三变量的函数，但是 $f$ 仅在 $S$ 上有定义，所以 $f$ 只能计算切向梯度，
$\bm x_{u_0}, \bm x_{u_1}$ 是 $S$ 上的切向量，
在 $S$ 上定义标架 $\bm n, \bm x_{u_0}, \bm x_{u_1}$，那么 $f$ 的切向梯度可以写作：
$$
\partial_{u_0} f \tilde{\bm x}_{u_0} + \partial_{u_1} f \tilde{\bm x}_{u_1}
$$
其中 $\tilde{\bm x}_{u_0}, \tilde{\bm x}_{u_1}$ 是 $\bm x_{u_0}, \bm x_{u_1}$ 的
对偶基。其计算方式为：
$$
\tilde{\bm x}_{u_0} = \bm g \bm x_{u_0}, \quad
\tilde{\bm x}_{u_1} = \bm g \bm x_{u_1}
$$
其中 $\bm g$ 是 $S$ 上的度量矩阵:
$$
\bm g = g_{ij}\bm x_{u_i}\otimes \bm x_{u_j}, \quad 
g_{ij} = \bm x_{u_i}\cdot \bm x_{u_j}
$$

## 梯度的推前
令 $\phi : M\to N$ 是一个微分同胚，$f$ 是 $M$ 上的一个标量场，
对于一个标量场 $f$，在 $N$ 上的推前为 $\tilde{f} := \phi_* f$，定义如下度量:
- 在 $M$ 上的度量为 $\bm g_{ab}^M$，与之适配的梯度算子为 $\nabla_a^M$，
- 在 $N$ 上的度量为 $\bm g_{ab}^N$，与之适配的梯度算子为 $\nabla_i^N$，
- 映射 $\phi$ 将 $g_{ab}^M$ 推前到 $N$ 上，得到 $N$ 上另一个度量 $\tilde{g}_{ab}^N$，
    与 $\tilde{g}_{ab}^N$ 适配的梯度算子为 $\tilde{\nabla}_a^N$，

$f$ 的梯度为 $\nabla_a^M f$，$\tilde{f}$ 的梯度为 $\nabla_a^N \tilde{f}$，
那么 $\phi_*\nabla^M_a f$ 满足:
$$
\phi_*\nabla^M_a f(v^b) = \nabla^M_a f(\phi^* v^b) = \phi^*v^b(f) = v^b(\phi_* f) 
= \nabla^N_b \phi_*f(v^b)
$$
所以 $\phi_*\nabla^M_a f = \nabla^N_a \phi_*f$，即标量场梯度的推前为推前的梯度。
证明的关键在于任意梯度 $\nabla_a$ 满足 $\nabla_a f(v^b) = v^b(f)$。
那么对于一般的张量场 $T_{ab}$ 是否有这样的性质呢？答案是否定的，以度量张量为例，
$$
\nabla_{a}^M g_{bc}^M = 0
$$
若 $g_{ab}^M$ 梯度的推前等于推前的梯度，那么
$$
0 = \phi_*(\nabla_{a}^M g_{bc}^M) = \nabla_{a}^N \phi_*(g_{bc}^N) 
= \nabla_{a}^N \tilde{g}_{bc}^N
$$
那么 $\tilde{g}_{bc}^N$ 就是 $\nabla_a^N$ 适配的度量，这显然与 $\nabla_a^N$ 
的任意性矛盾。

但是如果我们固定 $N$ 上的度量为 $\tilde{g}_{ab}^N$，梯度算子为 $\tilde{\nabla}_a^N$，
那么上述结论就会成立。首先在 $M$ 定义坐标系 $\{x_i\}$，坐标基矢为 
$\{(\partial x_i)^a\}$， 
对偶基矢为 $\{(\mathrm{d}x_i)_a\}$，那么 $N$ 上就可以被推前为坐标系 $\{y_i\}$，
坐标基矢为 $\{(\partial y_i)^a\}$，对偶基矢为 $\{(\mathrm{d}y_i)_a\}$，

那么对于任意的张量场 $T_{b_1\cdots b_k}^{a_1\cdots a_l}$:
$$
\nabla_a^M T_{b_1\cdots b_k}^{a_1\cdots a_l} = 
(\mathrm{d} T_{\nu_1\cdots \nu_k}^{\mu_1\cdots \mu_l})_a
(\mathrm{d}x_{\mu_1})^{a_1}\cdots (\mathrm{d}x_{\mu_l})^{a_l}
(\partial x_{\nu_1})^{b_1}\cdots (\partial x_{\nu_k})^{b_k}
$$
$$
\tilde{\nabla}_a^N \phi_* T_{b_1\cdots b_k}^{a_1\cdots a_l} = 
(\mathrm{d} T_{\nu_1\cdots \nu_k}^{\mu_1\cdots \mu_l})_a
(\mathrm{d}y_{\mu_1})^{a_1}\cdots (\mathrm{d}y_{\mu_l})^{a_l}
(\partial y_{\nu_1})^{b_1}\cdots (\partial y_{\nu_k})^{b_k}
$$
注意 $T_{\nu_1\cdots \nu_k}^{\mu_1\cdots \mu_l}$ 是一个标量场，满足推前的梯度等于梯度的推前，
所以 $\phi_*\nabla^M_a T = \tilde{\nabla}_a^N \phi_* T$。

推前的梯度等于梯度的推前等价与推前算子与梯度算子的交换，
现在我们讨论那些特殊的微分算子能与推前拉回算子交换，
这个问题其实等价于讨论什么样的微分算子 $\mathrm{d}$ 满足
$$
\mathrm{d} \phi^* - \phi^* \mathrm{d} = 0
$$
假设现在 $N$ 上的坐标基矢就是 $M$ 上的坐标基矢的推前，即 $(\partial y_i)^a$
记 $\nabla_a^M$ 的克氏符为 $\Gamma_{ab}^c$，那么对于任意的张量场 $T_{b_1\cdots b_k}^{a_1\cdots a_l}$，
$$
(\nabla_a^M - \tilde{\nabla}_a^M) T_{b_1\cdots b_k}^{a_1\cdots a_l} = 
\sum_{i}\Gamma_{ad}^{a_i} T_{b_1\cdots b_k}^{a_1\cdots a_{i-1}da_{i+1}\cdots a_l} -
\sum_{j}\Gamma_{ab_j}^{c} T_{b_1\cdots b_{j-1}cb_{j+1}\cdots b_k}^{a_1\cdots a_l}
$$
现在仅考虑 $T$ 是 (0, k) 张量场的情况，那么
$$
(\nabla_a^M - \tilde{\nabla}_a^M) T_{b_1\cdots b_k}=
-\sum_{j}\Gamma_{ab_j}^{c} T_{b_1\cdots b_{j-1}cb_{j+1}\cdots b_k}
$$
现在我们只考虑梯度的反对称部分:
$$
(\nabla_a^M - \tilde{\nabla}_a^M)_{[a} T_{b_1\cdots b_k]}=
-\sum_{j}\Gamma_{[ab_j}^{c} T_{b_1\cdots b_{j-1}|c|b_{j+1}\cdots b_k]}
= -\sum_{j}\Gamma_{[(ab_j)}^{c} T_{b_1\cdots b_{j-1}c b_{j+1}\cdots b_k]}
= 0
$$
其中第二个等号是因为 $\Gamma_{ab}^c$ 是对称的。
所以任意一个梯度算子的反对称部分相等，可以与推前算子交换。



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














