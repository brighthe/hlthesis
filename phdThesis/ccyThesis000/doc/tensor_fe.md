# 张量有限元

## 一、Hu-Zhang 元

胡张元是一种 $H(\mathbf{div})$ 协调的二阶对称张量有限元。

- **格林公式**
$$
(\mathbf{div}\boldsymbol{\sigma}, \boldsymbol{v})_{\Omega} =
-(\boldsymbol{\sigma}, \varepsilon(\boldsymbol{v}))_{\Omega} + \langle
\boldsymbol{\sigma}\cdot\boldsymbol{n}, \boldsymbol{v} \rangle_{\partial\Omega}
$$
将 $\Omega$ 剖分为网格 $\mathcal{T}_h$，等式左边等于:
$$
\begin{aligned}
\sum_{K\in\mathcal{T}_h}(\mathbf{div}\boldsymbol{\sigma}, \boldsymbol{v})_{K} 
& = 
\sum_{K\in\mathcal{T}_h}-(\boldsymbol{\sigma}, \varepsilon(\boldsymbol{v}))_{K} + \sum_{K\in\mathcal{T}_h}\langle
\boldsymbol{\sigma}\cdot\boldsymbol{n}, \boldsymbol{v} \rangle_{\partial K}\\
& = 
\sum_{K\in\mathcal{T}_h}-(\boldsymbol{\sigma}, \varepsilon(\boldsymbol{v}))_{K}
+ \sum_{F\in\mathcal{F}_h^i}\langle
[\![\boldsymbol{\sigma}\cdot \boldsymbol{n}_F]\!], \boldsymbol{v} \rangle_{F}
+ \sum_{F\in\mathcal{F}_h^b}\langle
\boldsymbol{\sigma}\cdot\boldsymbol{n}_F, \boldsymbol{v} \rangle_{F}
\end{aligned}
$$
等式右边等于:
$$
\begin{aligned}
\sum_{K\in\mathcal{T}_h}-(\boldsymbol{\sigma}, \varepsilon(\boldsymbol{v}))_{K}
+ \sum_{F\in\mathcal{F}_h^b}\langle
\boldsymbol{\sigma}\cdot\boldsymbol{n}_F, \boldsymbol{v} \rangle_{F}
\end{aligned}
$$
所以，$\boldsymbol{\sigma} \in H(\mathbf{div})$ 的充要条件是
$[\![\boldsymbol{\sigma}\cdot\boldsymbol{n}_F]\!] = 0, \forall F\in\mathcal{F}_h^i$。

- **Hu-Zhang 元的定义**

1. $K$ : 三角形或四面体
2. $\mathcal{P}$ : $\mathbb{P}_k(K, \mathbb{S})$
3. $\mathcal{N}$ :
$$
\begin{aligned}
 \boldsymbol{\sigma}(\delta) \quad &\forall \delta \in \mathcal{V}(K),\\
 (\boldsymbol{n}_i^T \boldsymbol{\sigma} \boldsymbol{n}_j, q)_F \quad &\forall q
 \in \mathbb{P}_{k+r-d-1}(F), F \in \mathcal{F}^r(K), \\
 &\quad i, j = 1, ..., r, r = 1, ..., d-1,\\
 (\boldsymbol{t}_i^T \boldsymbol{\sigma}\boldsymbol{n}_j, q) \quad &\forall q \in \mathbb{P}_{k+r-d-1}(F), F\in \mathcal{F}^r(K), \\
 &\quad i = 1, ..., d-r, j = 1, ..., r, r = 1, ..., d-1.\\
 (\boldsymbol{t}_i^T \boldsymbol{\sigma} \boldsymbol{t}_j, q)_F \quad &\forall q
 \in \mathbb{P}_{k+r-d-1}(F), F \in \mathcal{F}^r(K), \\
 &\quad i, j = 1, ..., d-r, r = 0, ..., d-1\\
\end{aligned}
$$

$\mathcal{T}_h$ 上的胡张元空间为:
$$
\begin{aligned}
\Sigma_h := \{\boldsymbol{\sigma} \in L^2(\Omega, \mathbb{S}) :
& \boldsymbol{\sigma}|_K \in \mathbb{P}_k(K, \mathbb{S}), \forall K \in
\mathcal{T}_h, \\
& \boldsymbol{\sigma} \textrm{ single vlaue at the first three types Dofs}\}
\end{aligned}
$$

## 二、$H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调有限元
- **格林公式**
$$
\begin{aligned}
(\mathrm{div}\mathbf{div}\boldsymbol{\sigma}, v)_{\Omega} = &
(\boldsymbol{\sigma}, \nabla^2 v)_{\Omega} 
-\sum_{F\in \partial \Omega}\sum_{e\in \partial F} (\boldsymbol{n}_{F, e} \boldsymbol{\sigma} \boldsymbol{n}, v)_{e}\\
& - \sum_{F\in \partial \Omega} [(\boldsymbol{n^T\sigma n}, \partial_{\boldsymbol{n}} v)_{F}
-(\boldsymbol{n}^T \mathbf{div}\boldsymbol{\sigma} +  \mathrm{div}_F(\boldsymbol{\sigma n}), v)_{F}]
\end{aligned}
$$
根据格林公式定义三种迹算子：
$$
\begin{aligned}
[\![\mathrm{tr}^1_F(\boldsymbol{\sigma})]\!] &:= \boldsymbol{n}^T_{\partial T}
\boldsymbol{\sigma} \boldsymbol{n}_{\partial T},\\
[\![\mathrm{tr}^2_F(\boldsymbol{\sigma})]\!] &:= 
(\boldsymbol{n}^T_{\partial T} \mathbf{div}(\boldsymbol{\sigma}) + 
\mathrm{div}_F(\boldsymbol{\sigma} \boldsymbol{n}_{\partial T}))|_F,\\
[\![\mathrm{tr}_e(\boldsymbol{\sigma})]\!] &:= \sum_{F\in\partial T, e\in
\partial F} \boldsymbol{n}_{F, e}^T \boldsymbol{\sigma} \boldsymbol{n}_{\partial
T}.
\end{aligned}
$$
与 Hu-Zhang 元类似，$\boldsymbol{\sigma} \in H(\mathrm{div}\mathbf{div})$ 
的充要条件是下面三种跳量为零：
$$
\begin{aligned}
[\![\mathrm{tr}^1_F(\boldsymbol{\sigma})]\!]|_F
&:= \boldsymbol{n}^T_{\partial T_1} \boldsymbol{\sigma} \boldsymbol{n}_{\partial
T_1} - \boldsymbol{n}^T_{\partial T_2} \boldsymbol{\sigma}
\boldsymbol{n}_{\partial T_2} = 0, \\
[\![\mathrm{tr}^2_F(\boldsymbol{\sigma})]\!]|_F &:= 
(\boldsymbol{n}^T_{\partial T_1} \mathbf{div}(\boldsymbol{\sigma}) 
+ \mathrm{div}_F(\boldsymbol{\sigma} \boldsymbol{n}_{\partial T_1}))|_F
+ (\boldsymbol{n}^T_{\partial T_2} \mathbf{div}(\boldsymbol{\sigma})
+ \mathrm{div}_F(\boldsymbol{\sigma} \boldsymbol{n}_{\partial T_2}))|_F = 0,\\
[\![\mathrm{tr}_e(\boldsymbol{\sigma})]\!]|_e 
&:= \sum_{T\in \omega_e}\sum_{F\in\partial T, e\in \partial F}
\boldsymbol{n}_{F, e}^T \boldsymbol{\sigma} \boldsymbol{n}_{\partial T} = 0.
\end{aligned}
$$
那么构造 $H(\mathrm{div}\mathbf{div}, \mathbb{S})$
协调有限元的关键是构造如下空间：
$$
\mathbb{B} := \{\boldsymbol{\sigma} \in \Pi_k :
\boldsymbol{\sigma} \textrm{ statify the jump conditions}\}
$$
其中：
$$
\Pi_k = \{\boldsymbol{\sigma} \in \mathbb{P}_k(K, \mathbb{S})\ \forall K \in
\mathcal{T}_h\}
$$

### （1）$H(\mathrm{div}\mathbf{div}, \mathbb{S}) \cap H(\mathbf{div}, \mathbb{S})$ 协调有限元
为了满足上面的三种跳量为零，Chen-Huang 提出了一种 
$H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调且 $H(\mathbf{div}, \mathbb{S})$
协调的有限元，其要求对于
$\boldsymbol{\sigma} \in \Pi_k$ 满足：
- 对于任意的 $F \in \Delta_1(\mathcal{T}_h)$ : 
  $\boldsymbol{\sigma n_F}$ 和
  $\boldsymbol{n}_F^T\mathbf{div}(\boldsymbol{\sigma})$ 在 $F$ 上连续，

相比于 $H(\mathbf{div}, \mathbb{S})$ 协调有限元，这种有限元多了一个 
$\boldsymbol{n}_F^T\mathbf{div}(\boldsymbol{\sigma})$ 在 $F$ 上连续的条件。
可以对 Hu-Zhang 元进行改进，使其满足上面的条件。

1. $K$ : $d$ 维单纯形
2. $\mathcal{P}$ : $\mathbb{P}_k(K, \mathbb{S})$
3. $\mathcal{N}$ :
$$
\begin{aligned}
 \boldsymbol{\sigma}(\delta) \quad &\forall \delta \in \mathcal{V}(K),\\
 (\boldsymbol{n}_i^T \boldsymbol{\sigma} \boldsymbol{n}_j, q)_F \quad &\forall q
 \in \mathbb{P}_{k+r-d-1}(F), F \in \mathcal{F}^r(K), \\
 &\quad i, j = 1, ..., r, r = 1, ..., d-1,\\
 (\boldsymbol{t}_i^T \boldsymbol{\sigma}\boldsymbol{n}_j, q) \quad &\forall q \in \mathbb{P}_{k+r-d-1}(F), F\in \mathcal{F}^r(K), \\
 &\quad i = 1, ..., d-r, j = 1, ..., r, r = 1, ..., d-1.\\
 (\boldsymbol{n}^T \mathbf{div}\boldsymbol{\sigma}, q) \quad &\forall q \in
 \mathbb{P}_{k-1}(F), F\in \mathcal{F}^1(K).\\
 (\mathrm{div}\mathbf{div}\boldsymbol{\sigma}, q) \quad &\forall q \in \mathbb{P}_{k-2}(K)/\mathbb{P}_1(K).\\
 (\mathbf{div}\boldsymbol{\sigma}, \boldsymbol{q}) \quad &\forall \boldsymbol{q} \in
 (\mathbb{P}_{k-3}(K, \mathbb{K})/\mathbb{P}_0(K, \mathbb{K}))\boldsymbol{x}.\\
 (\boldsymbol{\sigma}, \boldsymbol{q}) \quad &\forall \boldsymbol{q} \in \mathrm{ker}(\cdot \boldsymbol{x})\cap \mathbb{P}_{k-2}(K, \mathbb{S}).\\
\end{aligned}
$$

$\mathcal{T}_h$ 上的 $H(\mathrm{div}\mathbf{div}, \mathbb{S}) \cap H(\mathbf{div}, \mathbb{S})$ 协调有限元空间为:
$$
\begin{aligned}
\Sigma_h^{\mathrm{div}\mathbf{div}, 0} := \{\boldsymbol{\sigma} \in L^2(\Omega, \mathbb{S}) :
& \boldsymbol{\sigma}|_K \in \mathbb{P}_k(K, \mathbb{S}), \forall K \in
\mathcal{T}_h, \\
& \boldsymbol{\sigma} \textrm{ single vlaue at the first four types Dofs}\}
\end{aligned}
$$

### （2）顶点连续的 $H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调有限元
前面对 $\boldsymbol{\sigma n}$ 的连续性要求可以降低为对 $\boldsymbol{\sigma n}$
的法向的连续性
为了满足上面的三种跳量为零，Chen-Huang 提出了一种仅有
$H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调性的有限元，其要求对于
$\boldsymbol{\sigma} \in \Pi_k$ 满足：
- 对任意的 $e \in \Delta (\mathcal{T}_h)$: 
    $\boldsymbol{n}_{e, i}^T \boldsymbol{\sigma} \boldsymbol{n}_{e, j}$ 
    在 $e$ 上连续，
- 对于任意的 $F \in \Delta_1(\mathcal{T}_h)$: $\boldsymbol{n}^T
    \mathbf{div}(\boldsymbol{\sigma})+\mathrm{div}_F(\boldsymbol{\sigma n})$ 在
    $F$ 上连续。

1. $K$ : $d$ 维单纯形
2. $\mathcal{P}$ : $\mathbb{P}_k(K, \mathbb{S})$
3. $\mathcal{N}$ :
$$
\begin{aligned}
 \boldsymbol{\sigma}(\delta) \quad &\forall \delta \in \mathcal{V}(K),\\
 (\boldsymbol{n}_i^T \boldsymbol{\sigma} \boldsymbol{n}_j, q)_F \quad &\forall q
 \in \mathbb{P}_{k+r-d-1}(F), F \in \mathcal{F}^r(K), \\
 &\quad i, j = 1, ..., r, r = 1, ..., d-1,\\
 (\boldsymbol{n}^T \mathbf{div}\boldsymbol{\sigma} +
 \mathrm{div}_F(\boldsymbol{\sigma n}) , q) \quad &\forall q \in
 \mathbb{P}_{k-1}(F), F\in \mathcal{F}^1(K).\\
 (\boldsymbol{t}_i^T \boldsymbol{\sigma}\boldsymbol{n}_j, q) \quad &\forall q \in \mathbb{P}_{k+r-d-1}(F), F\in \mathcal{F}^r(K), \\
 &\quad i = 1, ..., d-r, j = 1, ..., r, r = 1, ..., d-1.\\
 (\boldsymbol{\sigma}, \mathbf{def}{q}) \quad &\forall q \in \mathrm{ND}_{k-3}(K).\\ 
 (\boldsymbol{\sigma}, \boldsymbol{q}) \quad &\forall \boldsymbol{q} \in \mathrm{ker}(\cdot \boldsymbol{x})\cap \mathbb{P}_{k-2}(K, \mathbb{S}).\\
\end{aligned}
$$

$\mathcal{T}_h$ 上的 $H(\mathrm{div}\mathbf{div}, \mathbb{S}) \cap H(\mathbf{div}, \mathbb{S})$ 协调有限元空间为:
$$
\begin{aligned}
\Sigma_h^{\mathrm{div}\mathbf{div}, 1} := \{\boldsymbol{\sigma} \in L^2(\Omega, \mathbb{S}) :
& \boldsymbol{\sigma}|_K \in \mathbb{P}_k(K, \mathbb{S}), \forall K \in
\mathcal{T}_h, \\
& \boldsymbol{\sigma} \textrm{ single vlaue at the first three types Dofs}\}
\end{aligned}
$$

### （3）顶点不连续的 $H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调有限元
前面的两种有限元都要求 $\boldsymbol{\sigma}$
在任意维的子单形的法平面上连续，其原因是为了让
$\boldsymbol{n}^T_i\boldsymbol{\sigma}) \boldsymbol{n}_j = 0$ 在边上连续, 
从而满足 $[\![\mathrm{tr}_e(\boldsymbol{\sigma})]\!] = 0$。

## 四、$\mathbf{div}$ 有限元复形

## 五、$\mathrm{div}\mathbf{div}$ 有限元复形









