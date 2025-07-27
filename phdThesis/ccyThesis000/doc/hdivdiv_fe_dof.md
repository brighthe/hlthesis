# $H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调有限元的自由度和基函数


## 一、迹的定义
考虑 $\boldsymbol{\sigma} \in \mathbb{P}_k(T, \mathbb{S})$:
$$
\begin{aligned}
  \mathrm{tr}_F^1(\boldsymbol{\sigma}) & := \boldsymbol{n}^T_F \boldsymbol{\sigma}
  \boldsymbol{n}_F, \\
  \mathrm{tr}_F^2(\boldsymbol{\sigma}) & := \boldsymbol{n}_F^T
  \mathrm{div}(\boldsymbol{\sigma}) + \mathrm{div}_F(\boldsymbol{\sigma n}_F),
  \\
  \mathrm{tr}_e(\boldsymbol{\sigma}) & := \sum_{F\in \partial T, e \in \partial F}
  \boldsymbol{n}_{F,e}^T \boldsymbol{\sigma} \boldsymbol{n}_F, \\
\end{aligned}
$$

假设 $\boldsymbol{\sigma} = f \boldsymbol{M}$, 那么:
$$
\begin{aligned}
  \mathrm{tr}_F^1(\boldsymbol{\sigma}) & = f \boldsymbol{n}^T_F \boldsymbol{M}
  \boldsymbol{n}_F\\
  \mathrm{tr}_F^2(\boldsymbol{\sigma}) & = (\nabla f + \Pi_F \nabla f) \boldsymbol{M}\boldsymbol{n}_F\\
  \mathrm{tr}_e(\boldsymbol{\sigma}) & = \sum_{F\in \partial T, e \in \partial F}
  f \boldsymbol{n}_{F,e}^T \boldsymbol{M} \boldsymbol{n}_F
\end{aligned}
$$


## 一、二维情况



