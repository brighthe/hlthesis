# 有限元方法的笔记

## 勒贝格积分
**勒贝格控制收敛定理**：设 $f_n \in L^1(\Omega)$, $f_n \to f$
几乎处处收敛，且存在 $g \in L^1(\Omega)$, 使得 $|f_n| \le g$
几乎处处成立，则有：
$$
\lim_{n\to\infty} \int_{\Omega} f_n \, dx = \int_{\Omega} f \, dx.
$$

## 多项式理论

### 平均 Taylor 多项式

**cut-off 函数**：任取 $B\subset\Omega$ 定义 $B$ 上的函数:
$$
\psi(x) = \begin{cases}
e^{\frac{-1}{1 - \left(\frac{\|x - x_0\|}{r}\right)^2}}, 
& x \in B, \\
0, & x \notin B
\end{cases}
$$
$c = \int_{\mathbb{R}^d} \psi(x) \, dx$ 是一个常数。cut-off 函数定义为：
$\psi(x) = \frac{1}{c} \psi(x)$.
满足:
1. $supp(\phi) = \bar B$.
2. $\int_{\mathbb{R}^d} \phi(x) \, dx = 1$.

**Taylor 展开**：对于 $f \in C^m(\mathbb{R}^d)$, $x_0 \in \mathbb{R}^d$, $r > 0$, 定义：
$$
T_y^m f(x) = \sum_{|\alpha| \le m} \frac{1}{\alpha!} D^\alpha f(y)(x - y)^{\alpha}
$$
其中 $(x-y)^\alpha = \Pi_{i=0}^{d-1} (x_i - y_i)^{\alpha_i} = \sum_{\beta+\gamma
= \alpha}a(\beta, \gamma)x^\beta y^\gamma$.
其中 $a(\beta, \gamma) = (-1)^{|\gamma|}\frac{(\beta+\gamma)!}{\beta!\gamma!}$.
是一个常数。

**平均算子**：对于 $f \in C^m(\mathbb{R}^d)$, $x_0 \in \mathbb{R}^d$, $r > 0$, 定义：
$$
Q^m f(x) = \int_{B} T_y^m f(x) \phi(y) \, dy.
$$

那么有：
$$
\begin{aligned}
Q^m f(x) & = \sum_{|\alpha| \le m} \int_{B}\frac{1}{\alpha!} D^\alpha f(y)
(y - x)^\alpha \phi(y) \, dy.\\
& = \sum_{|\alpha| \le m}\sum_{\beta+\gamma = \alpha} \frac{1}{\alpha!}a(\beta, \gamma) x^\beta \int_{B} D^\alpha f(y)y^\gamma \phi(y) \, dy.\\
& = \sum_{|\alpha| \le m}\sum_{\beta+\gamma = \alpha} \frac{(-1)^{|\alpha|}}{\alpha!}a(\beta,
\gamma) x^\beta \int_{B} f(y) D^\alpha(y^\gamma \phi(y)) \, dy.\\
\end{aligned}
$$
所以 $Q^m f(x)$ 是一个 $m$ 次多项式。且按照上面的式子可将算子的定义扩展到 $f \in
L^1(\Omega)$.

**平均算子的性质**：
1. $Q^m$ 是一个 $L^1(\Omega)$ 到 $W^{k,\infty}(\Omega)$ 的有界线性算子。
2. $Q^m f(x) = f(x)$, 如果 $f$ 是一个 $m$ 次多项式。
3. $D^{\alpha} Q^m f(x) = Q^{m-|\alpha|} D^{\alpha} f(x), \forall f \in W^{|\alpha|,1}(\Omega)$.

   证明: $T_{d}^{k-l} = \{\alpha\in T_d^k, \alpha - \beta > 0\}\ \forall \beta
   \in T_{d}^l$, 所以当 $f \in C^m(\Omega)$：
   $$
   \begin{aligned}
   D^{\alpha} T_y^m f(x) & = \sum_{|\beta| \le m} \frac{1}{\beta!} D^{\beta} f(y) D^{\alpha}((x-y)^\beta).\\
    & = \sum_{|\beta| \le m, \beta \le \alpha} \frac{1}{(\beta-\alpha)!} 
    D^{\beta} f(y)(x-y)^{\beta-\alpha}.\\
    & = \sum_{|\gamma| \le m-|\alpha|} \frac{1}{\gamma!} D^{\gamma+\alpha} f(y)(x-y)^{\gamma}.\\
    & = T_y^{m-|\alpha|} D^{\alpha} f(x). 
   \end{aligned}
   $$
   所以：
   $$
   \begin{aligned}
    D^{\alpha} Q^m f(x) & = \int_{B} D^{\alpha} T_y^m f(x) \phi(y) \, dy.\\
    & = \int_{B} T_y^{m-|\alpha|} D^{\alpha} f(x) \phi(y) \, dy.\\
    & = Q^{m-|\alpha|} D^{\alpha} f(x).
   \end{aligned}
   $$
   若 $f \in W^{|\alpha|,1}(\Omega)$, 则 $f$ 可以由 $C^\infty(\Omega)$ 中的函数逼近，
   $$
   \mathfrak{ABCDEFGHIJKLMNOPQRSTUVWXYZ}
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
</br>
</br>
</br>
</br>
</br>








