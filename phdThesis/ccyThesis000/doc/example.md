
旋度算子定义如下
$$
\mathrm{rot}\mathbf{u} = \frac{\partial u_1}{\partial x} - \frac{\partial u_0}{\partial y},  
\quad \mathbf{rot}u = (\frac{\partial u}{\partial y}, -\frac{\partial u}{\partial x})
$$

$\Omega^{+} = (-1, 0)\times(0, 1)$, $\Omega^{-} = (0, 1)\times(-1, 0)$, 
$\Omega = (-1, 1)\times(-1, 1)$, $\Gamma = (-1, 1)\times\{0\}$

令
$$
\mathbf{u} = \left\{
\begin{aligned}
\begin{pmatrix}
2\cos(ax)\cos(cy)\\
2\sin(ax)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{+},\\
\begin{pmatrix}
\cos(bx)\cos(cy)\\
\sin(bx)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$
$$
\alpha = \left\{
\begin{aligned}
1,\quad  &\text{in } \Omega^{+},\\
1,\quad  &\text{in } \Omega^{-}.
\end{aligned}
\right., \quad\quad
\beta = \left\{
\begin{aligned}
1,\quad &\text{in } \Omega^{+},\\
\beta^-,\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$
那么:
$$
[\![\mathbf{u}\cdot \mathbf{t}]\!]|_\Gamma = 0
$$

$$
\mathrm{rot}\mathbf{u} = \left\{
\begin{aligned}
(2a+2c)\cos(ax)\sin(cy),\quad\quad &\text{in } \Omega^{+},\\
(b+c)\cos(bx)\sin(cy),\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$
当 $c = b-2a$ 时
$$
[\![\alpha\mathrm{rot}\mathbf{u}]\!]|_\Gamma = 0 
$$

<!--  
$$
\mathbf{rot}\alpha\mathrm{rot}\mathbf{u} = \left\{
\begin{aligned}
\begin{pmatrix}
\alpha^+(2ac+2c^2)\cos(ax)\cos(cy)\\
\alpha^+(2a^2+2ac)\sin(ax)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{+},\\
\begin{pmatrix}
\alpha^-(bc+c^2)\cos(bx)\cos(cy)\\
\alpha^-(b^2+bc)\sin(bx)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$-->
取 $a = \sqrt{\frac{\beta^+}{2}}$, $b = \sqrt{\frac{\beta^-}{2}}$,
对于不同 $\beta^-$ 相对误差结果如下:
![相对误差](/home/cbtxs/Figure_2.png)

| Dof   | 32         | 64         | 128        | 256        | 512        | 1024       |
|-------|------------|------------|------------|------------|------------|------------|
| $\frac{\| u - u_h\|_{0}}{h\|u\|_{0}}$         | 1.0135e+00 | 1.3776e+00 | 2.1270e+00 | 3.1661e+00 | 4.8213e+00 | 6.9390e+00 |
| Order | --         | -0.44      | -0.63      | -0.57      | -0.61      | -0.53      |
| $\frac{\|\mathrm{rot} u - \mathrm{rot} u_h\|_0}{h\|u\|_{0}}$ | 3.5999e+00 | 8.0405e+00 | 1.7845e+01 | 3.8505e+01 | 8.2377e+01 | 1.6885e+02 |
| Order | --         | -1.16      | -1.15      | -1.11      | -1.1       | -1.04      |

## 没有间断系数的情况

$\Omega^{+} = (-1, 0)\times(0, 1)$, $\Omega^{-} = (0, 1)\times(-1, 0)$, 
$\Omega = (-1, 1)\times(-1, 1)$, $\Gamma = (-1, 1)\times\{0\}$

令
$$
\mathbf{u} = \left\{
\begin{aligned}
\begin{pmatrix}
2\cos(ax)\cos(cy)\\
2\sin(ax)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{+},\\
\begin{pmatrix}
\cos(bx)\cos(cy)\\
\sin(bx)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$
$$
\alpha = 1, \beta = \beta
$$
那么:
$$
[\![\mathbf{u}\cdot \mathbf{t}]\!]|_\Gamma = 0
$$

$$
\mathrm{rot}\mathbf{u} = \left\{
\begin{aligned}
(2a+2c)\cos(ax)\sin(cy),\quad\quad &\text{in } \Omega^{+},\\
(b+c)\cos(bx)\sin(cy),\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$
当 $c = b-2a$ 时
$$
[\![\alpha\mathrm{rot}\mathbf{u}]\!]|_\Gamma = 0 
$$

<!--  
$$
\mathbf{rot}\alpha\mathrm{rot}\mathbf{u} = \left\{
\begin{aligned}
\begin{pmatrix}
\alpha^+(2ac+2c^2)\cos(ax)\cos(cy)\\
\alpha^+(2a^2+2ac)\sin(ax)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{+},\\
\begin{pmatrix}
\alpha^-(bc+c^2)\cos(bx)\cos(cy)\\
\alpha^-(b^2+bc)\sin(bx)\sin(cy)
\end{pmatrix},\quad\quad &\text{in } \Omega^{-}.
\end{aligned}
\right.
$$-->
取 $a = \sqrt{\frac{\beta^+}{2}}$, $b = \sqrt{\frac{\beta^-}{2}}$,
对于不同 $\beta^-$ 相对误差结果如下:
![相对误差](/home/cbtxs/Figure_2.png)











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
</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>

