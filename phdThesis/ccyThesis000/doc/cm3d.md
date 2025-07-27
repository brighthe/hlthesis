# 三维 $C^m$ 光滑元

1. ${\rm sym}(\bm\Lambda^{\bm\alpha})(\bm{t}^{\bm\beta}_{\delta_i})=
   \left\{\begin{aligned}
   & \frac{(-1)^{\alpha_i}\bm\beta!\alpha_i!}{(\bm\alpha_{\to
   i}-\bm\beta)!|\bm\alpha|!} \quad & \bm{\beta-\alpha_{\to i}}\geq0\\ &
   0 \quad & {\rm other} \end{aligned} \right.$

**证明**: 以 $i=0$ 为例， 因为.
所以当 $\alpha_1>\beta_0$ 或 $\alpha_2>\beta_1$ 时值为 $0$。
$\bm{\Lambda}_i(\bm t_{0, j})= \delta_{ij}, i, j = 1, 2, 3$, 
$\bm{\Lambda}_0(\bm t_{0, i})= -1, i = 1, 2, 3$, 
所以:

$$
\begin{aligned}
{\rm sym}(
\bm\Lambda_0^{\alpha_0}\otimes
\bm\Lambda_1^{\alpha_1}\otimes
\bm\Lambda_2^{\alpha_2}\otimes
\bm\Lambda_3^{\alpha_3})
((\bm t_{0, 1})^{\beta_0}\otimes
(\bm t_{0, 2})^{\beta_1}\otimes
(\bm t_{0, 3})^{\beta_2})
= \frac{(-1)^{\alpha_0}}{|\bm \alpha|!}
\begin{pmatrix}
\alpha_1\\\beta_0
\end{pmatrix}
\begin{pmatrix}
\alpha_2\\\beta_1
\end{pmatrix}
\begin{pmatrix}
\alpha_0\\\alpha_0
\end{pmatrix}=
\frac{(-1)^{\alpha_0}\bm\beta!\alpha_0!}{(\bm\alpha_{\to 0}-\bm\beta)!|\bm\alpha|!}
\end{aligned}
$$

3. $\phi^{(k-r)e_i}(\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{v_{i}}^{\bm\gamma_{\to i}}} ) = \left\{
\begin{aligned}
&\frac{(-1)^{\alpha_i-k+r}k!\bm{\gamma}_{\to i}!}{(k-r)!\bm{\alpha}_{\to i}!(\bm{\gamma}_{\to i}-\bm{\alpha_{\to i}})!}
\quad
&\bm{\gamma}_{\to i}-\bm{\alpha}_{\to i}\geq0\\
&0 & {\rm other}
\end{aligned}
\right.$

证明：
$$
\begin{aligned}
\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{v_{i}}^{\bm\gamma_{\to i}}} 
= & \nabla^rB^{\alpha}(\bm N^{\bm\gamma_{\to{i}}}_{v_i})\\
= & \sum_{\bm{\beta} \in\mathbb{T}_r^d, \bm{\alpha-\beta}\geq0}
\frac{r!k!}{(k-r)!\bm\beta!}
{\rm sym}(\bm\Lambda^{\bm{\beta}})(\bm N^{\bm\gamma_{\to i}}_{v_i})
B^{\bm{\alpha-\beta}}\\
= & 
\sum_{\bm\beta\in\mathbb{T}_r^d,\ \bm\alpha-\bm\beta\geq0,\ \bm\beta_{\to i}-\bm\gamma_{\to i}\geq0}
\frac{(-1)^{\beta_i}\bm\gamma_{\to i}!k!}
{(k-r)!\bm\beta_{\to i}!(\bm\beta_{\to i}-\bm\gamma_{\to i})!} B^{\bm\alpha-\bm\beta}
\end{aligned}
$$
显然当 $\bm\alpha-\bm\beta = (k-r)\bm e_i$ 且 $\bm\beta_{\to i}-\bm\gamma_{\to i}\geq0$ 
时:
$$
\phi^{(k-r)\bm e_i}
(\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{v_{i}}^{\bm\gamma_{\to i}}})
= 
\frac{(-1)^{\beta_i}\bm\gamma_{\to i}!k!}
{(k-r)!\bm\beta_{\to i}!(\bm\beta_{\to i}-\bm\gamma_{\to i})!} B^{\bm\alpha-\bm\beta}
$$
否则：
$$
\phi^{(k-r)\bm e_i}
(\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{v_{i}}^{\bm\gamma_{\to i}}})
 = 0
$$
所以上式得证。(注意：当 $\bm\alpha-\bm\beta = (k-r)\bm e_i$ 时，$\bm\alpha_{\to i} = \bm\beta_{\to i}$)

4. 
${\rm sym}(\bm \Lambda^{\bm\beta})(\bm N_{e_{ij}}^{\bm \gamma}) = 
\left\{
\begin{aligned}
\frac{\bm\beta!}{|\bm \beta|!}\frac{\bm\gamma!}{(\bm\gamma-\bm\beta_{\to i, j})!
\bm\beta_{\to i, j}!}
\sum_{\bm\sigma\in\mathbb T_{\beta_0}^1, \bm\gamma-\bm\sigma-\bm\beta_{\to{0,1}}\geq0}&\frac{(\bm\gamma-\bm\beta_{\to{i,j}})!}{\bm{\sigma}!(\bm\gamma-\bm\beta_{\to{i,j}}-\bm\sigma)!} \bm E_{ij}^{\bm \sigma\gets(\bm \gamma-\bm \sigma-\bm \beta_{\to i, j })}
& \quad \bm \gamma-\bm \beta_{\to i, j }\geq0\\
&0
& {\rm other}
\end{aligned} 
\right.
$

证明：$\nabla\lambda_k \cdot\bm N_{e_{ij}} = (1, 0),\ \nabla\lambda_l \cdot\bm N_{e_{ij}} = (0, 1),$ 以边 $\{0, 1\}$ 为例，
$$
\begin{aligned}
{\rm sym}(
\bm\Lambda_0^{\beta_0}\otimes
\bm\Lambda_1^{\beta_1}\otimes
\bm\Lambda_2^{\beta_2}\otimes
\bm\Lambda_3^{\beta_3})
& ((\bm N_{e_{01}, 0})^{\gamma_0}\otimes
(\bm N_{e_{01}1})^{\gamma_1}
= \\
& \frac{\bm\beta!}{|\bm \beta|!}\frac{\bm\gamma!}{(\bm\gamma-\bm\beta_{\to i, j})!
\bm\beta_{\to i, j}!}
\sum_{\bm\sigma\in\mathbb T_{\beta_0}^1, \bm\gamma-\bm\sigma-\bm\beta_{\to{0,1}}\geq0}\frac{(\bm\gamma-\bm\beta_{\to{i,j}})!}{\bm{\sigma}!(\bm\gamma-\bm\beta_{\to{i,j}}-\bm\sigma)!} \bm E_{ij}^{\bm \sigma\gets(\bm \gamma-\bm \sigma-\bm \beta_{\to i, j })}
\end{aligned}$$ 

5. 
$
\begin{aligned}
\phi^{\bm \gamma_{\{i,j\}}}(\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{e_{ij}}^{\bm\gamma_{\to i, j}}} ) = \left\{
\begin{aligned}
\frac{k!}{(k-r)!} \sum_{\bm\sigma\in\mathbb T_{\beta_0}^1,
\bm\gamma_{\to i,j }-\bm\sigma-\bm\beta_{\to{i,j}}\geq0}
& \frac{(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}})!}{\bm{\sigma}!(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}}-\bm\sigma)!}
\bm E_{ij}^{\bm \sigma\gets(\bm \gamma_{\to i, j}-\bm \sigma-\bm \beta_{\to i, j })} 
& \bm\gamma_{\to i,j}-\bm{\beta}_{\to i, j}\geq 0\\
& 0 & {\rm other}
\end{aligned}
\right.
\end{aligned}
$

证明:
$$
\begin{aligned}
\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{e_{i, j}}^{\bm\gamma_{\to i, j}}} 
= & \nabla^rB^{\alpha}(\bm N^{\bm\gamma_{\to{i, j}}}_{e_{i, j}})\\
= & \sum_{\bm{\beta} \in\mathbb{T}_r^d, \bm{\alpha-\beta}\geq0}
\frac{r!k!}{(k-r)!\bm\beta!}
{\rm sym}(\bm\Lambda^{\bm{\beta}})(\bm N^{\bm\gamma_{\to i, j}}_{e_{i, j}})
B^{\bm{\alpha-\beta}}\\
= & 
\sum_{\bm\beta\in\mathbb{T}_r^d,\ \bm\alpha-\bm\beta\geq0,\ \bm\gamma_{\to i, j}-\bm\beta_{\to i, j}\geq0}
B^{\bm \alpha-\bm\beta}
\frac{k!}{(k-r)!} \sum_{\bm\sigma\in\mathbb T_{\beta_0}^1,
\bm\gamma_{\to i, j}-\bm\sigma-\bm\beta_{\to{i,j}}\geq0}&\frac{(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}})!}{\bm{\sigma}!(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}}-\bm\sigma)!}
\bm E_{ij}^{\bm \sigma\gets(\bm \gamma_{\to i, j}-\bm \sigma-\bm \beta_{\to i, j })}
\end{aligned}
$$

因此当 $\bm\alpha-\bm\beta=\bm\gamma_{\{ij\}}$ 且 $\bm\gamma_{\to
i,j}-\bm\beta_{\to i,j}\geq 0$ 时，
$$
\phi^{\bm \gamma_{\{i,j\}}}(\frac{\partial^r B^{\bm \alpha}}{\partial \bm{N}_{e_{ij}}^{\bm\gamma_{\to i, j}}} )
 = 
\frac{k!}{(k-r)!} \sum_{\bm\sigma\in\mathbb T_{\beta_0}^1,
\bm\gamma_{\to i, j}-\bm\sigma-\bm\beta_{\to{i,j}}\geq0}
\frac{(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}})!}{\bm{\sigma}!(\bm\gamma_{\to i, j}-\bm\beta_{\to{i,j}}-\bm\sigma)!}
\bm E_{ij}^{\bm \sigma\gets(\bm \gamma_{\to i, j}-\bm \sigma-\bm \beta_{\to i, j })}
$$
得证(注意：当 $\bm\alpha-\bm\beta=\bm\gamma_{\{ij\}}$ 时，
$\bm\alpha_{\to i,j} = \bm\beta_{\to i,j}$)。














</br>
</br>
</br>
</br>
</br>
</br>
</br>
</br>
