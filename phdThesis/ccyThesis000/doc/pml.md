# 张量角度的 PML 边界条件理论
## 一、映射与与推前拉回映射
考虑两个微分流形 $M, N$, $\phi$ 是一个 $M$ 到 $N$ 的一个微分同胚, 在 $M$
上有一个坐标系 $\{x_i\}_{i=0}^{n-1}$, 其坐标基矢量为 $\{\partial x_i\}_{i=0}^{n-1}$, 
对偶基矢量为 $\{ \mathrm{d} x_i\}_{i=0}^{n-1}$, 
$N$ 上有一组坐标系 $\{y_i\}_{i=0}^{n-1}$ 其坐标基矢量为 $\{ \partial y_i\}_{i=0}^{n-1}$, 
对偶基矢量为 $\{ \mathrm{d} y_i\}_{i=0}^{n-1}$, 
现在考虑推前拉回映射的张量表示.

### 1. 函数的拉回映射
对于任意 $f \in C^{\infty}(N)$, 其拉回为 $\phi^* f  = f \circ\phi\in C^{\infty}(M)$,
两个流形上的梯度为:
$$
\nabla_y f = \frac{\partial f}{\partial y_i} \mathrm{d}y_i, 
\quad \nabla_x \phi^* f = \frac{\partial f}{\partial x_i} \mathrm{d}x_i
$$

### 2. 向量场的推前映射
记 $\phi_* \partial x_i = F_{ji} \partial y_j$ 则:
$$
\nabla_x \phi^* f \cdot \partial x_i = \nabla_y f \cdot \phi_* \partial x_i
$$

即:
$$
\begin{aligned}
\frac{\partial f}{\partial x_i}
& = \frac{\partial f}{\partial y_j} \mathrm{d}y_jF_{ki} \partial y_k\\
& = \frac{\partial f}{\partial x_l} \frac{\partial x_l}{\partial y_j} F_{ji}
\end{aligned}
$$

所以 $F_{ji} = \frac{\partial y_j}{\partial x_i}$, 定义张量 $\bm F = \frac{\partial
y_j}{\partial x_i} \partial y_j \otimes \mathrm{d}x_i$, 那么对于任意光滑向量场 
$\bm v = v_i \partial x_i, v_i \in C^{\infty}(F)$, 
$$
\phi_* \bm v = \bm F\bm v
$$

### 3. 对偶矢量场(1-形式)的拉回映射
记 $\phi^* \mathrm{d}y_i = A_{ji} \mathrm{d}x_j$ 则:
$$
\mathrm{d}y_i(\phi_*\partial x_j) = (\phi^* \mathrm{d}y_i)( \partial x_j)
$$

即:
$$
\begin{aligned}
F_{kj} \delta_{ki} & =  A_{mi} \delta_{mj}
\end{aligned}
$$

所以 $A_{ji} = F_{ij} = \frac{\partial y_i}{\partial x_j}$, 
 那么对于任意光滑对偶矢量场 
$\bm \omega = \omega_i \mathrm{d} y_i, v_i \in C^{\infty}(N)$, 
$$
\phi^* \bm \omega = \bm F^T\bm \omega
$$
其中 $\bm F^T = \frac{\partial y_j}{\partial x_i} \mathrm{d}x_i \otimes \partial y_j$, 

### 3. 2-形式的 Hodge star 的拉回映射
$$
\phi^*( \mathrm{d}y_i \wedge \mathrm{d}y_j) = (\phi^* \mathrm{d}y_i) \wedge 
(\phi^* \mathrm{d}y_j) = F_{il}F_{jm} \mathrm{d}x_l \wedge \mathrm{d}x_m
$$
对于任意 2-形式 $\bm \omega = \omega_{ij} \mathrm{d}y_i \wedge \mathrm{d}y_j$,
我们讨论其 Hodge star 的拉回映射的张量表示:
#### 3.1 $n=2$
$$
\star\bm \omega = \frac{1}{2}\in_{ij} \omega_{ij}
$$
$$
\begin{aligned}
\star \phi^* \bm \omega 
& = \frac{1}{2}\in_{lmn} F_{il}F_{jm} \omega_{ij} \mathrm{d}x_n\\
& = \frac{1}{2}|\bm F|\in_{ij}\omega_{ij}\\
& = |\bm F|\star \bm \omega
\end{aligned}
$$
所以:
$$
\phi^* (\star \bm \omega) =  |\bm F|^{-1}\star \phi^*\bm \omega
$$
若 $\bm \omega = \mathrm{d}\bm \eta = \frac{ \partial \eta_i}{ \partial y_k}
\mathrm{d} y_k \wedge \mathrm{d}y_i$
那么有:
$$
\phi^* (\star \mathrm{d}\bm \eta) = |\bm F|^{-1} \star \phi^* \mathrm{d}\bm
\eta = |\bm F|^{-1} \star \mathrm{d}\phi^*\bm\eta
$$

#### 3.2 $n=3$
$$
\star\bm \omega = \frac{1}{2}\in_{ijk} \omega_{ij} \mathrm{d}y_k
$$
$$
\begin{aligned}
\star \phi^* \bm \omega 
& = \frac{1}{2}\in_{lmn} F_{il}F_{jm} \omega_{ij} \mathrm{d}x_n\\
& = \frac{1}{2}|\bm F|\in_{ijk}F^{-1}_{nk} \omega_{ij} \mathrm{d}x_n
\end{aligned}
$$

其中 $F^{-1}_{nk} = \frac{ \partial x_n}{ \partial y_k},$
令 $\bm G = \frac{ \partial x_n}{ \partial y_k} \mathrm{d}x_n \otimes \partial y_k$
则有:
$$
\star\phi^* \bm \omega = |\bm F|\bm G \star \bm \omega
, \quad 
\phi^* (\star \bm \omega) =  \bm A\star \phi^*\bm \omega
$$
其中 $\bm A = |\bm F|^{-1}\bm F^T \bm G^{-1}$.
若 $\bm \omega = \mathrm{d}\bm \eta = \frac{ \partial \eta_i}{ \partial y_k}
\mathrm{d} y_k \wedge \mathrm{d}y_i$
那么有:
$$
\phi^* (\star \mathrm{d}\bm \eta) = \bm A \star \phi^* \mathrm{d}\bm
\eta = \bm A \star \mathrm{d}\phi^*\bm\eta
$$

## 二、 PML 层的方程
考虑 PML 层 $M$, 某个映射 $\phi$ 将 $M$ 映射到 $N = \phi(M)$, 在 $N$ 上有如下时谐 
Maxwell 方程:
$$
\nabla\times \mu^{-1}\nabla\times \bm E - k^2 \bm E = 0
$$
其中 $k^2 = \omega^2\epsilon$, 将其写为外微分形式为($\nabla \times \bm E = \star \mathrm{d}\bm E$ ):
$$
\star \mathrm{d}\mu^{-1}\star \mathrm{d} \bm E - k^2 \bm E = 0
$$

### 1. $n=2$
将该方程的拉回到 $M$ 得到:
$$
\begin{aligned}
0 = \phi^*(\star \mathrm{d}\mu^{-1}\star \mathrm{d} \bm E - k^2 \bm E) 
& = \phi^*(\star( \mathrm{d}\mu^{-1}\star \mathrm{d}\bm E)) - k^2\phi^*\bm E\\
& = |\bm F|^{-1} \star \mathrm{d}\mu^{-1}\phi^* (\star \mathrm{d}\bm E) - k^2
\phi^* \bm E\\
& = |\bm F|^{-1} \star \mathrm{d}\mu^{-1}|\bm F|^{-1}
\star\mathrm{d}\phi^* \bm E - k^2 \phi^* \bm E\\
\end{aligned}
$$
记 $\tilde{\bm E} = \phi^*\bm E$, 整理上式得到:
$$
\star \mathrm{d}\mu^{-1}|\bm F|^{-1}
\star\mathrm{d}\tilde{\bm E} - k^2 |\bm F| \tilde{\bm E} = 0
$$
上式写成 $\mathbb{R}^2$ 空间矢量场论的形式为:
$$
\nabla\times \mu^{-1}|\bm F|^{-1} \nabla \times \tilde{\bm E} - k^2|\bm F| \tilde{\bm E} = 0
$$

其中 $\bm F$ 是映射 $\phi$ 的 Jacobi 矩阵
### 2. $n=3$
将该方程的拉回到 $M$ 得到:
$$
\begin{aligned}
0 = \phi^*(\star \mathrm{d}\mu^{-1}\star \mathrm{d} \bm E - k^2 \bm E) 
& = \phi^*(\star( \mathrm{d}\mu^{-1}\star \mathrm{d}\bm E)) - k^2\phi^*\bm E\\
& = \bm A \star \mathrm{d}\mu^{-1}\phi^* (\star \mathrm{d}\bm E) - k^2
\phi^* \bm E\\
& = \bm A \star \mathrm{d}\mu^{-1}\bm A 
\star\mathrm{d}\phi^* \bm E - k^2 \phi^* \bm E\\
\end{aligned}
$$
记 $\tilde{\bm E} = \phi^*\bm E$, 整理上式得到:
$$
\star \mathrm{d}\mu^{-1}\bm A 
\star\mathrm{d}\tilde{\bm E} - k^2 \bm A^{-1} \tilde{\bm E} = 0
$$
上式写成 $\mathbb{R}^3$ 空间矢量场论的形式为:
$$
\nabla\times \mu^{-1}\bm A \nabla \times \tilde{\bm E} - k^2\bm A^{-1} \tilde{\bm E} = 0
$$

其中 $\bm A = |\bm F|^{-1} \bm F^T \bm F$, $\bm F$ 是映射 $\phi$ 的 Jacobi 矩阵

### 3. 圆形 PML 层
考虑圆形区域 $\Omega = \{(x, y) : x^2 + y^2 < r_0^2\}$, 构造 PML 区域 
$\Omega_{P} = \{(x, y) : r_0^2 \leq x^2 + y^2 < r_1^2\}$

定义拉伸函数:
$$
f(r) = r + \omega^{-1} s ( \frac{r-r_0}{r_1-r_0})^p(1- \mathrm{i})
$$
$$
f(r) = r + \mathrm{i} \frac{s}{\omega\sqrt{\mu \epsilon}} ( \frac{r-r_0}{r_1-r_0})^p
$$
映射 $\phi$
$$
\left\{
\begin{aligned}
  \tilde{r} & = f(r)\\
  \tilde{\theta} & = \theta\\
\end{aligned}
\right.
$$
将 $\Omega_{P}$ 映射到 $\tilde{\Omega}_P = \phi(\Omega_P)$ 笛卡尔坐标为:
$$
\left\{\begin{aligned}
x = r\cos(\theta)\\
y = r\sin(\theta)
\end{aligned}\right.
, \quad 
\left\{\begin{aligned}
\tilde{x} = \tilde{r}\cos(\tilde{\theta})\\
\tilde{y} = \tilde{r}\sin(\tilde{\theta})
\end{aligned}\right.
$$
那么映射 $\phi$ 的 Jacobi 矩阵为：
$$
\begin{aligned}
\boldsymbol{F} 
& = \frac{\partial(\tilde{x}, \tilde{y})}{\partial(x, y)}
= \frac{\partial(\tilde{x}, \tilde{y})}{\partial(\tilde{r}, \tilde{\theta})}
\frac{\partial(\tilde{r}, \tilde{\theta})}{\partial(r, \theta)}
\frac{\partial(r, \theta)}{\partial(x, y)}\\
& = 
\begin{pmatrix}
\cos(\tilde{\theta}) &-\tilde{r}\sin(\tilde{\theta})\\
\sin(\tilde{\theta}) & \tilde{r}\cos(\tilde{\theta})\\
\end{pmatrix}
\begin{pmatrix}
f'(r) & 0\\
0     & 1
\end{pmatrix}
\begin{pmatrix}
\cos(\theta) & \sin(\theta)\\
-\frac{1}{r}\sin(\theta) & \frac{1}{r}\cos(\theta)
\end{pmatrix}
\end{aligned}
$$
Jacobi 行列式为:
$$
|\bm F| = f(r)f'(r)r^{-1}
$$

### 4. 球形 PML 层
考虑球形区域 $\Omega = \{(x, y, z) : x^2 + y^2 + z^2 < r_0^2\}$, 构造 PML 区域 
$\Omega_{P} = \{(x, y, z) : r_0^2 \leq x^2 + y^2 + z^2 < r_1^2\}$

定义拉伸函数:
$$
f(r) = r + \omega^{-1} s ( \frac{r-r_0}{r_1-r_0})^p(1- \mathrm{i})
$$
$$
f(r) = r + \mathrm{i} \frac{s}{\omega\sqrt{\mu \epsilon}} ( \frac{r-r_0}{r_1-r_0})^p
$$
映射 $\phi$
$$
\left\{
\begin{aligned}
  \tilde{r} & = f(r)\\
  \tilde{\theta} & = \theta\\
  \tilde \varphi & = \varphi
\end{aligned}
\right.
$$
将 $\Omega_{P}$ 映射到 $\tilde{\Omega}_P = \phi(\Omega_P)$ 笛卡尔坐标为:
$$
\left\{\begin{aligned}
x & = r\cos(\theta)\sin(\varphi)\\
y & = r\sin(\theta)\sin(\varphi)\\
z & = r\cos(\varphi)
\end{aligned}\right.
, \quad 
\left\{\begin{aligned}
\tilde{x} & = \tilde{r}\cos(\tilde{\theta})\sin(\tilde \varphi)\\
\tilde{y} & = \tilde{r}\sin(\tilde{\theta})\sin(\tilde\varphi)\\
z & = \tilde r\cos(\tilde\varphi)
\end{aligned}\right.
$$
那么映射 $\phi$ 的 Jacobi 矩阵为：
$$
\begin{aligned}
\boldsymbol{F} 
& = \frac{\partial(\tilde{x}, \tilde{y}, \tilde z)}{\partial(x, y, z)}
  = \frac{\partial(\tilde{x}, \tilde{y}, \tilde z)}
         {\partial(\tilde{r}, \tilde{\theta},\tilde \varphi)}
    \frac{\partial(\tilde{r}, \tilde{\theta},\tilde \varphi)}
         {\partial(r, \theta, \varphi)}
    \frac{\partial(r, \theta, \varphi)}
         {\partial(x, y, z)}\\
& = 
\begin{pmatrix}
\cos(\tilde{\theta})\sin(\tilde{\varphi}) & 
-\tilde{r}\sin(\tilde{\theta})\sin(\tilde{\varphi}) & 
\tilde{r}\cos(\tilde\theta)\cos(\tilde{\varphi})\\

\sin(\tilde{\theta})\sin(\tilde{\varphi}) & 
\tilde{r}\cos(\tilde{\theta})\sin(\tilde{\varphi}) & 
\tilde{r} \sin(\tilde\theta)\cos(\tilde\varphi)\\

\cos(\tilde{\varphi}) & 
0 & 
-\tilde{r} \sin(\tilde\varphi)\\
\end{pmatrix}
\begin{pmatrix}
f'(r) & 0 & 0\\
0     & 1 & 0\\
0     & 0 & 1\\
\end{pmatrix}\\

& \times 
\begin{pmatrix}
\cos(\theta)\sin(\varphi) & \sin(\theta)\sin(\varphi) & \cos(\varphi)\\
-\frac{\sin(\theta)}{r\sin(\varphi)} & \frac{\cos(\theta)}{r\sin(\varphi)} & 0\\
\frac{1}{r}\cos(\theta)\cos(\varphi) & \frac{1}{r}\sin(\theta)\cos(\varphi) &  -\frac{1}{r}\sin(\varphi)
\end{pmatrix}\\
& = 
\begin{pmatrix}
f'(r)\cos(\tilde{\theta})\sin(\tilde{\varphi}) & 
-\tilde{r}\sin(\tilde{\theta}) & 
\tilde{r}\cos(\tilde\theta)\cos(\tilde{\varphi})\\

f'(r)\sin(\tilde{\theta})\sin(\tilde{\varphi}) & 
\tilde{r}\cos(\tilde{\theta})& 
\tilde{r} \sin(\tilde\theta)\cos(\tilde\varphi)\\

f'(r)\cos(\tilde{\varphi}) & 
0 & 
-\tilde{r} \sin(\tilde\varphi)\\
\end{pmatrix}
\\
& \times
\begin{pmatrix}
\cos(\theta)\sin(\varphi) & \sin(\theta)\sin(\varphi) & \cos(\varphi)\\
-\frac{1}{r} \sin(\theta)& \frac{1}{r}\cos(\theta) & 0\\
\frac{1}{r}\cos(\theta)\cos(\varphi) & \frac{1}{r}\sin(\theta)\cos(\varphi) &  -\frac{1}{r}\sin(\varphi)
\end{pmatrix}\\
\end{aligned}
$$
Jacobi 行列式为:
$$
|\bm F| = f(r)f'(r)r^{-1}
$$
## 三、算例  





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





