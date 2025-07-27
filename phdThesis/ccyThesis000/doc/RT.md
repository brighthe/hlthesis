# RT 元分解
定义 $\mathbb{P}_k$ 为 $k$ 次多项式空间，$\mathbb{H}_k$
为 $k$ 次齐次多项式空间. 对于 $n$ 维单形 $\tau$, 
$$RT_k(\tau) = \{v \in (\mathbb{P}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau) :
\bm{v}\cdot\bm{n_f}\in \mathbb{P}_k(f)\ \forall f \in \partial\tau\}$$

对 $(\mathbb{P}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau)$ 作如下分解:
$$
(\mathbb{P}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau) = 
(\mathbb{P}_{k-1}(\tau))^n\oplus(\mathbb{H}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau)
$$
其中 $(\mathbb{P}_{k-1}(\tau))^n$ 的基函数可以由 $t-n$ 分解得出。

对于 $(\mathbb{H}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau)$, 我们令
$\tau$ 的顶点为 $\{\bm{p}_i, \bm{p}_1, ..., \bm{p}_n\}$, 面为
$\{{f}_i, {f}_1, ..., {f}_n\}$ 满足 $f_i\cap\bm{p}_i = \empty$.
定义 $B_k$ 为$\mathbb{H}_k$ 的 Bernstein 基函数（或者是Lagrange插值基函数）。
那么:
$$
(\mathbb{H}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau) = 
\mathrm{span}\{\cup_{i=0}^n(\bm{x-p}_i)B_k\} 
$$
证明：维数相等，且$(\bm{x-p}_i)B_k \in 
(\mathbb{H}_k(\tau))^n\oplus\bm{x}\mathbb{H}_k(\tau)$ 所以等式成立。

另外可以证明当 $i\not = j$, $(\bm{x-p}_i)B_k \cdot \bm{n}_{f_j} = 0$.

