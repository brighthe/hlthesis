本文中我们使用递归的方法定义了任意维任意阶任意 $m\in \mathbb{N}$ 的 $H^m$
协调虚单元。一维的 $H^m$ 协调虚单元本质上是 $C^{m-1}$ 连续的有限元:
$$
V_k^m(K):=\{v\in H^m(K):v^{(2m)}\in\mathbb{P}_{k-2m}(K)\}=\mathbb{P}_{\max\{k,2m-1\}}(K).
$$
对于 $n$ 维多面体 $K$，假设 $n-1$ 维多面体上的 $H^l$
协调虚单元已经被定义 $l=1, 2, ..., m$，定义 $K$ 上空间 
$$
\begin{aligned}
\widetilde{V}_k^m(K)&:=\big\{ v\in H^m(K):  (-\Delta)^mv\in \mathbb P_{k}(K),  (\nabla^jv)|_{\mathcal S_K^r}\in H^1(\mathcal S_K^r;\mathbb S_{n}(j)) \textrm{ for } 0\leq j\leq m-1, \\
&\quad\;\;\frac{\partial^{|\alpha|}v}{\partial \nu_F^{\alpha}}\Big|_F\in V_{k-|\alpha|}^{m-|\alpha|}(F)\quad\forall~F\in\mathcal F^{r}(K), 1\leq r\leq n-1, \alpha\in A_r, \textrm{ and } |\alpha|\leq m-1\big\}.    
\end{aligned}
$$
定义如下自由度：
$$
\begin{aligned}h_{K}^{j}\nabla^{j}v(\delta)&\quad \forall \delta\in\mathcal{F}^{n}(K),
j=0,\ldots,m-1,\\
\frac{h_{K}^{|\alpha|}}{|F|}\left(\frac{\partial^{|\alpha|}v}
{\partial\nu_{F}^{\alpha}},q\right)_{F}&\quad\forall
q\in\mathbb{P}_{k-2m+|\alpha|}(F),F\in\mathcal{F}^{r}(K), r=1,\ldots\\
&\quad \quad \alpha\in A_{r}, \mathrm{and} |\alpha|\leq m-1,\\
\frac{1}{|K|}(v,q)_{K}&\quad \forall q\in\mathbb{P}_{k-2m}(K).\end{aligned}
$$
$H^m$
那么可以定义到 $\widetilde{V}_k^m(K)$ 到 $\mathbb{P}_k(K)$ 的 $H^m$ 投影算子
$\Pi_k^m$:
$$
\begin{aligned}
(\nabla^{m}\Pi_{k}^{K}v,\nabla^{m}q)_{K}& =(\nabla^{m}v,\nabla^{m}q)_{K}\quad\forall
q\in\mathbb{P}_{k}(K),\\
\sum_{\delta\in\mathcal{F}^{n}(K)}(\nabla^{j}\Pi_{k}^{K}v)(\delta) & =\sum_{\delta\in\mathcal{F}^{n}(K)}(\nabla^{j}v)(\delta),\quad
j=0,\ldots,m-1\end{aligned}
$$
最终，我们定义 $H^m$ 协调虚单元空间
$$
V_k^m(K):=\{v\in\widetilde{V}_k^m(K):(v-\Pi_k^Kv,q)_K=0\quad\forall q\in\mathbb{P}_{k-2m}^\perp(K)\}.
$$
我们证明了逆不等式，范数等价性等。定义了如下稳定化项：
$$
S_{K}(w,v):=\sum_{r=1}^{n}\sum_{F\in\mathcal{F}^{r}(K)}\sum_{|\alpha|\leq m-1}h_{K}^{r+2|\alpha|-2m}\left(Q_{k-2m+|\alpha|}^{F}\frac{\partial^{|\alpha|}w}{\partial\nu_{F}^{\alpha}},Q_{k-2m+|\alpha|}^{F}\frac{\partial^{|\alpha|}v}{\partial\nu_{F}^{\alpha}}\right)_{F}
$$
给出了多重调和方程的 $H^m$ 协调虚单元方法，数值实验验证了方法的有效性。


