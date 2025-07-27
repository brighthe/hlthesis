# Complex from Complex 阅读笔记

## 一、复形简介
1. 光滑 de Rham 复形
    $$
    0\longrightarrow
    C^\infty\Lambda^0\xrightarrow{d^0}
    C^\infty\Lambda^1\xrightarrow{d^1}
    \cdots\xrightarrow{d^{n-1}}
    C^\infty\Lambda^n\longrightarrow0.
    $$
    三维情况的光滑复形：
    $$
    0\longrightarrow
    C^\infty(\Omega)\xrightarrow{\mathrm{grad}}
    C^\infty(\Omega;\mathbb{R}^3)\xrightarrow{\mathrm{curl}}
    C^\infty(\Omega;\mathbb{R}^3)\xrightarrow{\mathrm{div}}
    C^\infty(\Omega)\longrightarrow0,
    $$
2. Sobolev de Rham 复形，对于任意的实数$q$
    $$
    0\longrightarrow
    H^q\Lambda^0\xrightarrow{d^0}H^{q-1}
    \Lambda^1\xrightarrow{d^1}
    \cdots\xrightarrow{d^{n-1}}
    H^{q-n}\Lambda^n\longrightarrow0.
    $$
    Sobolev 复形中的微分算子是有界算子因此整体称为有界 Hilbert 复形。
3. $L^2$ de Rham 复形
    $$
    0\longrightarrow
    L^2\Lambda^0\xrightarrow{d^0}
    L^2\Lambda^1\xrightarrow{d^1}
    \cdots
    \xrightarrow{d^{n-1}}
    L^2\Lambda^n\longrightarrow
    0.
    $$
    $L^2$de Rham 复形中微分算子是无界的，但是是几乎闭的，也是稠密定义的。定义domain：
    $$
    H \Lambda^k=\{u\in L^2\Lambda^k: d^ku\in L^2\Lambda^{k+1}\}.
    $$
    那么得到一个 $L^2$de Rham 复形的 domain 复形
    $$
    0\longrightarrow H\Lambda^0\xrightarrow{d^0}
    H\Lambda^1\xrightarrow{d^1}\cdots\xrightarrow{d^{n-1}}
    H\Lambda^n\longrightarrow0.
    $$
    这个复形是一个有界 Hilbert 复形。
    - 对于 $L^2$ 上的无界算子 $D$，定义图范数 
      $\Vert u\Vert_{D}^2 := \Vert u\Vert^2 + \Vert Du\Vert^2$.
    - 稠密算子的对偶算子是闭算子。
    - Sobolev 复形和 $L^2$ de Rham 复形中的微分算子都是闭算子，这个性质在 FEEC
      中非常重要，这可以推导出：
      - Poincare 不等式
      - Hodge 分解
      - Hodge Laplacian 边值问题的适定性
    - 以上所有 de Rham 复形的同调空间都是有限维的，而且彼此同构，所以可以选择
      $C^\infty$同调来表示所有的同调空间。
    - 如果 $\Omega$ 是可缩的，那么以上所有复形的上同调空间为 $0$.
    - $L^2$domain 复形中，所有空间可以做正则分解。
    - $L^2$de Rham 复形满足紧性。
4. Hilbert 复形：
    $$
    \cdots\to Z^{k-1}\xrightarrow{D^{k-1}}Z^{k}\xrightarrow{D^{k}}Z^{k+1}\to\cdots,
    $$
    其中 $Z^i$ 是一系列的空间 $D^i$ 是 $Z^i$ 到 $Z^{i+1}$ 的线性算子，
    其中 $D^i$ 的定义域是 $Z^i$ 的某个子空间，且满足 
    $\mathcal{R}(D^k)\subset \mathcal{N}(D^{k+1})$，
    满足这个条件成这一系列空间为 **复形**,
    若 $Z^k$ 是 Hilbert 空间, 且 $D^k$ 是闭的且稠密定义的，
    则称之为 **Hilbert 复形**。
    - 复形的同调空间定义为 $\mathcal{H}^k := \mathcal{N}(D^k)/\mathcal{R}(D^{k-1})$. 
    - 算子 $D^k$ 的核空间一定是闭的，但是像空间不一定是闭的。 所以不能写:
      $$
      \mathcal{H}^k = \mathcal{N}(D^k)\cap\mathcal{R}(D^{k-1})^{\perp}.
      $$
      所以 $\mathcal{H}^k$ 不一定是 Hilbert 空间(不一定完备)。
    - $\mathcal{H}^k$ 是 $D^{k-1}$ 的像空间在 $D^k$
        的核空间的补空间。补空间有很多定义，只要满足:
        $$
        \mathcal{N}(D^k) = \mathcal{R}(D^{k-1})\oplus\mathcal{H}^k.
        $$
    - 定义空间 $H^q\otimes \mathbb{E}$ 为 $H^q$ 上取值于 $\mathbb{E}$ 的函数空间。
        其范数定义为 $\Vert u\Vert_{q}$, $q=0$ 时忽略 $q$。
    - 一个 $\mathcal{R}(D)$ 是闭的的充分条件: $\mathcal{H}^k$ 是有限维的.
    - 对于一个算子 $T: X\to Y$, 其对偶算子的定义域为: 
      $$ D(T^*) = \{w\in Y: \exists c_w > 0, s.t. |\langle Tu, w\rangle_Y|\leq
      c_w\Vert u\Vert_X, \forall u\in D(T)\}. 
      $$
      即 $D(T^*)$ 中的元素 $w$ 可以定出来一个 $D(T)$ 的有界线性泛函：
      $$
      L_w(u) := \langle Tu, w\rangle_Y = \langle u, T^*w\rangle_X.
      , \forall u\in D(T).
      $$
      即 $T^*w$ 根据 $X$ 的内积定义了一个 $D(T)$ 的有界线性泛函。


































































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











