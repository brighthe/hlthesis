# Literature research
### 1. Partially Discontinuous Nodal Finite Elements for H(curl)and H(div)
**总结**：本文基于宏元，构造了四面体上顶点和边上连续的 $H(curl)$ 协调的有限元和 $H(div)$
协调的有限元，并分别给出了基函数，其中基函数是有一个标量乘以某个向量的形式定义的，该向量与四面体的某个子单形的几何有关。

**疑问**:  为什么文章中说 : We cannot further relax the vertex and edge
continuity, which is required to ensure the existence of nodal bases. 

### 2. GEOMETRIC DECOMPOSITIONS AND LOCAL BASES FOR SPACES OF FINITE ELEMENT DIFFERENTIAL FORMS
**总结**：本文基于微分形式的几何分解给出了 BDM 元，RT 元，第一类第二类 Nedelec
元的基函数, 但是没有给出基函数的对偶自由度.

### 3. Computational Bases for $RT_k$ and $BDM_k$ on Triangles
**总结**：本文基于 Poila 映射, 给出了三角形上 $RT$ 元和 $BDM$ 元的基函数,
首先给出参考单元上的基函数，然后使用协变或逆变的 poila
映射将参考单元的基函数映射的真实单元上，得到真实单元上的基函数，

### 4. Nodal finite element de Rham complexes
**总结** 构造了具有更高连续性的 de Rham 复形，包括如: 在顶点连续的 $H(div)$ 元, 
在顶点上 $C^1$ 连续在边上 $C^0$ 连续的 $H(curl)$ 元,
并给出了节点基函数的构造。其中基函数是由标量的节点基函数乘以一个向量构造的。

### 5. EFFICIENT ASSEMBLY OF H(div) AND H(curl) CONFORMING FINITE ELEMENTS
**总结**：本文基于 Poila 映射给出了单纯性网格上 $H(div)$ 元和 $H(curl)$
元的基函数，并给出了矩阵的快速组装方法。

### 6. Bernstein–Bézier bases for tetrahedral finite elements
**总结** 本文给出了四面体网格上 de rahm 复形对应的有限元基函数包括 RT 元，BDM 元 $H(curl), H(div), L2$ 的基函数，且能保证非梯度的
$H(curl)$ 协调的基函数的旋度是 $H(div)$ 的基函数的一部分, 非旋度的 $H(div)$
协调的基函数的散度是 $L^2$ 有限元的基函数的一部分.

### 7. A Bernstein–Bzier Basis for Arbitrary Order Raviart–Thomas Finite Elements
**总结** 基于 Bernstein 多项式给出了 RT 元的基函数.

### 8. A Well-Conditioned Hierarchical Basis for Triangular H(curl)-Conforming Elements 
### 9. ON THE CONSTRUCTION OF WELL-CONDITIONED HIERARCHICAL BASES FOR TETRAHEDRAL H(curl)-CONFORMING N ́ED ́ELEC ELEMEN
### 10. On the Construction of Well-Conditioned Hierarchical Bases for H(div)-Conforming Rn Simplicial Elements

### 11. Hierarchic finite element bases on unstructured tetrahedral meshes

**总结** 建立了 $H(div), H(curl)$ 的分层基函数.

### 12. FreeFEM

**总结** 只有 RT0, RT1, RT2 空间

### 13. Fenics, MFEM 
**总结** 使用 poila 映射构造 Hcurl, Hdiv 元的基函数

### 14. LOW-COMPLEXITY FINITE ELEMENT ALGORITHMS FOR THE DE RHAM COMPLEX ON SIMPLICES
**总结** 接着 2 做了快速算法.


1, 4 都是节点基函数

在物理仿真应用中，物理量的类型可能是标量(密度)，向量(速度，电场)或张量(应力，应变)，
对于向量形的

H1协调的向量型有限元， Hdiv 元和 Hcurl 元在流体仿真和电磁仿真领域广泛使用，
Langrange 元，BDM 元和第二类 Nedelec 元(简写为 $ND^2$)是常见的 H1 协调，Hdiv
协调， Hcurl 协调的有限元, 且它们有一个共同点, 他们的形函数空间都是 $\mathbb
P_k^n$, 其中 $k$ 是多项式次数，$n$ 是空间维数, 但是由于空间的迹的定义不同，
他们的基函数也不相同。

本文提出了单纯形网格上 Langrange 元，BDM 元和第二类 Nedelec 元的节点基函数，
这三类有限元分别是H1， Hdiv 和 Hcurl 协调的有限元，他们的形函数空间都是 $\mathbb
P_k^n$, 其中 $k$ 是多项式次数，$n$ 是空间维数, 但是由于空间的迹的定义不同，
他们的基函数也不相同。 

对于 Lagrange 有限元，常用的节点型的基函数具有形式简单，
容易求值等优点，而现有的 BDM 元和 $ND^2$ 元的基函数却复杂很多,
其经典的构造方式是使用 poila 映射，先在标准单元上构造基函数，
然后使用协变(保切向，BDM 元)或逆变(保法向，$ND^2$ 元)
poila 映射，将基函数映射到实际单元上得实际单元上的基函数。这方面的工作可以在[5,
3] 中找到详细的介绍。同时这也是一些开源软件如 MFEM，Fenics 的实现方式。

在[2]中 Arnold 提出了多项式空间以及各阶微分形式的几何分解, 并给出了基于
Bernstein 多项式的基函数。由于 Bernstein 基函数相关的快速算法的提出， [6, 7]
中同样提出了基于 Bernstein 多项式的基函数，以及相关矩阵的快速组装，
无矩阵组装算法。 特别地, 在 [6] 中，构造的基函数有一些有趣的性质: $H(curl)$
协调的有限元的基函数的旋度恰好是 $H(div)$ 协调的有限元空间的基函数的一部分, 
$H(div)$ 协调的有限元的基函数的旋度恰好是 $L^2$ 有限元空间的基函数的一部分。
另外, 关于 $H(curl)$协调的分层基函数可以在[8, 9, 10, 11] 中找到.

前述的基函数构造方式比较复杂，因此一些学者开始尝试做一些简单的构造方式，在[5]
中，H 等人使用标量的节点有限元方法乘以某个向量的方式，
给出了在顶点和边上都连续的 $H(div), H(curl)$
协调有限元，但是在 [1] 中，经过数值验证发现使用这种 $H(curl)$ 协调的有限元离散时谐 
Maxwell 方程会出现伪解，因此一种基于宏元的在边和顶点上完全连续的 $H(curl)$
协调有限元被提出。但是，对于 $BDM$ 元和 $ND^2$ 元的节点基函数还不存在. 


























