# 张量有限元方法的研究计划

## 1. 研究背景
有限元方法是一种求解偏微分方程的数值方法[1]，于20世纪50年代由Courant等人提出，
有限元方法的基本思想是将复杂求解区域分解为简单几何单元
（如三角形、四边形、四面体、六面体等）的组成的网格，在网格上建立具备合适连续性的
分片多项式函数组成的有限元空间，
然后基于 Galerkin 方法，将偏微分方程离散为代数方程组，通过求解代数方程组得到数值解。
在应用方面因为有限元方法可以在复杂的物理区域上建立数值模型，适用于各种复杂的边界条件，
所以在结构力学、传热学、流体力学、电磁场计算等领域都有广泛应用。

有限元方法因为其完备且优雅的数学基础，在理论方面也得到广泛研究，特别是最近关于
有限元外微分的研究，使用 Hilbert 复形将一类的问题联系起来，通过研究 Hilbert
复形给出一类问题的高效数值方法。 有限元外微分的主要思想是 构造与 Hilbert
复形对应的有限元复形, 通过它们之间的上链投影算子，
建立了离散有限元复形的稳定性，收敛性理论，
且有限元外微分技术易于推广到高次元情况，为构造稳定且精确的有限元方法提供了理论框架。
基于有限元外微分，各种微分复形及其对应的有限元复形被提出，如 de Rham 复形、
弹性复形、Hessian 复形、divdiv 复形等[1]。

### （1）de Rham 复形

$L^2$ de Rham 复形是最常见的一种 Hilbert 复形，
考虑 3 维空间的 de Rham 复形:
$$
\mathbb{R} \xrightarrow{\subset} H^1(\Omega) \xrightarrow{\nabla} H(\text{curl}, \Omega)
\xrightarrow{\nabla \times} H(\text{div}, \Omega) \xrightarrow{\nabla \cdot}
L^2(\Omega) \rightarrow 0 
$$
其中 $\Omega$ 是 $\mathbb{R}^3$ 中的多面体区域.
de Rham 复形可用于电磁问题、流体力学问题的数值模拟, 
将 $\Omega$ 剖分为网格 $\mathcal{T}$，构造 
$\mathcal{T}$ 有限元复形即构造具有对应协调性的有限元，根据如下格林公式，
$$
\begin{aligned}
\int_{\Omega} \nabla u \cdot \boldsymbol{v} \, dx & = 
- \int_{\Omega} u \nabla\cdot \boldsymbol{v} \, dx 
+ \int_{\partial \Omega} u \boldsymbol{v}\cdot \boldsymbol{n} \, ds\\
\int_{\Omega} \nabla \times \boldsymbol{u}\cdot \boldsymbol{v} \, dx & =
\int_{\Omega} \boldsymbol{u} \cdot \nabla \times \boldsymbol{v} \, dx +
\int_{\partial \Omega} (\boldsymbol{n}\times \boldsymbol{u}) \cdot \boldsymbol{v} \, ds\\
\int_{\Omega} \nabla \cdot \boldsymbol{u} v \, dx & =
- \int_{\Omega} \boldsymbol{u} \cdot \nabla v \, dx +
\int_{\partial \Omega} \boldsymbol{u} \cdot \boldsymbol{n} v \, ds
\end{aligned}
$$
可知对于分片多项式函数 $u$ (向量函数 $\boldsymbol{u}$)，
$H^1$ 协调性要求函数在整个网格上连续,
$H(\text{curl})$ 协调性要求函数 $u$ 在两个相邻单元的面上切向分量连续,
$H(\text{div})$ 协调性要求函数在两个相邻单元的边上法向分量连续。
常用的 $H^1$ 协调的有限元有 Lagrange 有限元、Hermite 有限元等，
$H(\text{curl})$ 协调的有限元有第一类 Nedelec 有限元、第二类 Nedelec 有限元，
$H(\text{div})$ 协调的有限元有 RT 有限元、BDM 有限元。

$H(\text{curl})$ 和 $H(\text{div})$ 协调的有限元属于向量型有限元，
其自由度相比于标量的 Lagrange 有限元更加复杂，其基函数的构造也更加困难。
传统的向量型有限元的基函数构造方式是在参考单元上构造基函数，然后
对于 $H(\text{curl})$ 协调有限元使用协变 Piola 变换， 
对于 $H(\text{div})$ 协调有限元使用逆变 Piola 变换，
将参考单元上的基函数映射到实际单元上，从而得到实际单元上的基函数。
在[1]中，基于对 $k$ 次多项式空间的几何分解，提出了一种不需要映射的基函数构造方法，
使用 Lagrange 有限元基函数乘以不同的向量，对函数施加不同的连续性，
得到 BDM 元和第二类 Nedelec 元的基函数。

### （2）弹性复形

三维空间的弹性复形：
$$
\boldsymbol{RM} \xrightarrow{\subset} H^1(\Omega, \mathbb{R}^3) 
\xrightarrow{\mathrm{sym}\nabla} H(\mathbf{inc}, \Omega, \mathbb{S})
\xrightarrow{\mathbf{inc}} H(\mathbf{div}, \Omega, \mathbb{S})
\xrightarrow{\mathbf{div}} L^2(\Omega, \mathbb{R}^3) \rightarrow 0
$$
其中 $H(\mathbf{div}, \mathbb{S})$ 协调的有限元在固体力学的数值模拟中有重要应用， 
固体力学的常用方式是位移方法，以位移为未知量，通过位移的梯度计算应力与应变，但是
对于粘弹性问题，塑性问题，位移法要么不可行，要么计算效果不好，因此将应力和位移一起
作为未知量的混合方法是一种合适的选择, 其中位移属于 $L^2$ 空间，
应力属于 $H(\text{div})$ 空间，
在1970年就提出了混合方法，但是其对应的稳定的有限元方法一直没有解决，
直至[]中，根据有限元外微分技术给出了一种稳定的有限元离散方法。
关键是构造 $H(\mathbf{div}, \mathbb{S})$ 协调的有限元，根据格林公式，
$$
\begin{aligned}
  (\nabla\cdot \boldsymbol{\sigma}, v) & = -(\boldsymbol{\sigma}, \nabla v) +
  (\boldsymbol{\sigma}\cdot \boldsymbol{n}, v)_{\partial \Omega}\\
\end{aligned}
$$
$H(\mathbf{div}, \mathbb{S})$ 协调的有限元要求 $\boldsymbol{\sigma}$ 
在两个相邻单元的边上法向分量连续，在[1]
中胡俊教授提出了以 $\mathbb{P}_k(\mathbb{S})$ 为形函数空间的 
$H(\mathbf{div}, \mathbb{S})$ 
协调的有限元: Hu-Zhang 元，其中由于 $P_k(\mathbb{S})$
的光滑性与对称性，使得有限元函数必须在顶点上施加额外的连续性。在[1]中，
提出了一种 $H(\mathbf{inc}, \mathbb{S})$ 协调的有限元，与 Hu-Zhang 元，$H^1$
协调的 Neilan 元，一起组成第一个完整的弹性复形。

### （3）$\mathrm{div}\mathbf{div}$ 复形
$\mathrm{div}\mathbf{div}$ 复形如下:
$$
RT \xrightarrow{\subset} H^1(\Omega, \mathbb{R}^3) \xrightarrow{\mathrm{dev}\nabla} 
H(\text{sym curl}, \Omega, \mathbb{T}) \xrightarrow{\mathrm{sym\ curl}} 
H(\mathrm{div}\mathbf{div}, \Omega, \mathbb{S}) \xrightarrow{\mathrm{div}\mathbf{div}} L^2(\Omega) \rightarrow 0
$$
$\mathrm{div}\mathbf{div}$ 复形是最近的研究热点，其可以用于双调和方程的混合元求解，
线性 Einstein-Bianchi 方程的求解。在双调和方程混合元方法中 
未知量的 $Hessian$ 值也被当做未知量，其属于 $H(\mathrm{div}\mathbf{div},
\mathbb{S})$ 空间，根据格林公式，
$$
(\mathrm{div}\mathbf{div}\boldsymbol{\sigma}, v) = 
-(\boldsymbol{\sigma}, \nabla^2 v) +
(\boldsymbol{\sigma}\cdot \boldsymbol{n}, \nabla v)_{\partial \Omega}
+ (\mathrm{div}\boldsymbol{\sigma}\cdot \boldsymbol{n}, v)_{\partial \Omega}
$$
[1]中提出了一种 $H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调且 $H(\mathbf{div}, \mathbb{S})$ 协调的有限元，其满足 $\boldsymbol{\sigma}\cdot \boldsymbol{n}$
以及 $\mathrm{div}\boldsymbol{\sigma}\cdot \boldsymbol{n}$ 在两个相邻单元的面上连续。
[2]中将 $\boldsymbol{\sigma}\cdot \boldsymbol{n}$ 分解为 
$\boldsymbol{t}^T\boldsymbol{\sigma}\boldsymbol{n}$ 和
$\boldsymbol{n}^T\boldsymbol{\sigma}\boldsymbol{n}$，对格林公式右端第二项进行分解，
$$
\begin{aligned}
  (\boldsymbol{\sigma}\cdot \boldsymbol{n}, \nabla v)_F & = 
  (\boldsymbol{n}^T\boldsymbol{\sigma}\boldsymbol{n}, 
  \partial_{\boldsymbol{n}} v)_F -
  (\mathbf{div}_F(\boldsymbol{\sigma n}), v)_F + 
  \sum_{e\in\partial F}(\boldsymbol{n}_{F,e}^T\boldsymbol{\sigma}\boldsymbol{n}, v)_{e}
\end{aligned}
$$
将[1]中的连续性降低为 $\boldsymbol{\sigma}$ 在面，边的法平面上连续，
$\mathbf{div}_F(\boldsymbol{\sigma n})+\boldsymbol{n}^T \mathbf{div}\boldsymbol{\sigma}$ 
在面上连续，得到了另一个 $H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 协调的有限元。
特别的，基于这个有限元，[1] 中定义了对应的 $H^1$ 协调的有限元，$H(\text{sym
curl}, \mathbb{T})$ 协调的有限元，共同组成了第一个完整的
$\mathrm{div}\mathbf{div}$ 有限元复形, 并将其应用于线性 Einstein-Bianchi 方程，
提出了线性 Einstein-Bianchi 方程对偶公式的有限元离散格式。

## 2. 研究目标
在前面提到的 $H(\mathbf{inc}, \mathbb{S})$、$H(\mathbf{div}, \mathbb{S})$、
$H(\mathrm{sym\ curl}, \mathbb{T})$、$H(\mathrm{div}\mathbf{div}, \mathbb{S})$
协调的有限元都是张量型有限元，其自由度相比于向量型的有限元还要更加复杂，
其基函数的构造也更加困难。对于向量型的有限元，尚且可以通过协变，或逆变 Piola 变换，
将参考单元上的基函数映射到实际单元上，但是对于张量型的有限元，
由于其复杂的连续性要求，这样的映射构造起来非常困难。
对于 Hu-Zhang 元, 
在[1](Nodal fem)，使用了 Lagrange
基函数乘以不同的对称向量的方式，构造了其显式基函数， 
而且形式非常简单，实现起来也相对容易，
在 [1] 中对于 
$H(\mathrm{div}\mathbf{div}, \mathbb{S})\cap H(\mathbf{div}, \mathbb{S})$ 协调的有限元，
给出了二维的构造方法。
但是其他的张量型有限元目前还没有显式的基函数构造，这对将这些新型张量有限元
和基于它们设计出来的数值离散格式应用到实际问题中是一个障碍。

另一方面，由于协调性的要求，张量有限元经常会有额外的连续性要求，如 Hu-Zhang
元要求在顶点上连续，$H(\mathrm{div}\mathbf{div}, \mathbb{S})$ 
调的有限元要求边的法平面上连续，而且以此构造的有限元复形会有更多的连续性要求，
如 [1] 中的 $\mathbf{div}$
有限元复形中，$H(\mathbf{inc}, \mathbb{S})$ 
协调的有限元在顶点上 1 阶连续, $H^1$ 协调的有限元在顶点上 2 阶连续。类似的，
[2] 中的 $\mathrm{div}\mathbf{div}$ 有限元复形中，$H(\mathrm{sym\ curl}, \mathbb{T})$
协调的有限元在顶点上 1 阶连续，$H^1$ 协调的有限元在顶点上 2 阶连续。
这些额外的连续性对算法实现来说是一个挑战，在最近研究中，基于分布 Hilbert 复形，
构造分布有限元，降低了连续性要求，如[3]中 TDNNS
方法使用法向法向连续的对称张量离散弹性方程中的应力，得到了稳定的弹性方程离散格式。
这种方法自由度简单，易于实现。

基于上述背景，本研究的目标是：对于现有的张量型有限元，如 [1，2] 中的
$\mathbf{div}$ 复形，$\mathrm{div}\mathbf{div}$ 复形中的有限元，构造其显式基函数，
并将其应用到弹性方程，双调和方程，线性 Einstein-Bianchi 方程等问题中。另外，
研究这些方程的分布有限元离散格式，分析其稳定性，收敛性，并实现其算法。

## 3. 本人博士期间的研究情况
本人在博士期间主要的研究内容是虚单元方法和有限元方法。
虚单元方法其可以看做是有限元方法的一种扩展，可以定义在任意多边形，多面体区域上，
但是由于其基函数是某个方程的解，没有显示表达式，因此其计算必须投影到多项式空间，
这样就需要引入稳定化项。

在虚单元方面，本人主要在提高虚单元光滑性，去除虚单元方法中的稳定化项，
以及降低虚单元方法对网格的依赖性等方面进行了研究。
本人及合作者提出了任意维任意次任意 $m\in \mathbb{N}$ 的 $H^m$ 协调虚单元，
其不含有超光滑的自由度，
且对多项式次数 $k$ 的要求是 $k\geq m$，小于 $H^m$
协调有限元的要求 $k\geq
(m-1)2^d+1$。在[4]中，本人及合作者基于有限元外微分，构造了多边形多面体上的宏元复形,
将虚单元函数的梯度投影到宏元空间，可以去掉虚单元方法中的稳定化项。
在[5]中，最低次的
$H(\mathrm{curl})$ 协调虚单元与多边形剖分成的三角网格上最低次棱元同构，
基于这个事实，本人及合作者分析了各向异性网格上虚单元方法的收敛性。
此外，根据虚单元方法对网格单元的要求较低的特点，本人及合作者提出了一种移动界面问题的
虚单元方法，其优点是网格生成简单，且无需对上一个时间层的解在当前时间层的网格上
进行插值。

在有限元方法方面，本人主要在构造有限元基函数进行了研究。对于向量型有限元，
基于对 $(\mathbb{P}_k)^d$ 的几何分解，
本人及合作者使用 Lagrange 有限元基函数乘以不同的向量，对函数施加不同的连续性，
提出了一种不需要映射的 BDM 元和第二类 Nedelec 元基函数的构造方法。
此外，基于[1]中对光滑元的几何分解，本人及合作者提出了任意维任意次任意 
$r\in \mathbb{N}$ 的 $C^r$ 光滑有限元基函数的构造方法，
其构造过程仅需要求解一个下三角矩阵的逆，且自由度形式与 Lagrange 元统一，
易于实现。

根据以上情况，本人基本具备了研究张量有限元的基本知识和技能，有能力完成本研究的目标，
并且本人对有限元方法的研究有浓厚的兴趣，希望能够在这个领域做出一些有意义的工作。

## 4. 看的论文
1. 有限元外微分
2. complex from complex

## 4. 问题
1. 有限元外微分为什么能够提供稳定的有限元方法？其和 InfSup 条件有什么关系？
2. 为什么 divdiv 复形有这个版本：
    $$
    RT \xrightarrow{\subset} H^1(\Omega, \mathbb{R}^3) \xrightarrow{\mathrm{dev}\nabla} 
    H(\text{sym curl}, \Omega, \mathbb{T}) \xrightarrow{\mathrm{sym\ curl}} 
    H(\mathrm{div}\mathbf{div}, \Omega, \mathbb{S}) \xrightarrow{\mathrm{div}\mathbf{div}} L^2(\Omega) \rightarrow 0
    $$
    还有这个版本： distributional divdiv complex
    $$
    RT \xrightarrow{\subset} H^1(\Omega, \mathbb{R}^3) \xrightarrow{\mathrm{dev}\nabla} 
    H(\text{sym curl}, \Omega, \mathbb{T}) \xrightarrow{\mathrm{sym\ curl}} 
    L^2(\Omega, \mathbb{S}) \xrightarrow{\mathrm{div}\mathbf{div}} H^{-2}(\Omega) \rightarrow 0
    $$
    似乎还有这个版本：
    $$
    RT \xrightarrow{\subset} H^1(\Omega, \mathbb{R}^3) \xrightarrow{\mathrm{dev}\nabla} 
    H(\text{sym curl}, \Omega, \mathbb{T}) \xrightarrow{\mathrm{sym\ curl}} 
    H(\mathrm{div}\mathbf{div}^{-1}, \Omega, \mathbb{S}) 
    \xrightarrow{\mathrm{div}\mathbf{div}} H^{-1}(\Omega) \rightarrow 0
    $$
    他们有什么关系？
3. 构造 Hilbert 复形对应的有限元复形的时候，要求他们之间要有一个上链投影算子，
    但是看很多文章都没有提到这个上链投影算子，这个是有什么理论保证其存在性么？
















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




