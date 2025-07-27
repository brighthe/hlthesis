# 《二维 H(curl) 协调虚单元》编程数学文档
## 一、引言


## 二、问题描述

## 三、算法数学原理
### （1）基本知识
1. 三个恰当复形
$$
\mathbb{R} \xrightarrow{i} \mathbb{P}_r \xrightarrow{\nabla
}\left(\mathbb{P}_{r-1}\right)^2 \xrightarrow{\text { rot }} \mathbb{P}_{r-2}
\xrightarrow{o} 0
$$

$$
\mathbb{R} \stackrel{i}{\longrightarrow} \mathbb{P}_r \xrightarrow{\mathbf{ rot
}}\left(\mathbb{P}_{r-1}\right)^2 \xrightarrow{\text { div }} \mathbb{P}_{r-2}
\xrightarrow{o} 0
$$

$$
\mathbb{R} \xrightarrow{i} \mathbb{P}_r
\xrightarrow{\nabla}\left(\mathbb{P}_{r-1}\right)^3 \xrightarrow{\mathbf{ curl
}}\left(\mathbb{P}_{r-2}\right)^3 \xrightarrow{\text { div }} \mathbb{P}_{r-3}
\xrightarrow{o} 0
$$

2. 几个空间的定义
- $K$ 是一个 $d$ 维单纯形, 记 $\mathbb{P}_k(K)$ 的维数为 $\pi_{k, d}$
- $\mathcal{G}_{s-1}(K) = \nabla \mathbb{P}_s(K)$
- $K$ 是三角形时：$\mathcal{R}_{s-1}(K) = \mathbf{rot}\mathbb{P}_s(K)$
- $K$ 是四面体时：$\mathcal{R}_{s-1}(K) = \mathbf{curl}(\mathbb{P}_s(K))^3$
- $\mathcal{R}_s^{\perp}(K)$ 和 $\mathcal{G}_s^{\perp}(K)$ 分别是 $\mathcal{R}_s(K)$ 和 $\mathcal{G}_s(K)$ 在 $\mathbb{P}_k(K)$ 中的正交补空间。

3. 映射与空间的关系
- $\nabla$ 是 $\mathbb{P}_s(K)/\mathbb{R}$ 到 $\mathcal{G}_{s-1}(K)$ 的同构映射
- $\mathcal{R}_s(K) = \mathrm{ker}(\mathrm{div})\cap \mathbb{P}_s(K)$
- $\mathrm{div}$ 是 $\mathcal{R}_{s}^{\perp}(K)$ 到 $\mathbb{P}_{s-1}(K)$ 的同构映射

### （2）二维 $H(\mathrm{rot})$ 协调虚单元
1. $H(\mathrm{curl})$ 协调虚单元的定义
    $$
    \begin{array}{r}
    V_{2, k}^{\text {face }}(K):=\left\{\boldsymbol{v} \in H(\operatorname{div} ; K)
    \cap H(\operatorname{rot} ; K): \boldsymbol{v} \cdot \boldsymbol{n}_{\mid e} \in
    \mathbb{P}_k(e)\ \ \forall \text { edge } e \text { of } K,\right. \\ \text { grad
    } \left.\operatorname{div} \boldsymbol{v} \in \mathcal{G}_{k-2}(K), \text { and
    } \operatorname{rot} \boldsymbol{v} \in \mathbb{P}_{k-1}(K)\right\} .
    \end{array}
    $$
    **计算一个空间的维数的方法: 看有多少个条件可以唯一确定一个元素。**

2. 自由度定义
    - 边界自由度
        $$
        \int_{e} \boldsymbol{v} \cdot \boldsymbol{n}\ q \mathrm{d} s, \quad \forall
        q \in \mathbb{P}_{k}(e)
        $$
    - 内部自由度
        $$
        \int_{K} \boldsymbol{v} \cdot \boldsymbol{q} \mathrm{d} x, \quad \forall
        \boldsymbol{q} \in \mathcal{G}_{k-2}(K)
        $$
        $$
        \int_{K} \boldsymbol{v} \cdot \boldsymbol{q} \mathrm{d} x, \quad \forall
        \boldsymbol{q} \in \mathbb{G}_{k}^{\perp}(K)
        $$


## 四  、程序设计

## 五、测试数据与结果

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
