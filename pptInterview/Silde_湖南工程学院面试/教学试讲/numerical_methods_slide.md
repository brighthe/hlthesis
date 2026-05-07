---
marp: true
theme: numerical-methods
paginate: true
math: katex
---

<!-- _class: cover -->
<!-- _paginate: skip -->

<div class="kecheng">课程试讲：《数值分析》</div>

<div class="cover-main">

# 误差的基本概念

<div class="cover-bottom">
  <div class="info-row">
    <span>试讲人：何亮</span>
    <span>湘潭大学 · 数学与计算科学学院</span>
    <span>2026 年 4 月 29 日</span>
  </div>
</div>

---

<!-- _class: objectives -->

## 本节学习目标

- 理解数值计算误差的**来源与类型**
- 掌握**绝对误差**、**相对误差**与**有效数字**的定义与计算
- 理解误差在四则运算中的**传播规律**
- 识别数值计算中常见的**不稳定**情形，并了解基本规避思路

<div style="height: 5em;"></div>

---

<!-- _class: objectives -->
<!-- _paginate: hold -->

## 本节学习目标

- 理解数值计算误差的**来源与类型**
- 掌握**绝对误差**、**相对误差**与**有效数字**的定义与计算
- 理解误差在四则运算中的**传播规律**
- 识别数值计算中常见的**不稳定**情形，并了解基本规避思路

> **思考** 
为什么一个在纸面上看似正确的公式，到了计算机中实现时，却可能出现精度丢失，甚至得到 NaN 这样的异常结果？

---

<style scoped>
.flow { display:flex; justify-content:center; align-items:center; gap:0.6em; margin:2.5em 0; }
.nb { background:#EBF2FC; border:1.5px solid #1A4A8A; border-radius:10px; padding:0.5em 1.3em; font-weight:600; color:#1A4A8A; }
.ng { background:#E3F5F0; border:1.5px solid #0F8A6E; border-radius:10px; padding:0.5em 1.3em; font-weight:600; color:#0F8A6E; }
.aw { display:flex; flex-direction:column; align-items:center; gap:0.15em; }
.arr { color:#8A94A6; font-size:1.3em; line-height:1; letter-spacing:-0.1em; }
.lbl { font-size:0.55em; color:#8A94A6; white-space:nowrap; line-height:1.4; }
</style>

## 一、误差的来源与类型

<div class="flow">
  <span class="nb">实际问题</span>
  <span class="aw">
    <span class="lbl">模型误差</span>
    <span class="arr">→</span>
    <span class="lbl">观测误差</span>
  </span>
  <span class="nb">数学模型</span>
  <span class="aw">
    <span class="lbl">截断误差</span>
    <span class="arr">→</span>
    <span class="lbl">&nbsp;</span>
  </span>
  <span class="nb">数值方法</span>
  <span class="aw">
    <span class="lbl">舍入误差</span>
    <span class="arr">→</span>
    <span class="lbl">传播误差</span>
  </span>
  <span class="ng">计算结果</span>
</div>

---

<!-- _paginate: hold -->

<style scoped>
.flow { display:flex; justify-content:center; align-items:center; gap:0.6em; margin:2.5em 0; }
.nb { background:#EBF2FC; border:1.5px solid #1A4A8A; border-radius:10px; padding:0.5em 1.3em; font-weight:600; color:#1A4A8A; }
.ng { background:#E3F5F0; border:1.5px solid #0F8A6E; border-radius:10px; padding:0.5em 1.3em; font-weight:600; color:#0F8A6E; }
.aw { display:flex; flex-direction:column; align-items:center; gap:0.15em; }
.arr { color:#8A94A6; font-size:1.3em; line-height:1; letter-spacing:-0.1em; }
.lbl { font-size:0.55em; color:#8A94A6; white-space:nowrap; line-height:1.4; }
</style>

## 一、误差的来源与类型

<div class="flow">
  <span class="nb">实际问题</span>
  <span class="aw">
    <span class="lbl">模型误差</span>
    <span class="arr">→</span>
    <span class="lbl">观测误差</span>
  </span>
  <span class="nb">数学模型</span>
  <span class="aw">
    <span class="lbl">截断误差</span>
    <span class="arr">→</span>
    <span class="lbl">&nbsp;</span>
  </span>
  <span class="nb">数值方法</span>
  <span class="aw">
    <span class="lbl">舍入误差</span>
    <span class="arr">→</span>
    <span class="lbl">传播误差</span>
  </span>
  <span class="ng">计算结果</span>
</div>

| 误差类型 | 来源 |
|:------:|:---:|
| 模型误差 | 数学建模时的简化假设 |
| 观测误差 | 测量仪器精度限制 |
| **截断误差** | 无穷级数取有限项 |
| **舍入误差** | 浮点数有限位表示 |
| **传播误差** | 已有误差随运算过程积累 |

---

## 模型误差与观测误差

在五类误差中，模型误差和观测误差更多来自现实问题本身，通常不是由数值方法内部产生的误差。

<div class="cols">
<div class="col-left">

1. **模型误差（建模层面）**
由数学建模中的理想化和简化假设引起的偏差。例如：忽略空气阻力、将复杂材料简化为理想材料、人口预测中忽略突发事件影响。

2. **观测误差（数据层面）**
由测量仪器、环境波动和人工读数等因素引起的偏差。例如：长度测量中的刻度误差、温度测量中的波动误差、传感器数据中的随机扰动。

因此，后续将重点讨论数值计算过程中更需要控制的三类误差：截断误差、舍入误差和传播误差。

---

## 数值计算中重点关注的三类误差

下面从数学过程、有限精度表示和算法传播三个层面，介绍数值计算中的三类核心误差：

<div class="cols">
<div class="col-left">

1. **截断误差 (数学层面)**
用有限过程逼近无限过程时产生的误差。例如：Taylor 展开取有限项、连续问题离散化。

2. **舍入误差 (硬件层面)**
计算机只能进行有限精度表示与运算，因此实数进入机器后会产生舍入误差。

3. **传播误差 (算法层面)**
已有误差会在后续运算中继续传递，并可能累积或放大。当出现相近数相减、近零数作除数时，这种现象尤其明显。

---

<!-- _class: errors -->

## 三类误差的典型例子

<div class="intro">下面通过简单例子，进一步理解<strong>截断误差</strong>、<strong>舍入误差</strong>和<strong>传播误差</strong>。</div>

<div class="three-cols">

<div class="err-card">

### 截断误差 </span>

<div class="desc">用有限项近似无限过程时产生误差。</div>

<div class="formula">

$$e^{-x} \approx 1 - x + \frac{x^2}{2} - \frac{x^3}{6}$$

</div>

<div class="conclude">只保留 Taylor 展开的前几项，因此会产生<strong>截断误差</strong>。</div>

</div>

<div class="err-card">

### 舍入误差</span>

<div class="desc">有些实数不能被计算机精确表示，只能存储为近似值。</div>

<div class="formula">

$$\frac{1}{3} \approx 0.33333$$

$$\sqrt{2} \approx 1.41421,\ \pi \approx 3.14159$$

</div>

<div class="conclude">有限精度表示与舍入处理带来的误差就是<strong>舍入误差</strong>。</div>

</div>

<div class="err-card">

### 传播误差</span>

<div class="desc">近似值继续参与运算时，误差会继续传递。</div>

<div class="formula">

$$\frac{1}{3} \approx 0.33333$$

$$0.33333 \times 3 = 0.99999 \neq 1$$

</div>

<div class="conclude">已有误差随计算过程继续传递，这就是<strong>传播误差</strong>。</div>

</div>

</div>

> **结论**：数值计算不仅要关注“如何近似”，还要关注“误差如何产生与传播”。

---

## 二、近似数的误差和有效数字

<h3 style="margin-bottom:0.1em;">绝对误差</h3>

$$E(x) = x - x^*$$

若 $|x - x^*| \leq \eta$，则称 $\eta$ 为近似数 $x^*$ 与准确数 $x$ 的**绝对误差限**。

<!-- **含义：** 绝对误差反映近似值与准确值相差多少。 -->

</div>

<h3 style="margin-bottom:0.1em;">相对误差</h3>

$$E_r(x) = \frac{x - x^*}{x} = \frac{E(x)}{x},\quad x\neq0$$

若存在正数 $\delta$，使得 $|E_r(x)| \leq \delta$，则称 $\delta$ 为近似数 $x^*$ 与准确数 $x$ 的**相对误差限**。

<!-- **含义：** 相对误差反映误差相对于准确值的大小。 -->

> **思考**：绝对误差相同，是否意味着近似程度也相同？

---

<!-- _class: two-col -->

## 思考：绝对误差相同，近似程度一定相同吗？

<div class="intro">考察两个绝对误差均为 <strong>0.1</strong> 的近似——它们的“近似质量”是否一致？</div>

<div class="cols">
<div class="col-left">

### 例 1：$x = 1000,\ x^* = 999.9$

**绝对误差**

$$x - x^* = 0.1$$

**相对误差**

$$\frac{x - x^*}{x} = \frac{0.1}{1000} = 10^{-4}$$

</div>
<div class="col-right">

### 例 2：$x = 1,\ x^* = 0.9$

**绝对误差**

$$x - x^* = 0.1$$

**相对误差**

$$\frac{x - x^*}{x} = \frac{0.1}{1} = 10^{-1}$$

</div>
</div>

> **结论**：绝对误差相同，**近似质量未必相同**；判断近似程度还需结合准确数本身的**数量级**——引入**相对误差**的动机。

---

## 有效数字

设精确值 $x$ 的近似值为

$$x^* = \pm 0.a_1 a_2 \cdots a_n \times 10^m$$

其中 $a_i(i=1,2,\cdots,n)$ 是 0 到 9 中的某一整数，并且 $a_1\neq 0$。

- 若 $|x - x^*| \leq \frac{1}{2} \times 10^{m-k}$，则称 $x^*$ 具有 $k$ 位**有效数字**。

- 若 $x^{*}$ 的每一位都是有效数字，那么 $x^{*}$ 称为具有 $n$ 位有效数字的**有效数**。


**例：** 用舍入法取 $e = 2.7182818\cdots$ 的近似数 $2.7183$

$$|e - 2.7183| = 0.0000181\cdots \leq 0.5 \times 10^{-4}$$

所以 $2.7183$ 有 **5 位有效数字**。

---

## 三、初始误差在运算中的传播

设 $y = f(x_1, x_2, \ldots, x_n)$，当各输入量存在微小误差时，由 Taylor 一阶展开可得：

$$E(y) \approx \sum_{i=1}^{n} \frac{\partial f}{\partial x_i} E(x_i)$$

$$E_r(y) \approx \sum_{i=1}^{n} \frac{\partial f}{\partial x_i} \cdot \frac{x_i}{f} \cdot E_r(x_i),\quad f(x_1, x_2, \ldots, x_n)\neq0$$

**公式告诉我们什么？**

- 绝对误差是一阶近似下各输入误差**加权累积**的结果
- 相对误差不仅与输入误差有关，还与**函数形式**有关
- 当结果接近 0 或出现**相近数相减**时，相对误差可能显著放大

---

<!-- _class: two-col -->

## 加减法的误差传播

当两个相近数作减法时，计算结果的精度可能显著下降。

<div class="cols">
<div class="col-left">

### ❌ 直接计算（不稳定）

用四位浮点数计算：

$$1 - \cos 2° \approx 1 - 0.9994 = 0.0006$$

**仅有 1 位有效数字！**

> **问题根源**：
相近数相减会造成有效数字严重丢失。

</div>
<div class="col-right">

### ✅ 等价变换（稳定）

利用恒等式 $1 - \cos x = 2\sin^2\!\frac{x}{2}$：

$$1 - \cos 2° \approx 0.6092 \times 10^{-3}$$

**恢复 4 位有效数字！**

> **原则：**
通过等价变换避免相近数相减。

---

<!-- _class: two-col -->

## 乘除法的误差传播

在除法运算中，若分母很小，输入中的微小误差可能被显著放大。

<div class="cols">
<div class="col-left">

### ❌ 直接计算（不稳定）

用四位浮点数计算：

$$\frac{\sin 1°}{1 - \cos 1°} \approx \frac{0.01745}{1 - 0.9998} = \frac{0.01745}{0.0002} = 87.25$$

**结果偏差很大！**

> **问题根源：**
分母很小，且由两个相近数相减得到，误差在除法中被进一步放大。

</div>
<div class="col-right">

### ✅ 等价变换（稳定）

利用恒等式 $\dfrac{\sin x}{1 - \cos x} = \dfrac{1 + \cos x}{\sin x},$

$$\frac{\sin 1°}{1 - \cos 1°} \approx \frac{1 + 0.9998}{0.01745} = \frac{1.9998}{0.01745} \approx 114.6$$

**结果更接近真实值！**

> **原则：**
通过等价变换，避免使用过小且不稳定的分母。

---

<!-- _class: two-col -->

<style scoped>
.nb {
  background: var(--brand-light);
  border: 1.5px solid var(--brand);
  border-radius: var(--radius);
  padding: .6em 1.2em;
  font-weight: 600;
  color: var(--brand);
  font-size: .92em;
}
</style>

## 四、浮点运算与机器数的局限性

计算机的存储空间是有限的，底层数据由**规格化浮点机器数**表示，记为集合 $\mathbb{F}$。
当 $x, y \in \mathbb{F}$ 时，它们的精确实数运算结果 $x \pm y,\ x \times y,\ x / y$ **未必仍属于** $\mathbb{F}$，因此计算机必须对结果进行舍入，返回一个新的机器数。

**例如**：执行加法 $x + y$ 时，先将小阶数向大阶数对齐：

$$x + y = 0.003127 \times 10^{-4} + 0.4153 \times 10^{-4}$$

$$= 0.418427 \times 10^{-4}$$

尾数膨胀到 6 位，若机器尾数只保留 4 位，则该结果超出机器表示容量，

必须截断舍入为 $0.4184 \times 10^{-4}$，故 $x + y \notin \mathbb{F}$。

> **核心数学事实**：机器数集合 $\mathbb{F}$ 对四则运算一般不封闭，因此计算机运算通常伴随着舍入。

---

<!-- _class: two-col -->

<style scoped>
  section > p   { margin: 0.45em 0 !important; }
  .cols         { margin: 0.3em 0 !important; }
  .cols .katex-display { margin: 1.0em 0 !important; }
  .cols h3      { margin: 0 0 1.5em !important; }
</style>

## 浮点运算中的非结合现象

实数加法满足结合律 $(a+b)+c = a+(b+c)$。但在**有限精度浮点运算**下，每一步运算后都可能发生舍入，该等式**未必严格成立**。

**例如** 取 $a = 10^{16},\ b = -10^{16},\ c = 1$，比较两种运算顺序：

<div class="cols">
<div class="col-left">

### 顺序一：$(a+b)+c$

$$a + b = 10^{16} - 10^{16} = 0$$

$$(a+b) + c = 0 + 1 = 1$$

</div>
<div class="col-right">

### 顺序二：$a+(b+c)$

$$b + c = -10^{16} + 1 \;\xrightarrow{\text{舍入}}\; -10^{16}$$

$$a + (b+c) = 10^{16} - 10^{16} = 0$$

</div>
</div>

**原因**：$b + c$ 中的 $1$ 相对于 $10^{16}$ 量级太小，在尾数对齐时被舍入吞掉——即典型的**大数吃小数**现象。

> **结论**：浮点运算一般不满足严格的**结合律**，**运算顺序可能直接影响最终结果**。

---

## 五、数值案例

**考察函数：**

$$f(x, y) = \frac{|x| - |y|}{x - y}$$

接下来，我们通过下面三个问题来分析这个函数：

- 当 $x$ 与 $y$ 很接近时，若直接用**有限精度浮点运算**计算，会出现什么问题？

- 能否通过**等价变换**改写这个函数，从而避免不稳定的减法结构？

- 两个在数学上恒等的公式，在**计算机有限精度运算**下，结果会有什么差别？

---

## 原始函数的不稳定性分析

**原始函数：** $f(x, y) = \dfrac{|x| - |y|}{x - y}$

基于多元函数误差传播公式，可将其相对误差写为

$$\footnotesize  E_r(f) = \left[\frac{\text{sign}(x)}{|x|-|y|} - \frac{1}{x-y}\right]xE_r(x) + \left[\frac{-\text{sign}(y)}{|x|-|y|} + \frac{1}{x-y}\right]yE_r(y)$$

下面主要考虑 $x \approx y > 0$ 的情形，将分母 $x-y$ 代入误差传播公式，可得

$$\small E_r(x-y) \approx \frac{x}{x-y}E_r(x) - \frac{y}{x-y}E_r(y),\quad \left|\frac{x}{x-y}\right|, \left|\frac{y}{x-y}\right|\to\infty$$

> **结论**： 原始函数同时包含“相近数相减”和“分母过小”，直接计算容易产生严重的误差放大。

---

## 重构函数：化减法为加法

对原函数分子分母同乘 $(|x| + |y|)$，并约去公因子 $(x - y)$，可得重构函数 $g(x,y)=\dfrac{x+y}{|x|+|y|}$

基于多元函数误差传播公式，可将其相对误差写为

$$\footnotesize E_r(g) = \left[\frac{-\text{sign}(x)}{(|x| + |y|)} + \frac{1}{(x + y)}\right]xE_r(x) + \left[\frac{-\text{sign}(y)}{(|x| + |y|)} + \frac{1}{(x + y)}\right]yE_r(y)$$

减法被成功转化为同号相加。当 $x \approx y > 0$ 时，对分母 $x+y$ 有

$$\small E_r(x+y) \approx \frac{x}{x+y}E_r(x) + \frac{y}{x+y}E_r(y)\quad \frac{x}{x+y}, \frac{y}{x+y}\approx0.5$$

> **结论**： 等价变换将减法改写为加法，使误差系数保持有界，计算显著更稳定。

---

<style scoped>
.fig-row { display: grid; grid-template-columns: 1fr 1fr; gap: 1em; margin: 0.3em 0; }
.fig-box { display: flex; flex-direction: column; align-items: center; gap: 0; }
.fig-note-warn, .fig-note-ok { font-size: .80em; margin-top: -0.6em; }
.fig-note-warn { color: var(--warn); }
.fig-note-ok   { color: var(--accent); }
blockquote { margin-top: 0.2em; }
</style>

## 可视化验证：数学上的等价 ≠ 计算机底层的等价

<div class="fig-row">
<div class="fig-box">

![h:300](Figure_1.png)

<span class="fig-note-warn">⚠️ $f(x,y)$：在 $x \approx y$ 附近易出现数值震荡，甚至产生 NaN。</span>
</div>
<div class="fig-box">

![h:300](Figure_2.png)

<span class="fig-note-ok">✅ $g(x,y)$：化减为加，显著改善计算稳定性。</span>
</div>
</div>

> **注意**：在具体的数值实现中，为避免 $|x|+|y|$ 过小或为 0，通常在分母中加入极小常数 $\epsilon$ 作为数值保护。

---

<!-- _class: summary -->

## 课堂小结：从误差度量到稳定计算

<div class="intro">本节围绕“<strong>误差从哪里来、如何度量、怎样传播、如何规避</strong>”四个问题展开。</div>

<div class="summary-grid">

<div class="card">

<span class="num">01</span>

### 误差来源：误差不可避免

模型误差、观测误差、截断误差、舍入误差与传播误差共同影响计算结果；数值分析重点关注**截断误差**、**舍入误差**和**传播误差**。

</div>

<div class="card">

<span class="num">02</span>

### 误差度量：近似需要标准

**绝对误差**刻画偏差大小，**相对误差**刻画相对精度，**有效数字**反映近似数中可靠数字的位数。

</div>

<div class="card">

<span class="num">03</span>

### 误差传播：运算可能放大误差

初始误差会随计算过程继续传递；**相近数相减**、**近零数作除数**时，误差可能被显著放大。

</div>

<div class="card">

<span class="num">04</span>

### 稳定计算：数学等价 ≠ 数值等价

有限精度浮点运算下，数学上等价的公式可能表现出不同的**数值稳定性**；实际计算中应尽量选择更稳定的等价表达式。

</div>

</div>

> 数值分析的核心任务之一，就是在<strong>有限精度计算环境下</strong>，构造既<strong>准确</strong>又<strong>稳定</strong>的计算方法。

---

<!-- _class: cover -->
<!-- _paginate: skip -->

<div class="kecheng">课程试讲：《数值分析》</div>

<div class="cover-main">

# 感谢聆听

<div class="subtitle">试讲到此结束，欢迎各位老师批评指正</div>



