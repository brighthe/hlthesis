大连理工大学郭旭院士团队（[郭旭院士个人主页](https://faculty.dlut.edu.cn/2000011087/ "null")）长期依托工业装备结构分析国家重点实验室，深耕计算力学、结构优化以及工业级科学计算软件开发。团队打破了传统基于像素/体素（Pixel/Voxel）的隐式拓扑优化范式，主导开发了**显式拓扑优化理论（MMC/MMV）**，并在**工业级软件工程化（SIPESC）** 与 **问题无关机器学习（PIML）** 赋能的超大规模结构优化方面取得了具有国际影响力的突破。

# 研究方向一：显式拓扑优化 (MMC/MMV 框架)

《Doing Topology Optimization Explicitly and Geometrically—A New Moving Morphable Components Based Framework》
《A new topology optimization approach based on moving morphable components (MMC) and the ersatz material model》
《Explicit structural topology optimization based on moving morphable components (MMC) with curved skeletons》
《Explicit three dimensional topology optimization via Moving Morphable Void (MMV) approach》
《A moving morphable void (MMV)-based explicit approach for topology optimization considering stress constraints》
《Optimal design of shell-graded-infill structures by a hybrid MMC-MMV approach》
《A comprehensive review of explicit topology optimization based on moving morphable components (MMC) method》


## 核心文献与发展脉络

- **奠基与开山之作（确立显式几何描述与等效材料模型范式）：**
    
	1. _Doing Topology Optimization Explicitly and Geometrically—A New Moving Morphable Components Based Framework_
	
    2. _A new topology optimization approach based on moving morphable components (MMC) and the ersatz material model_
    
- **形态与维度的拓展 (曲线组件与三维空间孔洞演化)：** 
	
	3. _Explicit structural topology optimization based on moving morphable components (MMC) with curved skeletons_ 
	
	4. _Explicit three dimensional topology optimization via Moving Morphable Void (MMV) approach_
    
- **物理约束攻坚与复杂结构应用 (解决局部奇异性与多孔填充)：** 
	
	5. _A moving morphable void (MMV)-based explicit approach for topology optimization considering stress constraints_ 
	
	6. _Optimal design of shell-graded-infill structures by a hybrid MMC-MMV approach_
    
- **权威总结：** 
	7. _A comprehensive review of explicit topology optimization based on moving morphable components (MMC) method_

## 核心研究方向二：问题无关机器学习 (PIML)

问题无关机器学习（Problem-Independent Machine Learning, PIML）是郭旭院士团队近年来在“人工智能赋能计算力学与结构优化”方向形成的标志性方法体系。其核心并不是用神经网络直接端到端预测最终拓扑构型，而是面向拓扑优化中最耗时的有限元分析环节，学习可局部复用、可嵌入有限元框架的分析算子，从而在保持力学一致性的前提下显著降低大规模结构分析与设计优化的计算成本。

在传统密度法或显式拓扑优化中，每一次迭代都需要在当前材料分布下求解离散平衡方程

$$
\mathbf{K}(\rho)\mathbf{u}=\mathbf{f},
$$

其中 $\rho$ 表示设计变量或材料分布，$\mathbf{K}(\rho)$ 为依赖材料分布的整体刚度矩阵，$\mathbf{u}$ 为位移向量，$\mathbf{f}$ 为外载荷向量。对于高分辨率拓扑优化问题，有限元自由度和设计变量规模可同时达到千万级甚至十亿级，结构响应分析往往成为整个优化流程的主要瓶颈。PIML 的关键思想是将机器学习对象从“具体优化问题的最终解”转移到“局部有限元分析对象”，例如扩展多尺度有限元（EMsFEM）中的多尺度形函数、子结构数值形函数或凝聚刚度矩阵。典型映射可抽象为

$$
\mathcal{G}_{\theta}: \rho_s(\mathbf{x}) \longmapsto \{N_i^{\mathrm{ms}}(\mathbf{x})\}_{i=1}^{n_s},
$$

或在子结构框架下写成

$$
\mathcal{G}_{\theta}: \big(\rho_s,\, \mathcal{T}_s\big) \longmapsto \left(\boldsymbol{\Phi}_s,\, \widetilde{\mathbf{K}}_s\right),
$$

其中 $\rho_s$ 表示粗单元或子结构内部的材料分布，$\mathcal{T}_s$ 表示局部单元或子结构的几何信息，$\boldsymbol{\Phi}_s$ 表示局部数值形函数或降阶基，$\widetilde{\mathbf{K}}_s$ 表示凝聚后的子结构刚度矩阵。由于这些局部对象主要由控制方程、材料分布和局部几何决定，而不直接依赖具体边界条件、载荷形式或整体设计域，因而具备“问题无关”的可复用特征。

### 方法本质与技术内核

PIML 的方法论可以概括为“局部算子学习 + 多尺度有限元/子结构分析 + 在线拓扑优化”。其关键特征包括：

1. **学习对象局部化。** 传统机器学习拓扑优化方法多试图学习“边界条件/载荷/体积分数 $\rightarrow$ 最优拓扑”的全局映射，泛化能力高度依赖训练样本分布。PIML 则学习粗单元或子结构内部的局部数值形函数、凝聚刚度或等效分析算子，使训练对象从问题级别下降到单元/子结构级别。

2. **有限元框架内嵌。** PIML 并不替代力学控制方程，而是将学习得到的局部分析对象嵌入 EMsFEM 或子结构有限元框架中，再通过整体装配与平衡方程求解获得结构响应。因此，它本质上是一种机器学习增强的有限元分析加速方法，而不是脱离物理方程的黑箱预测器。

3. **离线训练与在线复用分离。** 早期 PIML 通过离线样本训练建立“材料分布 $\rightarrow$ 多尺度形函数/凝聚刚度”的映射；在线优化阶段只需进行快速神经网络推断和整体方程求解，从而减少每轮有限元分析的代价。

4. **从监督数据驱动走向 data-free 力学约束学习。** 2024 年 JMPS 论文进一步提出 mechanics-based data-free PIML，将多尺度形函数视为坐标的连续函数，并通过最小势能原理构造力学约束损失函数，避免了大规模预生成训练数据集的成本。这一工作是 PIML 从“监督学习加速器”走向“力学原理约束算子学习”的关键转折。

5. **面向超大规模与高性能并行。** PIML 的局部子结构推断天然适合并行化，可与矩阵自由实现、多重网格预条件、GPU/异构计算和工业级有限元软件平台结合，为十亿级设计变量拓扑优化提供可行路径。

### 核心文献与发展脉络

- **前导工作：显式拓扑优化与机器学习结合。**  
  _Machine learning-driven real-time topology optimization under moving morphable component-based framework_ 以 MMC 显式几何表示为基础，将机器学习用于实时拓扑优化预测。该工作尚未形成严格意义上的 PIML，但确立了“显式几何参数化 + 机器学习加速结构优化”的早期方向。

- **过渡工作：有限元分辨率映射与分析加速。**  
  _A mechanistic-based data-driven approach to accelerate structural topology optimization through finite element convolutional neural network (FE-CNN)_ 将机器学习用于高、低分辨率有限元系统之间的映射，是从“直接预测拓扑构型”转向“加速有限元分析过程”的重要过渡。

- **PIML 正式提出：EMsFEM 多尺度形函数学习。**  
  _Problem-independent machine learning (PIML)-based topology optimization—a universal approach_ 正式提出 PIML 概念。该文将 EMsFEM 中粗单元内部的多尺度形函数作为学习对象，建立局部材料分布到多尺度形函数的隐式映射，从而使训练模型能够在相同控制方程类型下跨边界条件、载荷和设计域复用。这是 PIML 方向的开山之作。

- **三维与子结构扩展：大规模线弹性结构分析。**  
  _A Problem-Independent Machine Learning (PIML) enhanced substructure-based approach for large-scale structural analysis and topology optimization of linear elastic structures_ 将 PIML 从 EMsFEM 进一步推广到子结构有限元框架中，强调对凝聚刚度矩阵和局部数值形函数的高效预测，并展示其在三维大规模结构分析和拓扑优化中的应用潜力。

- **复杂设计域扩展：等参单元 PIML。**  
  _Problem-independent machine learning-enhanced structural topology optimization of complex design domains based on isoparametric elements_ 针对规则网格 PIML 对复杂设计域适应性不足的问题，引入等参单元描述局部几何，将“单元形状 + 子结构内部材料分布”共同作为学习输入，从而使 PIML 从规则网格扩展到更一般的复杂设计域。

- **方法论升级：mechanics-based data-free PIML。**  
  _A mechanics-based data-free problem independent machine learning (PIML) model for large-scale structural analysis and design optimization_ 以力学约束损失函数替代预制监督数据集，将多尺度形函数表示为连续坐标函数，并通过算子学习框架实现 data-free 训练。该工作解决了早期 PIML 仍需大量离线训练样本的问题，是 PIML 方向最重要的方法论升级之一。

- **显式微结构应用：PIML + MMC 三维晶格复合结构。**  
  _Problem-Independent Machine Learning (PIML) enhanced 3D lattice composite structures optimization via moving morphable components approach_ 将 PIML 与 MMC 显式几何建模相结合，服务于三维梯度晶格复合结构优化。该工作说明 PIML 不仅可用于连续体高分辨率分析加速，也可以作为复杂显式微结构设计中的快速分析内核。

- **高性能并行扩展：PIML + HPC。**  
  _A high-performance parallel algorithm based on problem independent machine learning (PIML) for large-scale topology optimization_ 将 PIML 与并行计算、多重网格和矩阵自由实现相结合，进一步面向大规模拓扑优化的高性能求解。该方向标志着 PIML 从方法验证进入面向工程规模问题的并行算法阶段。

### 代表性论文清单

| 年份 | 类型 | 论文题目 | 期刊/载体 | DOI/编号 | 主要贡献 |
|---|---|---|---|---|---|
| 2019 | 前导工作 | _Machine learning-driven real-time topology optimization under moving morphable component-based framework_ | _Journal of Applied Mechanics_ | 10.1115/1.4041319 | MMC 显式拓扑优化与机器学习结合，服务于实时拓扑预测。 |
| 2021 | 过渡工作 | _A mechanistic-based data-driven approach to accelerate structural topology optimization through finite element convolutional neural network (FE-CNN)_ | arXiv / CoRR | 10.48550/arXiv.2106.13652 | 从直接预测拓扑转向有限元分辨率映射与分析加速。 |
| 2022 | PIML 创始 | _Problem-independent machine learning (PIML)-based topology optimization—a universal approach_ | _Extreme Mechanics Letters_ | 10.1016/j.eml.2022.101887 | 正式提出 PIML，学习 EMsFEM 多尺度形函数，实现问题无关复用。 |
| 2023 | 大规模扩展 | _A Problem-Independent Machine Learning (PIML) enhanced substructure-based approach for large-scale structural analysis and topology optimization of linear elastic structures_ | _Extreme Mechanics Letters_ | 10.1016/j.eml.2023.102041 | 将 PIML 推广到子结构框架，面向三维大规模线弹性分析与拓扑优化。 |
| 2024 | 复杂设计域 | _Problem-independent machine learning-enhanced structural topology optimization of complex design domains based on isoparametric elements_ | _Extreme Mechanics Letters_ | 10.1016/j.eml.2024.102237 | 引入等参单元，使 PIML 适应复杂几何与非规则设计域。 |
| 2024 | 方法论升级 | _A mechanics-based data-free problem independent machine learning (PIML) model for large-scale structural analysis and design optimization_ | _Journal of the Mechanics and Physics of Solids_ | 10.1016/j.jmps.2024.105893 | 提出 data-free PIML，通过力学约束损失函数学习连续多尺度形函数。 |
| 2025 | 显式微结构应用 | _Problem-Independent Machine Learning (PIML) enhanced 3D lattice composite structures optimization via moving morphable components approach_ | _Composite Structures_ | 10.1016/j.compstruct.2025.119330 | 将 PIML 与 MMC 结合，用于三维梯度晶格复合结构优化。 |
| 2026 | 并行高性能扩展 | _A high-performance parallel algorithm based on problem independent machine learning (PIML) for large-scale topology optimization_ | _Acta Mechanica Sinica_ | 10.1007/s10409-025-25942-x | 将 PIML 与并行计算、多重网格、矩阵自由实现结合，面向大规模并行拓扑优化。 |

### 与 MMC/MMV、SIPESC 和 SOPTX 的潜在结合点

从郭旭院士团队整体研究体系看，PIML 可以被理解为连接“显式拓扑优化理论”和“工业级高性能计算软件”的关键分析加速层。其与后续研究方向的结合主要体现在三个方面。

首先，PIML 与 MMC/MMV 具有天然互补性。MMC/MMV 通过显式组件或显式孔洞降低设计变量维度，并提供清晰的几何边界；PIML 则可以作为复杂显式结构的快速有限元分析内核，降低每次设计迭代中的响应分析成本。特别是在三维晶格、壳-填充结构和多尺度显式微结构设计中，MMC/MMV 负责几何表达，PIML 负责高效分析，两者可以形成“显式几何建模 + 快速力学分析”的统一框架。

其次，PIML 与 SIPESC 等国产工业软件平台的结合点在于工程级分析内核的加速。PIML 学习的是可装配、可复用的局部分析对象，因此更容易嵌入现有有限元软件体系，而不是额外构造一个完全独立的黑箱优化器。对于工业软件而言，这种“保留有限元装配与求解框架，只替换或加速局部分析算子”的方式更符合工程软件的稳定性、可验证性和模块化集成要求。

最后，PIML 与 SOPTX/FEALPy 的多后端架构高度契合。PIML 的离线训练、在线推断、局部子结构装配和灵敏度计算可分别映射到 PyTorch/JAX/NumPy 等不同后端；同时，子结构级推断天然适合 GPU 并行和分布式并行。对于 SOPTX 后续发展，可以将 PIML 设计为独立的 `analysis_accelerator` 或 `surrogate_fem` 模块，使其与材料插值、过滤器、优化器和有限元求解器解耦，从而形成面向超大规模拓扑优化的可扩展分析加速接口。

### 可重点关注的后续研究问题

1. **PIML 与高阶有限元/多分辨率拓扑优化结合。** 可以探索将 PIML 学习对象从低阶子结构形函数扩展到高阶有限元局部基函数或多分辨率密度单元内部的等效分析算子，以服务于高精度 MTOP 框架。

2. **PIML 与任意次胡张混合元结合。** 对近不可压缩材料和应力约束问题，可以研究面向应力-位移混合变量的 PIML 局部算子学习，使其不仅加速位移场分析，也能保持应力场的高精度与物理一致性。

3. **Data-free PIML 的误差估计与可信计算。** 现有 PIML 更强调效率提升和数值验证，后续若面向工业软件集成，需要进一步建立残差驱动误差估计、局部失效检测和自适应再训练机制，形成可认证的 PIML 分析框架。

4. **面向 SIPESC 的并行 PIML 软件实现。** 可以围绕矩阵自由有限元、多重网格预条件、局部神经算子推断和 GPU/CPU 异构并行，构建适合工业级 CAE 平台调用的 PIML 高性能分析模块。

总体而言，郭旭院士团队的 PIML 研究路线具有清晰的演进逻辑：从 MMC + ML 的实时拓扑优化愿景出发，经由 EMsFEM 多尺度形函数学习正式提出问题无关机器学习，再进一步发展到子结构化三维大规模分析、复杂设计域、data-free 力学约束算子学习、显式晶格结构设计以及并行高性能拓扑优化。其核心价值在于把机器学习从“替代优化过程的黑箱预测器”转化为“嵌入有限元体系的局部分析算子”，这也是该方向区别于一般深度学习拓扑优化工作的关键所在。
