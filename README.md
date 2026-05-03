# Phase Field Fracture Solver - AT2 Model

## 简介

本代码实现了一个基于AT2的断裂相场模型，用于模拟反平面剪切（Mode III）载荷下的断裂行为。

## 理论背景

### 控制方程

**位移平衡方程：** $$-\nabla \cdot \sigma = 0$$

其中应力 $\sigma = g(d) G \varepsilon$，$g(d) = (1-d)^2$ 为损伤退化函数。

**相场演化方程：** $$\frac{G_c}{l_0} d - G_c l_0 \Delta d = -g'(d) H$$

其中 $H = \max_{[0,t]}(\psi_{elastic})$ 是历史变量，$\psi_{elastic} = \frac{1}{2} G |\varepsilon|^2$ 是弹性应变能密度。

### AT2模型特点

- 使用二次型损伤函数：$(1-d)^2+10^{-10}$
- 正则化断裂能：$\psi_{frac} = \frac{G_c}{l_0} d^2 + G_c l_0 |\nabla d|^2$
- 适合模拟脆性断裂

## 文件结构

```
MATLAB/
├── PhaseFieldFractureSolver.m    # 核心求解器类
├── runRectangleFractureExample.m  # 带预制裂纹的算例
└── README.md                     # 本文件
```

## 快速开始

### 基本用法

```
% 创建网格
mesh = PhaseFieldFractureSolver.createMesh(nx, ny, Lx, Ly);

% 创建材料
mat = PhaseFieldFractureSolver.createMaterial(E, nu, Gc, l0);

% 设置边界条件
fixedDOFs = [...];      % 固定节点
prescribedDOFs = [...];  % 指定位移节点
prescribedValues = [...];% 边界值

% 求解
[u, d, mesh, history, converged] = PhaseFieldFractureSolver.solveStaggered(...
    mesh, mat, u0, d0, fixedDOFs, prescribedDOFs, prescribedValues, tolerance, maxIter);
```

### 运行示例

```
% 运行矩形域断裂算例
runRectangleFractureExample;

% 运行弹性验证
testElasticVerification;
```

## 输入参数

| 参数   | 说明           | 单位 |
| :----- | :------------- | :--- |
| E      | 弹性模量       | Pa   |
| nu     | 泊松比         | -    |
| Gc     | 临界能量释放率 | J/m² |
| l0     | 长度尺度参数   | m    |
| nx, ny | 网格节点数     | -    |
| Lx, Ly | 域尺寸         | m    |

### 参数选择建议

- 长度尺度：$l_0 \geq 2 \cdot \max(dx, dy)$
- 网格尺寸：建议 $dx, dy \leq l_0 / 2$
- 时间步长：根据收敛性调整，建议 $\Delta U \leq 1 \times 10^{-6}$ m

## 输出结果

- `u`: 位移场向量
- `d`: 相场（损伤）向量
- `mesh`: 更新后的网格结构体
- `history`: 收敛历史记录

## 可视化

代码支持以下可视化功能：

```
% 位移场云图
PhaseFieldFractureSolver.plotDisplacement(mesh, u);

% 相场云图
PhaseFieldFractureSolver.plotPhaseField(mesh, d);

% 应力云图
PhaseFieldFractureSolver.plotStress(mesh, mat, u, d);

% 保存VTK格式
PhaseFieldFractureSolver.saveToVTK(mesh, mat, u, d, 'output');
```

## 依赖

- MATLAB R2016b 或更高版本
- 无需额外工具箱

## 引用

如果使用本代码，请引用：

```
@book{shen2026fracture_cn,
    title={变分断裂相场法基础},
    author={沈泳星},
    publisher={清华大学出版社},
    year={2026}
}

@book{shen2026fracture_en,
    title={Variational Phase Field Method for Fracture: Fundamentals},
    author={Shen, Yongxing},
    publisher={Tsinghua University Press},
    year={2026},
    note={in Chinese}
}
```

## 许可

MIT License

## 作者

Yongxing Shen

## 联系方式

yongxing.shen@sjtu.edu.cn
