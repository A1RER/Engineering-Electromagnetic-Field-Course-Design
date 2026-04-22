# 电/磁偶极子空间场强数值仿真

**Electric / Magnetic Dipole Spatial Field Strength Numerical Simulation**

本项目提供电偶极子与磁偶极子空间场强的数值计算与可视化工具，同时提供 **MATLAB** 和 **C** 两种实现。支持静态近场与时谐场（Hertzian 偶极子）两种计算模式，包含解析校验模块，结果可导出为 CSV 供 Origin 等软件进一步分析。

---

## 目录

- [物理背景](#物理背景)
- [项目结构](#项目结构)
- [快速开始](#快速开始)
  - [MATLAB](#matlab)
  - [C](#c)
- [功能特性](#功能特性)
- [参数说明](#参数说明)
- [输出说明](#输出说明)
- [校验与精度](#校验与精度)
- [依赖与环境](#依赖与环境)
- [License](#license)

---

## 物理背景

### 静态近场

**电偶极子 E 场：**

$$\mathbf{E}(\mathbf{r}) = \frac{1}{4\pi\varepsilon_0} \cdot \frac{3\hat{r}(\mathbf{p}\cdot\hat{r}) - \mathbf{p}}{r^3}$$

**磁偶极子 H / B 场：**

$$\mathbf{H}(\mathbf{r}) = \frac{1}{4\pi} \cdot \frac{3\hat{r}(\mathbf{m}\cdot\hat{r}) - \mathbf{m}}{r^3}, \quad \mathbf{B} = \mu_0 \mathbf{H}$$

### 时谐场（Hertzian 偶极子，$e^{j\omega t}$ 约定）

$$\mathbf{E}(\mathbf{r}) = \frac{e^{-jkr}}{4\pi\varepsilon_0} \left[ \frac{k^2}{r}\,\mathbf{p}_\perp + \frac{1+jkr}{r^3}\,\mathbf{p}_\parallel \right]$$

其中 $\mathbf{p}_\perp$ 为横向（辐射）分量，$\mathbf{p}_\parallel$ 为纵向（反应）分量。磁偶极子时谐场与之具有对偶结构。

### 奇异点正则化

在 $r \to 0$ 处通过软化因子 $\varepsilon_s$ 避免数值奇点：

$$r^2 \;\longrightarrow\; r^2 + \varepsilon_s^2, \quad \varepsilon_s \approx 10^{-4} \times L$$

---

## 项目结构

```
.
├── matlab/
│   ├── main_dipole_sim.m          # 主程序：仿真控制、网格生成、可视化、数据导出
│   ├── electric_dipole_field.m    # 电偶极子 E 场矢量化计算函数
│   ├── magnetic_dipole_field.m    # 磁偶极子 H/B 场矢量化计算函数
│   └── verify_dipole.m            # 解析交叉校验与收敛性分析函数
├── c/
│   └── dipole_field.c             # C99 实现：三维网格场强计算与 CSV 导出
├── LICENSE                        # MIT License
└── README.md
```

---

## 快速开始

### MATLAB

**环境要求：** MATLAB R2019b 及以上

1. 打开 MATLAB，将工作目录切换至 `matlab/`。
2. 直接运行主程序：
   ```matlab
   run('main_dipole_sim.m')
   % 或在编辑器中按 F5
   ```
3. 程序将自动弹出以下可视化窗口：
   - **场强幅值热力图**（xz 平面，log 色阶）
   - **归一化矢量箭头图**（Quiver + 等幅线）
   - **极轴衰减特性曲线**（log-log，含 1/r³ 解析参考线）
   - **场线拓扑图**（Streamline）
   - **三维正交切面渲染**（3D Slice，N3 ≤ 50 时自动生成）
   - **收敛性分析图**（电偶极子模式下自动生成）
4. 计算结果自动导出至项目根目录 `dipole_field_data.csv`。

**切换偶极子类型或计算模式：** 编辑 `main_dipole_sim.m` 中 Section 1 的参数：

```matlab
DIPOLE_TYPE = 'electric';   % 'electric' | 'magnetic'
SIM_MODE    = 'static';     % 'static'   | 'timeharmonic'
```

---

### C

**环境要求：** GCC / MinGW，支持 C99

**编译：**

```bash
# 标准编译
gcc -O2 -std=c99 -o dipole_field c/dipole_field.c -lm

# 启用 OpenMP 并行加速（可选）
gcc -O2 -std=c99 -fopenmp -o dipole_field c/dipole_field.c -lm
```

**运行：**

```bash
# 使用默认参数（电偶极子，静态场，40³ 网格，±5 nm 计算域）
./dipole_field

# 完整参数：类型  模式  网格点数/轴  计算域半径(nm)
./dipole_field e s 40 5.0    # 电偶极子，静态场
./dipole_field m s 60 5.0    # 磁偶极子，静态场
./dipole_field e t 40 5.0    # 电偶极子，时谐场（1 GHz）
```

| 参数 | 可选值 | 说明 |
|------|--------|------|
| `type` | `e`（默认）/ `m` | 电偶极子 / 磁偶极子 |
| `mode` | `s`（默认）/ `t` | 静态近场 / 时谐场 |
| `N` | 整数，2–500（默认 40） | 每轴网格点数，总计 N³ 个点 |
| `L_nm` | 浮点数（默认 5.0） | 计算域半径 [nm] |

运行后输出 `dipole_field_data.csv`，同时在终端打印自检结果与进度。

---

## 功能特性

| 特性 | MATLAB | C |
|------|:------:|:-:|
| 电偶极子 E 场 | ✅ | ✅ |
| 磁偶极子 H / B 场 | ✅ | ✅ |
| 静态近场模式 | ✅ | ✅ |
| 时谐场模式 | ✅ | ✅ |
| 全矢量化计算（无循环） | ✅ | — |
| OpenMP 并行加速 | — | ✅（可选） |
| 奇异点正则化 | ✅ | ✅ |
| 解析校验（极轴 / 赤道 / 斜率拟合） | ✅ | ✅（自检）|
| 2D 热力图 / 矢量图 / 流线图 | ✅ | — |
| 3D 正交切面渲染 | ✅ | — |
| CSV 数据导出 | ✅ | ✅ |

---

## 参数说明

以下为 MATLAB 主程序中的关键可调参数（`main_dipole_sim.m` Section 1）：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `DIPOLE_TYPE` | `'electric'` | 偶极子类型 |
| `SIM_MODE` | `'static'` | 计算模式 |
| `p_vec` | `[0, 0, 1e-30]` C·m | 电偶极子矩矢量（支持任意方向） |
| `m_vec` | `[0, 0, 1e-20]` A·m² | 磁偶极子矩矢量 |
| `r_src` | `[0, 0, 0]` | 偶极子源点坐标 [m] |
| `freq` | `1e9` Hz | 时谐场频率 |
| `L` | `5e-9` m | 计算域半径（±5 nm） |
| `N2` | `120` | 2D 网格每边点数 |
| `N3` | `40` | 3D 网格每边点数 |
| `eps_s` | `L × 1e-4` | 奇异点软化因子 |

---

## 输出说明

### CSV 文件格式（`dipole_field_data.csv`）

```
x_nm, y_nm, z_nm, Fx, Fy, Fz, F_mag
```

| 列 | 单位 | 说明 |
|----|------|------|
| `x_nm`, `y_nm`, `z_nm` | nm | 观测点坐标 |
| `Fx`, `Fy`, `Fz` | V/m 或 A/m | 场矢量分量（实部） |
| `F_mag` | V/m 或 A/m | 场强模值 $\|\mathbf{F}\|$ |

该格式与 Origin、Python（pandas/numpy）等工具直接兼容。

---

## 校验与精度

`verify_dipole.m` 对 z 方向电偶极子执行三项独立校验：

| 校验项 | 参考公式 | 通过阈值 |
|--------|---------|---------|
| 极轴（z 轴）相对误差 | $E_z = 2p / (4\pi\varepsilon_0 r^3)$ | < 0.1% |
| 赤道面（x 轴）相对误差 | $E_z = -p / (4\pi\varepsilon_0 r^3)$ | < 0.1% |
| log-log 衰减斜率拟合 | 理论斜率 = −3.000 | 偏差 < 0.05 |

C 版本在 `main()` 入口处执行内联自检，于极轴 $r = L/2$ 处对比解析值，相对误差需小于 0.1%。

---

## 依赖与环境

### MATLAB
- MATLAB R2019b 及以上
- 无需任何额外工具箱

### C
- GCC 5.0+ 或 MinGW（Windows）
- C 标准库：`math.h`、`stdio.h`、`stdlib.h`
- OpenMP（可选，用于并行加速）

---

## License

本项目基于 [MIT License](LICENSE) 开源。

Copyright (c) 2026 Leslie Shen
