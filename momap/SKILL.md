# MOMAP Skill

**MOMAP (Molecular Materials Property Prediction Package)** - 分子材料性质预测软件包

## 功能概述

MOMAP 用于计算分子的光物理性质和传输性质：

### 光物理性质
- **荧光光谱** (Fluorescence Spectrum)
- **磷光光谱** (Phosphorescence Spectrum)
- **内转换速率** (Internal Conversion, IC)
- **系间交叉速率** (Intersystem Crossing, ISC)
- **辐射跃迁速率** (Radiative Rate)
- **振动分析** (Duschinsky Rotation Matrix)
- **电子-振动耦合** (Electron-Vibration Coupling, EVC)

### 传输性质
- **电荷传输计算** (Charge Transport)
- **转移积分** (Transfer Integral)
- **重组能** (Reorganization Energy)
- **随机游走** (Random Walk)

## 服务器信息

### 可用版本
```bash
module avail momap
# momap/2022B-openmpi
# momap/2024A-openmpi  ← 推荐使用最新版本
```

### 加载模块
```bash
module load momap/2024A-openmpi
```

## 基本使用流程

### 1. 准备量子化学计算文件

需要以下文件（来自 Gaussian/ORCA/Q-Chem）：
- **基态优化**: `gs.log`
- **激发态优化**: `es.log`
- **频率计算**: `gs_freq.log`, `es_freq.log`

### 2. 创建 MOMAP 输入文件

**momap.inp** 示例（荧光光谱计算）：
```bash
&control
    jobtype = 'fluor',        ! 计算类型
    qc_software = 'gaussian', ! 量子化学软件
    gs_file = 'gs.log',       ! 基态文件
    es_file = 'es.log',       ! 激发态文件
/
```

### 3. 运行 MOMAP

```bash
# 单核运行
momap < momap.inp > momap.out

# 并行运行 (MPI)
mpirun -np 8 momap < momap.inp > momap.out
```

## 常用计算类型

### 荧光光谱计算

**输入文件**:
```bash
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
    temperature = 298.0,      ! 温度 (K)
    spectrum_range = 1.0 3.5, ! 光谱范围 (eV)
    spectrum_step = 0.01,     ! 步长 (eV)
/
```

**输出**:
- `fluor_spectrum.dat` - 光谱数据
- `fluor_rate.dat` - 辐射跃迁速率

### 内转换 (IC) 速率

**输入文件**:
```bash
&control
    jobtype = 'ic',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 's1_freq.log',
    nacme_file = 'nacme.dat', ! 非绝热耦合矩阵元
/
```

### 系间交叉 (ISC) 速率

**输入文件**:
```bash
&control
    jobtype = 'isc',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 't1_freq.log',
    spin_orbit = .true.,
/
```

## 快速参考卡

| 计算类型 | jobtype | 所需文件 | 输出文件 |
|---------|---------|---------|---------|
| 荧光光谱 | `fluor` | gs_freq.log, es_opt.log | fluor_spectrum.dat |
| 磷光光谱 | `phos` | gs_freq.log, t1_opt.log | phos_spectrum.dat |
| 内转换 | `ic` | s0_freq.log, s1_freq.log, nacme.dat | ic_rate.dat |
| 系间交叉 | `isc` | s0_freq.log, t1_freq.log | isc_rate.dat |
| 电荷传输 | `transfer` | mol1.log, mol2.log, dimer.log | transfer_rate.dat |

**常用命令**:
```bash
module load momap/2024A-openmpi    # 加载模块
momap < momap.inp > momap.out      # 单核运行
mpirun -np 8 momap < momap.inp     # 并行运行
```

## 典型工作流程示例

### 计算荧光量子产率

```bash
# 1. 运行 Gaussian
g16 < gs_opt.com > gs_opt.log
g16 < gs_freq.com > gs_freq.log
g16 < es_opt.com > es_opt.log

# 2. 运行 MOMAP (荧光)
module load momap/2024A-openmpi
cat > momap_fluor.inp << 'EOF'
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
/
EOF

momap < momap_fluor.inp > fluor.out

# 3. 提取辐射速率
grep "Radiative rate" fluor.out
```

## 数据分析

### 绘制光谱
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fluor_spectrum.dat')
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity (a.u.)')
plt.title('Fluorescence Spectrum')
plt.savefig('fluor_spectrum.png')
```

## 参考文献

1. **MOMAP 论文**:
   - Shuai, Z. G.; et al. *J. Chem. Phys.* **2009**, 131, 224106.

2. **荧光理论**:
   - Niu, Y.; Peng, Q.; Shuai, Z. *Sci. China Chem.* **2016**, 59, 1361.

## 引用 MOMAP

如果在研究中使用 MOMAP，请引用：
```
MOMAP, Version 2024A, Hongzhiwei Technology (Shanghai) Co., Ltd and
Z.G. Shuai Group, 2024.
```

---

**官方资源**: http://www.momap.cn
