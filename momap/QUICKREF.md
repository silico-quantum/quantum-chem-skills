# MOMAP 快速参考卡

## 模块加载

```bash
module load momap/2024A-openmpi
```

## 常用命令

| 命令 | 说明 |
|------|------|
| `momap < momap.inp > momap.out` | 单核运行 |
| `mpirun -np 8 momap < momap.inp > momap.out` | 8 核并行运行 |
| `momap-ht < momap.inp > momap.out` | Herzberg-Teller 计算 |
| `momap-sos < momap.inp > momap.out` | Sum-over-states 计算 |

## 输入文件模板

### 荧光光谱 (Fluorescence)

```bash
&control
    jobtype = 'fluor',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 'es_opt.log',
    temperature = 298.0,
    spectrum_range = 3.0 6.0,
    spectrum_step = 0.01,
/
```

### 磷光光谱 (Phosphorescence)

```bash
&control
    jobtype = 'phos',
    qc_software = 'gaussian',
    gs_file = 'gs_freq.log',
    es_file = 't1_opt.log',
    spin_orbit = .true.,
/
```

### 内转换 (Internal Conversion)

```bash
&control
    jobtype = 'ic',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 's1_freq.log',
    nacme_file = 'nacme.dat',
/
```

### 系间交叉 (Intersystem Crossing)

```bash
&control
    jobtype = 'isc',
    qc_software = 'gaussian',
    gs_file = 's0_freq.log',
    es_file = 't1_freq.log',
    spin_orbit = .true.,
/
```

### 电荷传输 (Charge Transfer)

```bash
&control
    jobtype = 'transfer',
    qc_software = 'gaussian',
    molecule1 = 'mol1.log',
    molecule2 = 'mol2.log',
    dimer = 'dimer.log',
/
```

## 常用参数

### 控制参数

| 参数 | 类型 | 说明 | 默认值 |
|------|------|------|--------|
| `jobtype` | string | 计算类型 | 必需 |
| `qc_software` | string | 量子化学软件 | 'gaussian' |
| `gs_file` | string | 基态文件 | 必需 |
| `es_file` | string | 激发态文件 | 必需 |
| `temperature` | real | 温度 (K) | 298.0 |
| `spectrum_range` | real(2) | 光谱范围 (eV) | 1.0 5.0 |
| `spectrum_step` | real | 光谱步长 (eV) | 0.01 |
| `max_time` | real | 最大时间 (fs) | 1000.0 |
| `time_step` | real | 时间步长 (fs) | 0.5 |
| `convergence` | real | 收敛阈值 | 1.0E-8 |

### 高级参数

| 参数 | 类型 | 说明 | 默认值 |
|------|------|------|--------|
| `ht_effect` | logical | Herzberg-Teller 效应 | .false. |
| `ht_order` | integer | HT 展开阶数 | 1 |
| `sos_method` | logical | Sum-over-states 方法 | .false. |
| `nstates` | integer | SOS 激发态数 | 10 |
| `spin_orbit` | logical | 自旋轨道耦合 | .false. |
| `memory` | string | 内存限制 | '4GB' |

## 输出文件

| 文件名 | 内容 | 格式 |
|--------|------|------|
| `*_spectrum.dat` | 光谱数据 | Energy (eV), Intensity |
| `*_rate.dat` | 速率常数 | Temperature (K), Rate (s^-1) |
| `correlation.dat` | 相关函数 | Time (fs), Correlation |
| `*.out` | MOMAP 输出 | 文本 |

## Gaussian 输入模板

### 基态优化 + 频率

```bash
#P B3LYP/6-31G* Opt Freq

标题行

0 1
分子坐标
```

### 激发态优化 (TD-DFT)

```bash
#P TD(NStates=3) B3LYP/6-31G* Opt

标题行

0 1
分子坐标
```

### 三重态优化

```bash
#P UB3LYP/6-31G* Opt

标题行

0 3
分子坐标
```

### NACME 计算

```bash
#P TD(NStates=3) B3LYP/6-31G* Freq NoSymm

标题行

0 1
分子坐标
```

## 快速分析命令

### 提取速率

```bash
# 荧光速率
grep "Radiative rate" fluor.out

# IC 速率
grep "IC rate" ic.out

# ISC 速率
grep "ISC rate" isc.out
```

### 绘制光谱

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fluor_spectrum.dat')
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity (a.u.)')
plt.savefig('spectrum.png')
```

### 检查收敛性

```bash
# 查看相关函数
tail -n 20 correlation.dat

# 检查收敛标志
grep "CONVERGED" momap.out
```

## 常见错误

| 错误信息 | 原因 | 解决方法 |
|---------|------|---------|
| `Insufficient memory` | 内存不足 | `export MOMAP_MEM=8GB` |
| `Frequency not found` | 频率文件格式错误 | 检查 Gaussian 输出 |
| `Correlation not converged` | 时间不够 | 增加 `max_time` |
| `Negative frequency` | 有虚频 | 重新优化几何 |
| `File not found` | 文件路径错误 | 检查文件名和路径 |

## 量子产率计算

```python
# 提取速率常数
kr = 1.5e8   # 荧光速率 (s^-1)
kic = 5.0e7  # IC 速率 (s^-1)
kisc = 3.0e7 # ISC 速率 (s^-1)

# 计算量子产率
qy = kr / (kr + kic + kisc)

# 计算寿命
tau = 1.0 / (kr + kic + kisc)

print(f"量子产率: {qy:.4f} ({qy*100:.2f}%)")
print(f"寿命: {tau*1e9:.2f} ns")
```

## Slurm 作业模板

```bash
#!/bin/bash
#SBATCH --job-name=momap
#SBATCH --partition=X48Cv4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00

module load momap/2024A-openmpi
mpirun -np $SLURM_NTASKS momap < momap.inp > momap.out
```

## 兼容的量子化学软件

| 软件 | 输出文件 | 支持度 |
|------|---------|--------|
| Gaussian 16/09 | .log | ★★★★★ |
| ORCA | .out | ★★★★☆ |
| Q-Chem | .out | ★★★★☆ |
| DALTON | .out | ★★★☆☆ |
| TURBOMOLE | .log | ★★★☆☆ |

## 性能优化

### 并行计算
```bash
# MPI 并行
mpirun -np 16 momap < momap.inp > momap.out

# Slurm 作业
#SBATCH --ntasks-per-node=16
mpirun -np $SLURM_NTASKS momap < momap.inp
```

### 内存优化
```bash
# 限制内存使用
&control
    memory = '16GB',
/

# 或环境变量
export MOMAP_MEM=16GB
```

## 有用的一行命令

```bash
# 转换光谱单位 (eV → nm)
awk '{print 1239.8/$1, $2}' fluor_spectrum.dat > fluor_spectrum_nm.dat

# 找到光谱峰值
awk 'NR>1 && $2>max{max=$2; peak=$1} END{print "Peak:", peak, "eV"}' fluor_spectrum.dat

# 计算积分强度
awk '{sum+=$2*0.01} END{print "Integrated intensity:", sum}' fluor_spectrum.dat

# 提取所有速率
grep -E "(Radiative|IC|ISC) rate" *.out
```

---

**完整文档**: `/Users/molbot/.openclaw/workspace/skills/momap/SKILL.md`
**使用示例**: `/Users/molbot/.openclaw/workspace/skills/momap/EXAMPLES.md`
