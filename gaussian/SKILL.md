---
name: gaussian
version: 1.0.0
description: Gaussian 量子化学计算软件使用指南。支持 DFT、TDDFT、频率计算、几何优化等。
homepage: https://gaussian.com
metadata:
  category: computational_chemistry
  tags: [quantum_chemistry, DFT, TDDFT, frequency, optimization]
---

# Gaussian 量子化学计算技能

## 概述

Gaussian 是业界标准的量子化学计算软件，用于电子结构计算、几何优化、频率分析、激发态计算等。

## 安装位置

- **远程服务器:** marcus (124.16.75.110:8722 → marcus)
- **模块加载:** `module load gaussian/g16.c01-avx2` 或 `gaussian/g16.c01-avx`
- **版本:** Gaussian 16 Rev. C.01

## 连接方式

```bash
# 登录流程
ssh -p 8722 openclaw@124.16.75.110  # 密码: ucas@204
ssh marcus                           # 密码: ucas@204
```

## Gaussian 输入文件格式 (.gjf)

### 基本结构

```
%chk=filename.chk          # 检查点文件
%mem=60GB                  # 内存分配
%nproc=64                  # CPU 核心数
# 关键词行
标题行
电荷 自旋多重度
原子坐标 (笛卡尔坐标)
```

## 常用关键词

### 1. 计算类型

| 关键词 | 说明 |
|--------|------|
| `sp` | 单点能量计算 |
| `opt` | 几何优化 |
| `freq` | 频率分析 |
| `opt freq` | 优化 + 频率 |

### 2. 方法

| 方法 | 说明 |
|------|------|
| `b3lyp` | B3LYP 杂化泛函 |
| `pbe0` | PBE0 杂化泛函 |
| `m06-2x` | M06-2X 泛函 |
| `cam-b3lyp` | 范围分离泛函 |

### 3. 基组

| 基组 | 说明 |
|------|------|
| `6-31g*` | 添加极化函数 |
| `6-31g**` | 氢原子极化 |
| `6-311g**` | 三ζ基组 |
| `def2-tzvp` | def2 三ζ |

### 4. 激发态计算 (TDDFT)

```
# td=(nstates=10) b3lyp/6-31g*
```

- `nstates=N`: 计算 N 个激发态
- `root=N`: 优化第 N 个激发态

### 5. 溶剂效应

```
# b3lyp/6-31g* scrf=(smd,solvent=water)
```

## 提交作业

```bash
sbatch jobname.slurm jobname.gjf
```

## formchk 命令

```bash
formchk filename.chk filename.fchk
```

## 常见任务模板

### 基态优化 + 频率

```
%nproc=64
%chk=molecule.chk
%mem=60GB
# opt freq b3lyp/6-31g**

Title
0 1
C 0.0 0.0 0.0
```

### TDDFT 激发态

```
%nproc=64
%chk=molecule.chk
%mem=60GB
# opt freq td=(nstates=5) b3lyp/6-31g** geom=check

Title
0 1
```

## 输出文件

| 文件 | 说明 |
|------|------|
| `.log` | 主输出文件 |
| `.chk` | 检查点文件 |
| `.fchk` | 格式化检查点 |

