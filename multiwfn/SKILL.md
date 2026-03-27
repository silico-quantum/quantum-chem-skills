---
name: multiwfn
version: 1.0.0
description: Multiwfn 波函数分析工具。支持轨道分析、光谱绘制、拓扑分析等。
homepage: http://sobereva.com/multiwfn
metadata:
  category: computational_chemistry
  tags: [wavefunction_analysis, orbitals, spectrum, topology]
---

# Multiwfn 波函数分析技能

## 概述

Multiwfn 是一个强大的波函数分析工具，用于分子轨道分析、光谱绘制、电子密度拓扑分析等。

## 安装位置

- **本地:** `/opt/homebrew/bin/Multiwfn`
- **版本:** 3.8 (2026-Jan-7)
- **开发者:** Tian Lu (北京科音自然科学研究中心)

## 引用要求

使用 Multiwfn 必须引用:
1. Tian Lu, Feiwu Chen, J. Comput. Chem., 33, 580 (2012) DOI: 10.1002/jcc.22885
2. Tian Lu, J. Chem. Phys., 161, 082503 (2024) DOI: 10.1063/5.0216272

## 支持的输入文件

| 格式 | 说明 |
|------|------|
| .fchk / .fch | Gaussian 格式化检查点（推荐） |
| .chk | Gaussian 检查点 |
| .wfn | 波函数文件 |
| .molden | Molden 格式 |
| .xyz | 分子坐标 |
| .cub / .cube | Cube 文件 |

## 主菜单功能

### 0 - 查看分子结构和轨道

交互式 3D 分子可视化，可查看分子轨道。

### 1 - 点属性计算

计算空间中某点的各种属性（电子密度、梯度、拉普拉斯等）。

### 2 - 拓扑分析 (AIM)

分析电子密度的临界点 (CP)，用于化学键分析。

### 3 - 线属性图

沿某条线绘制属性变化（如键轴上的电子密度）。

### 4 - 平面图

在某个平面上绘制等高线图（如电子密度、轨道）。

### 5 - 空间网格数据

计算 3D 网格数据，生成 Cube 文件。

### 6 - 检查/修改波函数

查看轨道信息、能量、修改波函数。

### 7 - 布居分析

Mulliken、Hirshfeld、ADCH 等原子电荷计算。

### 8 - 轨道成分分析

分析分子轨道的原子/片段贡献。

### 9 - 键级分析

计算各种键级（Mayer、Wiberg 等）。

### 10 - DOS/PDOS

态密度、投影态密度、光电子能谱。

### 11 - 光谱绘制

绘制 IR、Raman、UV-Vis、ECD、NMR 光谱。

### 12 - 分子表面分析

计算分子表面的静电势、平均局部离子化能等。

### 18 - 激发态分析

电子激发分析，空穴-电子分析，电荷转移分析。

### 20 - 弱相互作用可视化

RDG、NCI 分析，可视化分子间相互作用。

## 常用功能编号速查

| 功能 | 菜单路径 |
|------|----------|
| 查看轨道能量 | 0 → q |
| UV-Vis 光谱 | 11 → 1 |
| IR 光谱 | 11 → 2 |
| 原子电荷 | 7 → 1/2/3 |
| 轨道成分 | 8 → 1 |
| 键级 | 9 → 1 |
| DOS | 10 → 1 |
| 分子表面 | 12 |
| RDG 分析 | 20 → 1 |
| 激发态分析 | 18 |

## 与 Gaussian 配合使用

1. 用 Gaussian 完成 DFT/TDDFT 计算
2. 用 formchk 转换 .chk → .fchk
3. 用 Multiwfn 分析 .fchk 文件

## 从远程服务器获取文件

```bash
# 下载 fchk 文件到本地
scp -P 8722 openclaw@124.16.75.110:/path/to/molecule.fchk ./
```

## 官方资源

- 官网: http://sobereva.com/multiwfn
- 英文论坛: http://sobereva.com/wfnbbs
- 中文论坛: http://bbs.keinsci.com/wfn

---

*创建日期: 2026-03-10*
