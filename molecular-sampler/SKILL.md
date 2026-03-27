---
name: molecular-sampler
version: 1.0.0
description: Sample monomers and multi-molecule complexes from Gaussian ONIOM or XYZ files. Extract all monomers, then sample dimers/trimers/tetramers/pentamers using distance-sorted nearest neighbors.
homepage: https://github.com/STOKES-DOT/molecular-sampler
metadata:
  category: chemistry
  emoji: 🧪
  tags:
    - chemistry
    - sampling
    - molecules
    - gaussian
    - oniom
---

# Molecular Sampler

从 Gaussian ONIOM 或 XYZ 文件中采样分子结构。

## 功能

1. **单体提取**：识别并保存所有独立分子
2. **多聚体采样**：基于距离排序，选择最近的相邻分子
   - 二聚体 (dimers): 2 个分子
   - 三聚体 (trimers): 3 个分子
   - 四聚体 (tetramers): 4 个分子
   - 五聚体 (pentamers): 5 个分子

## 使用方法

### 基本用法

```bash
python3 molecular_sampler.py <input_file> [options]
```

### 参数

- `input_file`: Gaussian GJF 或 XYZ 文件路径
- `--output-dir`: 输出目录 (默认: `./molecular_samples`)
- `--samples`: 每种多聚体的采样数量 (默认: 20)
- `--layer`: 选择层 ('H', 'L', 或 'all') (默认: 'L')

### 示例

```bash
# 采样 L 层分子
python3 molecular_sampler.py guest_monomer.gjf --layer L --samples 20

# 采样所有层分子
python3 molecular_sampler.py structure.xyz --layer all --output-dir ./my_samples
```

## 输出结构

```
molecular_samples/
├── monomers/           # 所有单体分子
│   ├── monomer_01.xyz
│   ├── monomer_02.xyz
│   └── ...
├── dimers/             # 二聚体 (20个样本)
├── trimers/            # 三聚体 (20个样本)
├── tetramers/          # 四聚体 (20个样本)
├── pentamers/          # 五聚体 (20个样本)
└── sampling_summary.txt
```

## 采样策略

### 分子识别
- 使用 Union-Find 算法识别连通分子
- 基于共价半径检测化学键
- 分子间距离 > 3Å 视为独立分子

### 多聚体采样
1. 计算所有分子中心点之间的距离
2. 对每个分子，按距离排序邻居
3. 选择最近的 N-1 个邻居形成 N 聚体
4. 从前 20 个分子开始采样，确保多样性

## XYZ 文件格式

标准 XYZ 格式：
```
<atom_count>
<comment>
<element> <x> <y> <z>
...
```

## 依赖

- Python 3.7+
- NumPy
- 标准库: `re`, `collections`, `os`, `argparse`

## 工作流程

1. **解析文件**：读取 GJF/XYZ 文件，提取原子坐标
2. **分子识别**：使用化学键检测识别独立分子
3. **距离计算**：计算所有分子中心点之间的欧氏距离
4. **邻居排序**：对每个分子按距离排序邻居列表
5. **采样**：选择最近的邻居形成多聚体
6. **输出**：生成 XYZ 文件和摘要

## 适用场景

- 蛋白质-配体复合物采样
- 晶体结构中的分子簇提取
- 分子动力学轨迹采样
- 量子化学计算体系准备

## 注意事项

- 确保输入文件中的分子间距合理 (>3Å)
- 大体系 (>10000 原子) 可能需要较长处理时间
- ONIOM 文件需包含正确的层标记 (H/L)

## 版本历史

- **v1.0.0** (2026-03-03): 初始版本
  - 支持单体的完整提取
  - 距离排序的多聚体采样
  - 标准 XYZ 输出格式
