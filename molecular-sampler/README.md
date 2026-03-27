# Molecular Sampler

🧪 **从 Gaussian ONIOM 或 XYZ 文件中采样分子结构**

## 快速开始

```bash
python3 /Users/molbot/.openclaw/skills/molecular-sampler/molecular_sampler.py <input_file> [options]
```

## 功能

1. **单体提取**：识别并保存所有独立分子（63 个 L 层分子）
2. **多聚体采样**：基于距离排序，选择最近的相邻分子
   - 二聚体 (dimers): 2 个最近分子
   - 三聚体 (trimers): 3 个最近分子
   - 四聚体 (tetramers): 4 个最近分子
   - 五聚体 (pentamers): 5 个最近分子

## 采样标准（v1.0.0）

### 分子识别
- **算法**: Union-Find 连通分量检测
- **键检测**: 基于共价半径 (1.3× tolerance)
- **最小分子**: 5 个原子

### 多聚体采样
- **策略**: 距离排序的最近邻居
- **样本数**: 每种类型 20 个
- **多样性**: 从前 20 个不同分子开始采样

### 输出格式
- **文件格式**: 标准 XYZ
- **命名**: `{type}_{number:02d}.xyz`
- **注释**: 包含分子数和原子数信息

## 使用示例

### 示例 1: 采样 L 层分子
```bash
python3 molecular_sampler.py guest_monomer.gjf \
  --layer L \
  --samples 20 \
  --output-dir ./molecular_samples
```

### 示例 2: 采样所有层
```bash
python3 molecular_sampler.py structure.xyz \
  --layer all \
  --samples 30
```

## 输出结构

```
molecular_samples/
├── monomers/           # 所有单体 (63 个)
│   ├── monomer_01.xyz
│   ├── monomer_02.xyz
│   └── ...
├── dimers/             # 二聚体 (20 个)
├── trimers/            # 三聚体 (20 个)
├── tetramers/          # 四聚体 (20 个)
├── pentamers/          # 五聚体 (20 个)
└── sampling_summary.txt
```

## 距离统计（示例）

- **最小距离**: 6.28 Å
- **最大距离**: 50.81 Å
- **平均距离**: 25.97 Å
- **中位数距离**: 25.95 Å

## 依赖

- Python 3.7+
- NumPy
- 标准库: `re`, `collections`, `os`, `argparse`

## 工作流程

1. 解析 GJF/XYZ 文件
2. 识别独立分子（化学键检测）
3. 计算分子中心点距离
4. 按距离排序邻居
5. 采样多聚体
6. 生成 XYZ 文件

## 创建日期

2026-03-03 - 为 Yuan Jiao 的分子采样工作创建

## 版本历史

- **v1.0.0** (2026-03-03): 初始版本
  - 完整的单体提取
  - 距离排序的多聚体采样
  - 标准 XYZ 输出
