# 二维势能面扫描 (2D Potential Energy Surface)

## 概述

二维势能面（2D PES）扫描是研究分子几何结构-能量关系的重要方法，特别适用于：
- 基态和激发态的平衡构型对比
- 几何弛豫效应分析
- Stokes 位移的物理起源
- 化学反应路径探索

## 理论背景

### 势能面 (PES)
势能面描述分子能量 E 作为原子核坐标 {R} 的函数：
```
E = E(R₁, R₂, ..., R₃N)
```

对于 N 原子分子，PES 是 3N-6 维（非线性）或 3N-5 维（线性）的超曲面。

### 二维扫描
固定两个关键坐标（如键长 r₁, r₂），计算能量：
```
E(r₁, r₂) = E(r₁, r₂, 其他坐标固定)
```

## 网格设计策略

### 1. 范围选择

**原则**：必须包含基态和激发态的平衡构型

**示例（苯分子）**：
```python
# 基态平衡: C-C ≈ 1.40 Å, C-H ≈ 1.08 Å
# 激发态平衡: C-C ≈ 1.80 Å, C-H ≈ 0.85 Å

# 推荐范围
cc_range = np.linspace(1.25, 1.80, 15)  # C-C 键长
ch_range = np.linspace(0.85, 1.20, 15)  # C-H 键长
```

### 2. 网格密度

| 密度 | 点数 | 计算时间 | 适用场景 |
|------|------|----------|----------|
| 粗略 | 5×5 = 25 | ~5 分钟 | 快速预览 |
| 中等 | 7×7 = 49 | ~10 分钟 | 常规计算 |
| 精细 | 10×10 = 100 | ~20 分钟 | 精确分析 |
| 超精细 | 15×15 = 225 | ~40 分钟 | 发表级结果 |

### 3. 基组选择

| 基组 | 精度 | 速度 | 适用场景 |
|------|------|------|----------|
| STO-3G | 低 | 快 | 快速测试 |
| 3-21G | 中 | 中 | 常规计算 |
| cc-pVDZ | 高 | 慢 | 精确计算 |

## 完整代码模板

### 基础版本（7×7 网格）

```python
import numpy as np
from pyscf import gto, dft, tdscf

# 网格定义
cc_range = np.linspace(1.25, 1.80, 7)
ch_range = np.linspace(0.85, 1.20, 7)

gs_energies = []
s1_energies = []

for cc in cc_range:
    for ch in ch_range:
        # 构建分子
        atoms = []
        for i in range(6):
            angle = i * np.pi / 3
            atoms.append(f"C {cc*np.cos(angle):.6f} {cc*np.sin(angle):.6f} 0.0")
            atoms.append(f"H {(cc+ch)*np.cos(angle):.6f} {(cc+ch)*np.sin(angle):.6f} 0.0")
        
        mol = gto.M(atom='\n'.join(atoms), basis='3-21g', verbose=0)
        
        # 基态计算
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        e_gs = mf.kernel()
        
        # 激发态计算
        td = tdscf.TDDFT(mf)
        td.nstates = 3
        td.kernel()
        e_s1 = td.e[0] * 27.2114
        
        gs_energies.append(e_gs)
        s1_energies.append(e_s1)

# 重塑为 2D 数组
gs_2d = np.array(gs_energies).reshape(7, 7)
s1_2d = np.array(s1_energies).reshape(7, 7)
```

## 可视化方法

### 等高线图（推荐）

```python
import matplotlib.pyplot as plt

CC, CH = np.meshgrid(cc_range, ch_range)

fig, ax = plt.subplots(figsize=(10, 8))
contour = ax.contourf(CC, CH, gs_2d.T, levels=25, cmap='viridis')
plt.colorbar(contour, ax=ax, label='Relative Energy (Ha)')
ax.contour(CC, CH, gs_2d.T, levels=12, colors='white', alpha=0.4, linewidths=0.5)

# 标记最小值
min_idx = np.unravel_index(gs_2d.argmin(), gs_2d.shape)
ax.plot(cc_range[min_idx[0]], ch_range[min_idx[1]], 'r*', markersize=20)

ax.set_xlabel('C-C Bond Length (Å)', fontsize=14)
ax.set_ylabel('C-H Bond Length (Å)', fontsize=14)
ax.set_title('Ground State Potential Energy Surface', fontsize=16, fontweight='bold')
ax.set_aspect('equal', adjustable='box')
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig('2d_pes.png', dpi=300)
```

## 结果分析

### 找到最小值

```python
# 基态最小值
min_gs = np.unravel_index(gs_2d.argmin(), gs_2d.shape)
gs_min_cc = cc_range[min_gs[0]]
gs_min_ch = ch_range[min_gs[1]]

# 激发态最小值
min_s1 = np.unravel_index(s1_2d.argmin(), s1_2d.shape)
s1_min_cc = cc_range[min_s1[0]]
s1_min_ch = ch_range[min_s1[1]]

print(f"基态 S₀: C-C = {gs_min_cc:.4f} Å, C-H = {gs_min_ch:.4f} Å")
print(f"激发态 S₁: C-C = {s1_min_cc:.4f} Å, C-H = {s1_min_ch:.4f} Å")
```

### Stokes 位移

```python
# 垂直吸收（基态几何下的激发能）
vertical_abs = s1_2d[min_gs]

# 垂直发射（激发态几何下的激发能）
vertical_em = s1_2d.min()

# Stokes 位移
stokes = vertical_abs - vertical_em

print(f"垂直吸收: {vertical_abs:.2f} eV ({1240/vertical_abs:.1f} nm)")
print(f"垂直发射: {vertical_em:.2f} eV ({1240/vertical_em:.1f} nm)")
print(f"Stokes 位移: {stokes:.3f} eV ({stokes*8065.5:.0f} cm⁻¹)")
```

## 实际案例：苯分子

### 计算结果（15×15 网格，3-21G 基组）

**基态 S₀ 最小值**：
- C-C = 1.4071 Å
- C-H = 1.0750 Å
- 与实验值（C-C ≈ 1.40 Å）一致

**激发态 S₁ 最小值**：
- C-C = 1.8000 Å
- C-H = 0.8500 Å
- 激发能 = 2.90 eV (427 nm)

**几何变化**：
- Δ(C-C) = +0.393 Å (+27.9%)
- Δ(C-H) = -0.225 Å (-20.9%)

**Stokes 位移**：
- 垂直吸收：5.60 eV (222 nm)
- 垂直发射：2.90 eV (427 nm)
- Stokes 位移：2.70 eV (21,751 cm⁻¹)

---

**创建日期**：2026-03-09
**适用版本**：PySCF v2.12.1
