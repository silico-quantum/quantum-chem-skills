# PySCF 发射光谱计算指南

## 概述

发射光谱（Fluorescence Spectrum）是分子从激发态回到基态时发射的光谱。与吸收光谱不同，发射光谱需要在**激发态优化后的几何结构**下计算。

## 物理过程

```
基态 (S₀) --[吸收]--> 激发态 (S_n)
  ↓                      ↓
基态几何              激发态几何 (弛豫)
                         ↓
                      [发射]
                         ↓
                      基态 (S₀)
```

## 关键概念

### 1. Stokes 位移
- **定义**: 吸收能量 - 发射能量
- **原因**: 激发态几何弛豫
- **公式**: ΔE_Stokes = E_abs - E_em

### 2. 垂直跃迁
- **垂直吸收**: 在基态几何下计算激发能
- **垂直发射**: 在激发态几何下计算发射能
- Franck-Condon 原理：电子跃迁比核运动快得多

### 3. 几何弛豫
激发态分子会经历:
1. 电子激发（~飞秒）
2. 振动弛豫（~皮秒）
3. 几何重排（~纳秒）
4. 荧光发射（~纳秒）

## PySCF 计算流程

### 标准流程

```python
from pyscf import gto, dft, tdscf
from pyscf.geomopt import geometric_solver

# 1. 定义分子（基态几何）
mol_gs = gto.M(atom='...', basis='cc-pvdz')

# 2. 基态 DFT 计算
mf_gs = dft.RKS(mol_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# 3. 计算吸收光谱（基态几何）
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()
absorption_energies = td_abs.e * 27.2114  # eV

# 4. 优化第一激发态几何
# 方法 1: 使用 TDA 激发态优化
from pyscf import tdscf
td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
mol_ex = geometric_solver.optimize(td_s1, state=0)

# 方法 2: 手动调整几何（近似）
# 激发态通常键长略长

# 5. 在激发态几何下做基态计算
mf_ex = dft.RKS(mol_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# 6. 计算发射光谱（激发态几何）
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()
emission_energies = td_em.e * 27.2114  # eV

# 7. 计算 Stokes 位移
stokes_shift = absorption_energies[0] - emission_energies[0]
```

## 完整示例：苯分子

```python
#!/usr/bin/env python3
"""
苯分子发射光谱计算示例
"""

import numpy as np
from pyscf import gto, dft, tdscf

# 1. 基态几何
benzene_gs = gto.M(
    atom='''
    C  0.000  1.396  0.000
    C  1.210  0.698  0.000
    C  1.210 -0.698  0.000
    C  0.000 -1.396  0.000
    C -1.210 -0.698  0.000
    C -1.210  0.698  0.000
    H  0.000  2.476  0.000
    H  2.147  1.239  0.000
    H  2.147 -1.239  0.000
    H  0.000 -2.476  0.000
    H -2.147 -1.239  0.000
    H -2.147  1.239  0.000
    ''',
    basis='cc-pvdz'
)

# 2. 基态 DFT
mf_gs = dft.RKS(benzene_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# 3. 吸收光谱
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()
abs_e = td_abs.e * 27.2114

print(f"吸收 S₁: {abs_e[0]:.2f} eV ({1240/abs_e[0]:.0f} nm)")

# 4. 激发态几何（近似：键长增加 2-3%）
benzene_ex = gto.M(
    atom='''
    C  0.000  1.430  0.000
    C  1.239  0.715  0.000
    C  1.239 -0.715  0.000
    C  0.000 -1.430  0.000
    C -1.239 -0.715  0.000
    C -1.239  0.715  0.000
    H  0.000  2.510  0.000
    H  2.190  1.255  0.000
    H  2.190 -1.255  0.000
    H  0.000 -2.510  0.000
    H -2.190 -1.255  0.000
    H -2.190  1.255  0.000
    ''',
    basis='cc-pvdz'
)

# 5. 激发态几何下的基态计算
mf_ex = dft.RKS(benzene_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# 6. 发射光谱
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()
em_e = td_em.e * 27.2114

print(f"发射 S₁→S₀: {em_e[0]:.2f} eV ({1240/em_e[0]:.0f} nm)")
print(f"Stokes 位移: {abs_e[0] - em_e[0]:.3f} eV")
```

## 激发态优化方法

### 方法 1: PySCF 内置优化（推荐）

```python
from pyscf.geomopt import geometric_solver

# TDA 激发态优化（更稳定）
td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
mol_ex = geometric_solver.optimize(td_s1, state=0)
```

### 方法 2: 使用外部优化器

```python
# 使用 PyBerny 优化器
from pyscf.geomopt import berny_solver

mol_ex = berny_solver.optimize(td_s1)
```

### 方法 3: 手动几何调整（快速近似）

激发态几何经验规则：
- **C-C 键**: 增长 2-4%
- **C-H 键**: 增长 1-2%
- **C=O 键**: 增长 3-5%

## 常见问题

### Q1: 激发态优化不收敛？

**解决方案:**
1. 使用 TDA 代替全 TDDFT（更稳定）
2. 增加初始猜测质量（用基态几何）
3. 降低收敛标准（`conv_tol=1e-4`）

```python
td_s1 = tdscf.TDA(mf_gs)
td_s1.conv_tol = 1e-4
```

### Q2: 振子强度为 0（禁阻跃迁）？

**原因:** 对称性禁阻

**检查:**
```python
# 查看对称性
print(f"分子对称性: {mol.symmetry}")

# 查看激发态对称性
for i, e in enumerate(td.e):
    print(f"S{i+1}: {e*27.2114:.2f} eV, f={td.oscillator_strength()[i]:.4f}")
```

### Q3: Stokes 位移过大/过小？

**正常范围:** 0.1-1.0 eV

**异常情况:**
- >1.0 eV: 检查几何优化是否正确
- <0.05 eV: 可能激发态几何变化很小

## 性能优化

### 1. 减少计算状态数

```python
# 只计算前 3 个态
td.nstates = 3  # 而不是 10
```

### 2. 使用密度拟合

```python
mf = dft.RKS(mol).density_fit(auxbasis='def2-universal-jkfit')
```

### 3. 粗网格（测试用）

```python
mf.grids.atom_grid = (50, 194)  # 默认 (75, 302)
```

## 结果验证

### 检查清单

- [ ] 基态 SCF 收敛
- [ ] TDDFT 收敛
- [ ] 激发态优化收敛
- [ ] Stokes 位移合理（0.1-1.0 eV）
- [ ] 振子强度非负
- [ ] 波长在合理范围（紫外/可见区）

### 对比实验数据

典型芳香分子的 Stokes 位移：
- 苯: ~0.1-0.3 eV
- 萘: ~0.2-0.4 eV
- 蒽: ~0.3-0.5 eV

## 可视化

```python
import matplotlib.pyplot as plt

# 高斯展宽
def gaussian(x, mu, sigma=0.15):
    return np.exp(-((x - mu)**2) / (2 * sigma**2))

energy_range = np.linspace(3, 8, 1000)

# 吸收光谱
abs_spectrum = sum(gaussian(energy_range, e) * f 
                   for e, f in zip(abs_e, td_abs.oscillator_strength()))

# 发射光谱
em_spectrum = sum(gaussian(energy_range, e) * f 
                  for e, f in zip(em_e, td_em.oscillator_strength()))

plt.plot(energy_range, abs_spectrum, 'b-', label='吸收')
plt.plot(energy_range, em_spectrum, 'r-', label='发射')
plt.xlabel('能量 (eV)')
plt.ylabel('强度')
plt.legend()
plt.show()
```

## 参考文献

1. **TDDFT 理论**
   - Casida, M. E. (1995). "Time-dependent density functional response theory for molecules"

2. **激发态优化**
   - Furche, F. & Ahlrichs, R. (2002). "Adiabatic time-dependent density functional methods"

3. **Stokes 位移**
   - Lakowicz, J. R. (2006). "Principles of Fluorescence Spectroscopy"

## 相关资源

- PySCF 文档: http://www.pyscf.org
- TDDFT 教程: `pyscf/examples/tdscf/`
- 几何优化: `pyscf/examples/geomopt/`
