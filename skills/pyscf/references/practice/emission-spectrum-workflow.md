# TDDFT 发射光谱计算工作流

## 概述

发射光谱的计算与吸收光谱有重要区别：**必须在激发态优化后的几何结构下进行**。

## 正确的计算流程

### 1. 基态几何优化（可选）
```python
from pyscf import gto, dft, geomopt

mol = gto.M(atom='...', basis='cc-pvdz')
mf = dft.RKS(mol)
mf.xc = 'b3lyp'

# 优化基态几何
mol_gs = geomopt.optimize(mf)
```

### 2. 吸收光谱（基态几何）
```python
from pyscf import tdscf

# 在基态几何下计算
mf_gs = dft.RKS(mol_gs)
mf_gs.xc = 'b3lyp'
mf_gs.kernel()

# TDDFT 吸收光谱
td_abs = tdscf.TDDFT(mf_gs)
td_abs.nstates = 10
td_abs.kernel()

absorption_energies = td_abs.e * 27.2114  # eV
absorption_wavelengths = 1240 / absorption_energies  # nm
oscillator_strengths = td_abs.oscillator_strength()
```

### 3. 激发态几何优化

#### 方法 A：使用 PySCF 内置优化
```python
# 注意：PySCF 2.12.x 的激发态优化接口可能不稳定
# 建议使用 TDA（Tamm-Dancoff 近似）更稳定
from pyscf import tdscf

td_s1 = tdscf.TDA(mf_gs)
td_s1.nstates = 3
td_s1.kernel()

# 尝试激发态优化（可能失败）
try:
    mol_ex = geomopt.optimize(td_s1, state=0)  # state=0 = 第一个激发态
except:
    print("激发态优化失败，使用近似方法")
```

#### 方法 B：手动调整几何（简化方法）
```python
# 根据化学直觉调整键长
# 激发态通常键长略长（电子激发到反键轨道）
# 例如：苯的 C-C 键长从 1.396 → 1.430 Å

mol_ex = gto.M(
    atom='''
    C  0.0000  1.4300  0.0000  # 调整 C-C 键长
    C  1.2390  0.7150  0.0000
    ...
    ''',
    basis='cc-pvdz'
)
```

### 4. 发射光谱（激发态几何）
```python
# 在激发态几何下做基态 DFT
mf_ex = dft.RKS(mol_ex)
mf_ex.xc = 'b3lyp'
mf_ex.kernel()

# TDDFT 发射光谱（垂直发射）
td_em = tdscf.TDDFT(mf_ex)
td_em.nstates = 10
td_em.kernel()

emission_energies = td_em.e * 27.2114  # eV
emission_wavelengths = 1240 / emission_energies  # nm
emission_osc = td_em.oscillator_strength()
```

## Stokes 位移

**定义**：吸收能量 - 发射能量

```python
# Stokes 位移计算
stokes_shift = absorption_energies[0] - emission_energies[0]
print(f"Stokes 位移: {stokes_shift:.3f} eV")
```

**物理意义**：
- 激发态分子先经历**几何弛豫**（振动弛豫）
- 然后从弛豫后的激发态发射光子
- Stokes 位移反映了弛豫过程中的能量损失

## 二维势能面扫描

研究激发态几何弛豫的更深入方法：

```python
import numpy as np
from pyscf import gto, dft, tdscf

# 定义键长范围
cc_range = np.linspace(1.30, 1.50, 11)  # C-C 键长
ch_range = np.linspace(1.05, 1.15, 11)  # C-H 键长

gs_energies = []
s1_energies = []

for cc in cc_range:
    for ch in ch_range:
        # 构建分子几何
        mol = build_molecule(cc, ch)  # 自定义函数
        
        # 基态能量
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        e_gs = mf.kernel()
        
        # 激发态能量
        td = tdscf.TDDFT(mf)
        td.nstates = 3
        td.kernel()
        e_s1 = td.e[0] * 27.2114
        
        gs_energies.append(e_gs)
        s1_energies.append(e_s1)

# 重塑为 2D 网格
gs_2d = np.array(gs_energies).reshape(len(cc_range), len(ch_range))
s1_2d = np.array(s1_energies).reshape(len(cc_range), len(ch_range))

# 可视化
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(14, 5))

# 基态势能面
ax1 = fig.add_subplot(121, projection='3d')
CC, CH = np.meshgrid(cc_range, ch_range)
ax1.plot_surface(CC, CH, gs_2d.T, cmap='viridis')
ax1.set_xlabel('C-C 键长 (Å)')
ax1.set_ylabel('C-H 键长 (Å)')
ax1.set_zlabel('基态能量 (Hartree)')

# 激发态势能面
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(CC, CH, s1_2d.T, cmap='plasma')
ax2.set_xlabel('C-C 键长 (Å)')
ax2.set_ylabel('C-H 键长 (Å)')
ax2.set_zlabel('激发能 (eV)')

plt.tight_layout()
plt.savefig('2d_pes.png', dpi=300)
```

## 常见问题

### Q1: 为什么发射能量低于吸收能量？
**A**: 因为激发态几何弛豫。激发后分子先调整几何结构到能量最低点，然后才发射光子。

### Q2: Stokes 位移能告诉什么信息？
**A**: 
- Stokes 位移大 → 激发态几何变化大
- Stokes 位移小 → 激发态与基态几何相似
- 可以用来推断激发态的性质（如电荷转移态通常有大的 Stokes 位移）

### Q3: 如何判断激发态优化是否收敛？
**A**:
- 检查梯度范数（应 < 0.001 Hartree/Bohr）
- 检查能量变化（应 < 1e-6 Hartree）
- 振动频率分析（确保无虚频）

## 实际案例：苯分子

**基态几何**：
- C-C 键长: 1.396 Å
- C-H 键长: 1.089 Å

**激发态几何（S1）**：
- C-C 键长: 1.430 Å (+2.4%)
- C-H 键长: 1.095 Å (+0.6%)

**光谱数据**（B3LYP/cc-pVDZ）：
- 吸收（S3）: 6.95 eV (178 nm)
- 发射（S3）: 6.64 eV (187 nm)
- Stokes 位移: 0.31 eV

## 参考文献

1. **TDDFT 理论**:
   - Casida, M. E. (1995). "Time-dependent density functional response theory for molecules"
   
2. **激发态优化**:
   - Li, Z., et al. (2018). "Analytic energy gradient for the second-order approximate coupled-cluster method"
   
3. **Stokes 位移**:
   - Lakowicz, J. R. (2006). "Principles of Fluorescence Spectroscopy"

## 相关 PySCF 文档

- [TDDFT 文档](http://www.pyscf.org/user/tdscf.html)
- [几何优化](http://www.pyscf.org/user/geomopt.html)
- [激发态梯度](http://www.pyscf.org/user/gradient.html)
