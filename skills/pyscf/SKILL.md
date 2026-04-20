# PySCF Skill

## 简介

PySCF 是模块化的纯 Python 量子化学程序库，支持：
- **波函数方法**：HF、MP2、CCSD、CCSD(T)、CASSCF、CASCI
- **DFT**：所有标准泛函（LDA、GGA、meta-GGA、hybrids）
- **激发态**：TDDFT、TDA、CIS、SF-DFT
- **几何优化**：能量最小化、过渡态
- **势能面**：内坐标扫描、 relaxed scan
- **谱分析**：UV-Vis 吸收、发射、CD 光谱

## 核心概念

### 分子定义
```python
from pyscf import gto
mol = gto.M(
    atom='H 0 0 0; F 0 0 1.0',  # 内联坐标
    # atom=open('h2o.xyz').read(),   # 或从文件读取
    basis='cc-pvdz',                # 基组
    charge=0,                      # 电荷
    spin=0,                        # 自旋 (0=闭壳层, 1=开壳层)
    verbose=4                       # 输出级别
)
```

### 基组选择
| 任务 | 推荐基组 |
|------|---------|
| 快速筛选 | sto-3g, 3-21G |
| 平衡精度 | 6-31G*, cc-pVDZ |
| 高精度 | cc-pVTZ, cc-pVQZ, def2-TZVP |
| 重原子 | LANL2DZ, Stuttgart RSC |

### 常用泛函
```python
# GGA
'pbe', 'blyp', 'bp86', 'pbe0'  # PBE0 是 hybrid

# meta-GGA
'tkwwl', 'scan', 'm06l'

# hybrid
'b3lyp', 'pbe0', 'm06', 'm06-2x', 'wb97x', 'wb97x-d'

# dispersion
'b3lyp-d3', 'wb97x-d3bj'
```

### 运行模式
```python
mf = scf.RHF(mol)
mf.kernel()           # 正式运行（保存中间结果）
mf.run()             # 同上
mf.converged         # 检查收敛
mf.e_tot             # 总能量

# 渐进式
mf = scf.RHF(mol).newton()  # Newton-Raphson 加速
mf = scf.RHF(mol).density_fit()  # 密度拟合加速
```

---

## 1. 波函数方法

### 1.1 Hartree-Fock (HF)

**闭壳层 RHF：**
```python
from pyscf import scf
mol = gto.M(atom='H 0 0 0; F 0 0 1.0', basis='cc-pvdz')
mf = scf.RHF(mol).run()
print('E(RHF) =', mf.e_tot)
```

**开壳层 ROHF/UHF：**
```python
mol = gto.M(atom='O 0 0 0; O 0 0 1.2', spin=2, basis='cc-pvdz')  # O2 triplet
mf = scf.UHF(mol).run()          # Unrestricted
mf = scf.ROHF(mol).run()          # Restricted open-shell
```

**DHF（Dirac-Hartree-Fock，用于相对论）：**
```python
from pyscf import scf
mol = gto.M(atom='Hg 0 0 0', basis='cc-pvdz')
mf = scf.DHF(mol).run()
```

### 1.2 MP2（二阶微扰理论）
```python
from pyscf import mp, scf, gto

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# 闭壳层 MP2
mmp = mp.MP2(mf).run()
print('E(MP2) =', mmp.e_tot)
print('E(corr) =', mmp.e_corr)         # 相关能
print('T2 amplitudes shape:', mmp.t2.shape)

# 开壳层 UMP2
mmp = mp.UMP2(mf).run()
```

### 1.3 CCSD / CCSD(T)

**CCSD（耦合簇单双）：**
```python
from pyscf import cc, scf, gto

mol = gto.M(atom='N 0 0 0; N 0 0 1.1', basis='cc-pvdz')
mf = scf.RHF(mol).run()

mycc = cc.CCSD(mf).run()
print('E(CCSD) =', mycc.e_tot)
print('T1 amplitudes:', mycc.t1.shape)
print('T2 amplitudes:', mycc.t2.shape)

# 计算单粒子密度矩阵
dm = mycc.make_rdm1()
```

**CCSD(T)（考虑三次项）：**
```python
eccsd_t = cc.CCSD(mf).ccsd_t()  # 返回修正值 ΔE
e_total = mycc.e_tot + eccsd_t
```

**开壳层 UCCSD：**
```python
mol = gto.M(atom='N 0 0 0', spin=1, basis='cc-pvdz')
mf = scf.UHF(mol).run()
mycc = cc.UCCSD(mf).run()
```

### 1.4 CASSCF / CASCI

**CASSCF（完全活化空间 SCF）：**
```python
from pyscf import mcscf, scf, gto

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='cc-pvdz')
mf = scf.RHF(mol).run()

# CAS(nelec, ncas) = CASSCF
# 10 电子, 8 个轨道活化空间
mc = mcscf.CASSCF(mf, 8, 10).run()
print('E(CASSCF) =', mc.e_tot)
print('CAS space orbitals:', mc.mo_coeff.shape)

# 自定义活性空间（手动选择轨道）
mc = mcscf.CASSCF(mf, 8, 10)
mo = mcscf.sort_mo(mc, mf.mo_coeff, [15,16,17,20,21,22,25,26])
mc.run(mo)
```

**CASCI：**
```python
mc = mcscf.CASCI(mf, 8, 10).run()
```

**CASSCF + 微扰修正（CASPT2）：**
```python
mc = mcscf.CASSCF(mf, 8, 10).run()
mcpt2 = mcscf.CASPT2(mc).run()
print('E(CASPT2) =', mcpt2.e_tot)
```

**NEVPT2：**
```python
mcpt2 = mcscf.NEVPT2(mc).run()
```

---

## 2. DFT 方法

### 2.1 基本 DFT 计算
```python
from pyscf import dft

mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')

# RKS = Restricted Kohn-Sham (闭壳层)
mf = dft.RKS(mol, xc='pbe').run()
print('E(PBE) =', mf.e_tot)

mf = dft.RKS(mol, xc='b3lyp').run()
print('E(B3LYP) =', mf.e_tot)

# UKS = Unrestricted Kohn-Sham (开壳层)
mol = gto.M(atom='O 0 0 0', spin=2, basis='cc-pvdz')
mf = dft.UKS(mol, xc='pbe').run()
```

### 2.2 泛函选择指南
```python
xc_options = {
    # LDA
    'lda': 'svwn5',           # VWN5 Slater-VWN (RPA)
    
    # GGA
    'pbe': 'pbe',              # Perdew-Burke-Ernzerhof
    'blyp': 'blyp',            # Becke-Lee-Yang-Parr
    'bp86': 'bp86',            # B88 + P86
    
    # meta-GGA
    'scan': 'scan',            # Strongly Constrained
    
    # hybrid
    'b3lyp': 'b3lyp',         # 最常用
    'pbe0': 'pbe0',            # PBE0 (25% HF exchange)
    'm06': 'm06',              # Minnesota 2006
    'm06-2x': 'm06-2x',        # 54% HF exchange
    'wb97x': 'wb97x',          # range-separated
    'wb97x-d3': 'wb97x-d3bj',  # with dispersion
    
    # dispersion correction
    'pbe-d3': 'pbe,pbe-d3bj',  # DFT-D3 BJ
    'b3lyp-d3': 'b3lyp,b3lyp-d3bj',
}
```

### 2.3 range-separated 泛函
```python
# ωB97X-D: 长程 corrected，适合 CT 态和电荷转移
mf = dft.RKS(mol, xc='wb97x-d').run()

# CAM-B3LYP: 长程 correction 版本
mf = dft.RKS(mol, xc='camb3lyp').run()
```

### 2.4 轨道可视化 / DOS 分析
```python
# 轨道系数
mo_coeff = mf.mo_coeff
mo_energy = mf.mo_energy

# 绘制轨道
from pyscf import dft
dft.gen_grid.with_rounding = False
dft.numint.NR_3Done = dft.numint._dot_ao_ao

# Mulliken population
from pyscf import lo
occ = lo.orth_ao(mol, 'minao')
dm = mf.make_rdm1()
mulliken = lo.mulliken(mol, dm, occ)
```

---

## 3. 激发态：TDDFT / TDA

### 3.1 TDDFT（线性响应）
```python
from pyscf import tddft, dft, scf, gto

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# TDDFT 计算
td = tddft.TDDFT(mf).run()
print('激发态能量 (eV):', td.e * 27.2114)
print('振子强度:', td.f)
print('S1 能量 (eV):', td.e[0] * 27.2114)

# 计算自然跃迁轨道 (NTO)
td.analyze()
```

### 3.2 TDA（单激发近似，Tamm-Dancoff）
```python
# TDA 比 TDDFT 快，但对 double excitation 不太准确
td = tddft.TDA(mf).run(nstates=5)  # 默认 10 个态
td.analyze()
```

### 3.3 UV-Vis 光谱
```python
import numpy as np
import matplotlib.pyplot as plt

# 计算 50 个激发态
td = tddft.TDDFT(mf).run(nstates=50)

# 展宽光谱
e = td.e * 27.2114  # eV
f = td.f
sigma = 0.1  # eV

energies = np.linspace(0, 10, 500)
spectrum = np.zeros_like(energies)
for ei, fi in zip(e, f):
    spectrum += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))

plt.plot(energies, spectrum)
plt.xlabel('Energy (eV)')
plt.ylabel('Oscillator strength')
plt.show()
```

### 3.4 开壳层 TDDFT
```python
mol = gto.M(atom='O 0 0 0', spin=2, basis='cc-pvdz')
mf = dft.UKS(mol, xc='pbe').run()
td = tddft.TDDFT(mf).run(nstates=5)
```

### 3.5 SF-DFT（态选择 DFT，适合 diradical）
```python
from pyscf import tddft
mf = dft.RKS(mol, xc='b3lyp').run()
td = tddft.TDDFT(mf).run(nstates=3)
# 或使用 sf王爷
from pyscf import df
td_sf = df.DFTStateTransfer(mf, ['4.0', '5.0'])  # 指定能隙范围
```

---

## 4. 几何优化

### 4.1 能量最小化
```python
from pyscf import dft, geomopt

mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# 标准优化
opt = geomopt.GeometryOptimizer(mf)
mol_opt = opt.run()

print('优化后坐标:')
for i in range(mol_opt.natm):
    print(mol_opt.atom_coord(i))

# 直接用 RKS 对象
mol_opt = dft.RKS(mol, xc='b3lyp').geomopt().run()
```

### 4.2 过渡态优化
```python
from pyscf import dft, geomopt

mol = gto.M(atom='H 0 0 0; Cl 0 0 1.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()

# 选择内坐标或直接指定
opt = geomopt.GeometryOptimizer(mf, 'gen炒股')
# 使用二聚体法找过渡态
opt = geomopt.DimerOptimizer(mf).run()
```

### 4.3 频率计算（验证极小点/过渡态）
```python
from pyscf import hessian

# 计算 Hessian
mol = gto.M(atom='O 0 0 0; H 0 0 1.8; H 0 0.8 2.7', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()
hess = hessian.RHF(mf).run()

# 频率
freq = hess.kernel()
print('频率 (cm⁻¹):', freq)
```

---

## 5. 势能面（PES）扫描

### 5.1 Relaxed 扫描（逐步优化）
```python
from pyscf import dft, relaxscan
import numpy as np

mol = gto.M(atom='C 0 0 0; H 0 0 1.1', basis='cc-pvdz', unit='Angstrom')
mf = dft.RKS(mol, xc='pbe').run()

# 扫描 C-H 键长
scan_coords = []
scan_energies = []

for r in np.linspace(0.9, 1.3, 10):
    mol_temp = gto.M(
        atom='C 0 0 0; H 0 0 %.4f' % r,
        basis='cc-pvdz'
    )
    mf_temp = dft.RKS(mol_temp, xc='pbe').run(chkfile=False)
    scan_coords.append(r)
    scan_energies.append(mf_temp.e_tot)

print('扫描完成')
for c, e in zip(scan_coords, scan_energies):
    print(f'r={c:.3f} E={e:.6f}')
```

### 5.2 1D 内坐标扫描
```python
from pyscf import mcscf

# 更复杂的扫描可用 as_scan
mol = gto.M(atom='N 0 0 0; N 0 0 1.5', basis='cc-pvdz')
mf = scf.RHF(mol).run()
mc = mcscf.CASSCF(mf, 4, 4).run()

# 扫描轨迹记录
scan_result = mc.as_scanner()
```

---

## 6. 光谱分析

### 6.1 UV-Vis 吸收光谱
```python
import numpy as np
import matplotlib.pyplot as plt

# 从 TDDFT 结果生成光谱
def uv_vis_spectrum(td, nstates=50, sigma=0.05, emin=0, emax=10):
    """生成展宽的 UV-Vis 光谱"""
    e = td.e[:nstates] * 27.2114  # eV
    f = td.f[:nstates]
    
    energies = np.linspace(emin, emax, 1000)
    absorption = np.zeros_like(energies)
    
    for ei, fi in zip(e, f):
        if fi > 0:
            absorption += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))
    
    return energies, absorption

# 使用
mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='b3lyp').run()
td = tddft.TDDFT(mf).run(nstates=50)

e, spec = uv_vis_spectrum(td, sigma=0.05)
plt.plot(e, spec)
plt.xlabel('Photon energy (eV)')
plt.ylabel('Absorption (arb. units)')
plt.show()
```

### 6.2 发射光谱（从 S1 态）
```python
from pyscf import tddft, dft

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = dft.RKS(mol, xc='pbe').run()
td = tddft.TDDFT(mf).run(nstates=5)

# S1 发射 = S1 态 → S0 垂直跃迁
td.emit()  # 计算发射（实验功能）
```

### 6.3 CD 光谱（圆二色谱）
```python
from pyscf import tddft, gdft

mol = gto.M(atom='C 0 0 0; N 0 0 1.2', basis='cc-pvdz')
mf = gdft.RKS(mol, xc='pbe').run()
td = tddft.TDDFT(mf).run(nstates=20)

# 旋光强度 (rotatory strength)
td.analyze()  # 包含 rotatory strength
```

---

## 7. 并行化与加速

### 7.1 密度拟合 (DF)
```python
from pyscf import df, scf

mol = gto.M(atom='C 0 0 0; O 0 0 1.2', basis='def2-tzvp')
mf = scf.RHF(mol).density_fit(auxbasis='def2-tzvp-ri').run()
```

### 7.2 多线程 BLAS
```python
import os
os.environ['OMP_NUM_THREADS'] = '16'

# 或在运行时
from pyscf.lib import num_threads
num_threads(16)
```

### 7.3 CHK 文件（断点续算）
```python
mf = scf.RHF(mol)
mf.chkfile = 'hf.chk'
mf.run()

# 重新加载
mf = scf.RHF(mol)
mf.init_guess = 'chk'
mf.run(chkfile='hf.chk')
```

---

## 8. 常用基组速查

| 基组 | 描述 | 适用场景 |
|------|------|---------|
| sto-3g | 最小基组，最快 | 快速测试 |
| 3-21G | 分裂基组 | 小体系筛选 |
| 6-31G* | 中等精度，有极化 | 常规有机分子 |
| cc-pVDZ | 相关一致性基组 | MP2/CCSD |
| cc-pVTZ | 高精度 | 精细计算 |
| def2-TZVP | 赝势基组 | 重原子 |
| LANL2DZ | 有效核势势 | 过渡金属 |

---

## 9. 常见问题

**Q: SCF 不收敛？**
```python
# 方法 1: 更换初始猜测
mf.init_guess = 'atom'    # 原子密度
mf.init_guess = 'hcore'   # 核哈密顿量
mf.diis = True             # DIIS 加速（默认开启）

# 方法 2: Newton-Raphson
mf = scf.RHF(mol).newton().run()

# 方法 3: 更换 SCF 方法
mf = scf.DHF(mol).run()  # Dirac-HF（重原子）
```

**Q: 内存不够？**
```python
# 使用密度拟合减少内存
mf = scf.RHF(mol).density_fit().run()

# 使用 direct SCF
mf.direct_scf = True
```

**Q: 如何输出所有轨道能量？**
```python
mf = scf.RHF(mol).run()
print('HOMO:', mf.mo_energy[mf.mo_occ > 0].max())
print('LUMO:', mf.mo_energy[mf.mo_occ == 0].min())
```

---

## Tools 脚本索引

| 脚本 | 功能 |
|------|------|
| `tools/scf.py` | HF/DFT 基础框架 |
| `tools/tddft.py` | TDDFT/TDA 计算 |
| `tools/dft.py` | DFT 泛函选择 |
| `tools/mp2.py` | MP2 计算框架 |
| `tools/ccsd.py` | CCSD/CCSD(T) 框架 |
| `tools/cascf.py` | CASSCF/CASPT2 框架 |
| `tools/geometry.py` | 几何优化 |
| `tools/spectrum.py` | 光谱生成与可视化 |
| `tools/pes.py` | 势能面扫描 |
| `tools/analysis.py` | 波函数分析（Mulliken, NBO, DOS） |
