# PySCF - Python量子化学库

## Description
当用户询问关于PySCF库的使用方法、API、代码示例、模块结构、SCF/TDDFT功能实现、JAX集成、输入输出格式等PySCF相关编程问题时触发此skill。

## 触发关键词
- PySCF
- pyscf库
- PySCF代码
- PySCF SCF
- PySCF TDDFT
- PySCF教程
- PySCF文档
- PySCF API
- PySCF示例
- PySCF CASSCF
- PySCF DFT
- PySCF安装

## PySCF简介

PySCF (Python-based Simulations of Chemistry Framework) 是一个开源的Python量子化学计算软件包，特点：
- **Python原生**: 易读易写，快速开发
- **高性能**: 基于NumPy/SciPy + C/C++优化
- **模块化**: 易于集成和扩展
- **全面**: 涵盖分子和周期体系的各种方法

## 核心模块结构

```
pyscf/
├── gto/          # 基组、分子/晶胞定义
├── scf/          # SCF方法（HF, KS-DFT）
├── dft/          # DFT泛函和网格
├── cc/           # 耦合簇（CCSD, CCSD(T)）
├── mp/           # 微扰理论（MP2）
├── ci/           # 组态相互作用
├── mcscf/        # 多组态SCF（CASSCF）
├── fci/          # 全CI求解器
├── tdscf/        # 含时SCF（TDDFT）
├── ao2mo/        # AO到MO积分变换
├── df/           # 密度拟合
├── grad/         # 梯度计算
├── geomopt/      # 几何优化
├── solvent/      # 溶剂效应（PCM, COSMO）
├── pbc/          # 周期边界条件
└── lib/          # C/C++扩展库
```

## 安装

```bash
# 基础安装
pip install pyscf

# 完整安装（包含所有功能）
pip install pyscf-full

# 从源码安装
git clone https://github.com/pyscf/pyscf.git
cd pyscf
pip install -e .
```

## 快速入门

### 1. 分子定义

```python
from pyscf import gto

# 水分子
mol = gto.M(
    atom='O 0 0 0; H 0 1 0; H 0 0 1',
    basis='cc-pvdz',
    charge=0,
    spin=0,  # 2S = nα - nβ
    symmetry=False,  # 或 'd2h', 'c2v'等
)

# 对称性分子
mol_c2 = gto.M(
    atom='C 0 0 0.625; C 0 0 -0.625',
    basis='cc-pvdz',
    symmetry='d2h'
)

# 打印分子信息
print(mol)
```

### 2. Hartree-Fock计算

```python
from pyscf import scf

# RHF（限制性闭壳层HF）
mf = scf.RHF(mol)
energy = mf.kernel()

# UHF（非限制性开壳层HF）
mol_o2 = gto.M(
    atom='O 0 0 0; O 0 0 1.2',
    basis='cc-pvdz',
    spin=2  # 三重态
)
uhf = scf.UHF(mol_o2)
energy = uhf.kernel()

# ROHF（限制性开壳层HF）
rohf = scf.ROHF(mol_o2)
energy = rohf.kernel()

# 结果访问
print(f'能量: {mf.e_tot:.6f} Hartree')
print(f'HOMO能量: {mf.mo_energy[mf.mo_occ>0][-1]:.6f}')
print(f'LUMO能量: {mf.mo_energy[mf.mo_occ==0][0]:.6f}')
```

### 3. Kohn-Sham DFT

```python
from pyscf import dft

# RKS（限制性KS-DFT）
mf = dft.RKS(mol)
mf.xc = 'b3lyp'  # 泛函选择
energy = mf.kernel()

# 常见泛函
mf.xc = 'pbe0'        # PBE0杂化
mf.xc = 'wb97x-d'     # ωB97X-D（范围分离+色散）
mf.xc = 'cam-b3lyp'   # CAM-B3LYP（范围分离）
mf.xc = 'scan'        # SCAN meta-GGA

# 自定义泛函
mf.xc = '.2*HF + .08*LDA + .72*B88, .81*LYP + .19*VWN'  # B3LYP配方

# 网格设置
mf.grids.atom_grid = (100, 770)  # 角度×径向网格
mf.grids.prune = None  # 不修剪

# 色散校正
mf.xc = 'wb97m_v'
mf.nlc = 'vv10'
mf.nlcgrids.atom_grid = (50, 194)
```

### 4. LR-TDDFT激发态

```python
from pyscf import tdscf

# TDDFT（线性响应）
td = tdscf.TDDFT(mf)
td.nstates = 6  # 计算激发态数
td.kernel()

# Tamm-Dancoff近似（TDA）
td_tda = tdscf.TDA(mf)
td_tda.nstates = 4
td_tda.kernel()

# 提取激发能和振子强度
for i, (e, f) in enumerate(zip(td.e, td.oscillator_strength)):
    print(f'态 {i+1}: {e*27.2114:.2f} eV, f = {f:.3f}')

# 自然跃迁轨道（NTO）分析
weights, nto = td.get_nto(state=0)
print(f'主导跃迁权重: {weights.max():.3f}')
```

### 5. 后HF方法

```python
from pyscf import mp, cc

# MP2
mp2 = mp.MP2(mf)
emp2 = mp2.kernel()[0]

# CCSD
ccsd = cc.CCSD(mf)
eccsd, t1, t2 = ccsd.kernel()

# CCSD(T)
e_ccsdt = ccsd.ccsd_t()

# DLPNO-CCSD（近似方法）
from pyscf.cc import ccsd
ccsd_dlpno = ccsd.DLPNOCCSD(mf, auxbasis='cc-pvtz-jkfit')
```

### 6. CASSCF

```python
from pyscf import mcscf

# CAS(6,8)：6轨道8电子
cas = mcscf.CASSCF(mf, 6, 8)
e_cas, ci_vec, mo, mo_occ = cas.kernel()

# 带密度拟合的CASSCF
cas_df = mcscf.DFCASSCF(mf, 6, 8, auxbasis='cc-pvtz-fit')
cas_df.kernel()

# 态平均CASSCF
cas_sa = mcscf.CASSCF(mf, 6, 8)
cas_sa.state_average_([0.5, 0.3, 0.2])  # 3个态

# NEVPT2
from pyscf import mrpt
e_nevpt2 = mrpt.NEVPT2(cas).kernel()
```

### 7. 几何优化

```python
from pyscf.geomopt.geometric_solver import optimize as geometric_opt

# HF/DFT优化
opt_mol = geometric_opt(mf)

# CASSCF优化
opt_mol_cas = geometric_opt(cas)

# 激发态优化
td.state_specific_(0)  # 第1激发态
opt_mol_ex = geometric_opt(td)
```

### 8. 溶剂效应

```python
# PCM（极化连续介质模型）
mf_pcm = mol.RHF().ddPCM()
mf_pcm.kernel()
mf_pcm.with_solvent.eps = 78.4  # 水

# COSMO
mf_cosmo = mol.RHF().ddCOSMO()
mf_cosmo.kernel()

# 点电荷（QM/MM）
from pyscf import qmmm
coords = [[0., 0., 2.], [0., 0., -2.]]
charges = [-0.5, -0.5]
mf_qmmm = qmmm.mm_charge(mf, coords, charges)
```

### 9. 密度拟合

```python
# 自动密度拟合
mf_df = mf.density_fit(auxbasis='def2-universal-jfit')

# 手动密度拟合
from pyscf import df
mf_df = df.density_fit(scf.RHF(mol), auxbasis='cc-pvtz-jkfit')
```

### 10. 周期边界条件

```python
from pyscf.pbc import gto as pbcgto, dft as pbcdft, scf as pbcscf

# 晶胞定义
cell = pbcgto.M(
    atom='''
    C 0. 0. 0.
    C 0.89 0.89 0.89
    ''',
    basis='gth-szv',
    pseudo='gth-pade',
    a=[[3.5, 0, 0], [0, 3.5, 0], [0, 0, 3.5]]
)

# k点采样
kpts = cell.make_kpts([2, 2, 2])

# PBC-DFT
mf_pbc = pbcdft.KRKS(cell, kpts).density_fit(auxbasis='weigend')
mf_pbc.xc = 'pbe'
mf_pbc.kernel()
```

## 高级功能

### 1. 自定义泛函

```python
from pyscf.dft import gen_grid

# 完全自定义
def my_xc(rho, spin=0, relativity=0, deriv=1, verbose=None):
    # rho: (N, 4) array [density, grad_x, grad_y, grad_z]
    ex, vrho, vgamma = ...  # 计算交换能和势
    ec, vrho, vgamma = ...  # 计算相关能和势
    return (ex + ec), (vrho[0] + vrho[0]), (vgamma[0] + vgamma[0])

mf._numint._xc = my_xc
mf.kernel()
```

### 2. 积分访问

```python
import numpy as np
from pyscf import ao2mo

# 1电子积分（AO基）
T = mol.intor('int1e_kin')  # 动能
V_nuc = mol.intor('int1e_nuc')  # 核吸引
H_core = T + V_nuc

# 2电子积分（AO基）
eri_ao = mol.intor('int2e_sph', aosym=4)  # 4重对称

# MO基积分
eri_mo = ao2mo.incore.full(eri_ao, mf.mo_coeff)
```

### 3. 布居分析

```python
# Mulliken布居
from pyscf import lo
pop = mf.mulliken_pop(mol, mf.make_rdm1())

# Lowdin布居
pop_lowdin = lo.vec_lowdin(mf.mo_coeff, mf.get_ovlp())
dm = mf.make_rdm1()
pop = mf.pop_analysis(mol, dm, mf.mo_coeff)

# NBO分析
from pyscf.tools import cubegen
cubegen.orbital(mol, 'h2o_homo.cube', mf.mo_coeff[:, mf.mo_occ>0][:, -1])
```

### 4. 梯度和频率

```python
# 能量梯度
g = mf.nuc_grad_method()
grad = g.kernel()

# 振动频率
from pyscf.hessian import rhf as rhf_hessian
h = rhf_hessian.Hessian(mf)
freq = h.kernel()
```

### 5. 与JAX集成（自动微分）

```python
# PySCF支持JAX用于自动微分
import jax
import jax.numpy as jnp

from pyscf import gto, scf
from pyscf.jax import scf as jax_scf

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol)

# 转换为JAX版本
mf_jax = jax_scf.RHF(mol)

# 自动微分梯度
grad_energy = jax.grad(lambda mol: mf_jax.kernel()[0])
```

## 实用工作流

### 工作流1：激发态计算

```python
from pyscf import gto, dft, tdscf

# 1. 分子和基态DFT
mol = gto.M(atom='C 0 0 0; O 0 0 1.1', basis='cc-pvdz')
mf = dft.RKS(mol)
mf.xc = 'cam-b3lyp'
mf.kernel()

# 2. LR-TDDFT激发态
td = tdscf.TDDFT(mf)
td.nstates = 10
td.kernel()

# 3. 分析激发态
print('\n激发态分析:')
for i in range(td.nstates):
    print(f'S{i+1}: {td.e[i]*27.2114:.2f} eV, '
          f'λ = {1240/(td.e[i]*27.2114):.1f} nm, '
          f'f = {td.oscillator_strength[i]:.3f}')

# 4. NTO分析
weights, nto = td.get_nto(state=0)
print(f'\n第1激发态NTO分析:')
print(f'主导跃迁权重: {weights.max():.3f}')
```

### 工作流2：几何优化+频率

```python
from pyscf import gto, dft, geomopt

# 1. 初始结构
mol = gto.M(atom='C 0 0 0; O 0 0 1.15', basis='def2-svp')

# 2. DFT优化
mf = dft.RKS(mol)
mf.xc = 'wb97x-d'
opt_mol = geomopt.geometric_solver.optimize(mf)

# 3. 频率计算
mf_opt = dft.RKS(opt_mol)
mf_opt.xc = 'wb97x-d'
mf_opt.kernel()

# 4. 频率分析
from pyscf.hessian import rks as rks_hessian
hess = rks_hessian.Hessian(mf_opt)
freq = hess.kernel()

# 检查虚频
if freq.min() < -10:
    print('警告：存在虚频，可能不是极小点')
```

### 工作流3：多参考计算（CASSCF）

```python
from pyscf import gto, scf, mcscf, mrpt

# 1. 分子和初始SCF
mol = gto.M(atom='C 0 0 0; C 0 0 1.2', basis='cc-pvdz')
mf = scf.RHF(mol)
mf.kernel()

# 2. 活性空间选择（6轨道8电子）
cas = mcscf.CASSCF(mf, 6, 8)

# 3. 自然轨道分析
no_occ, no_coeff = mcscf.cas_natorb(cas)
cas = mcscf.CASSCF(mf, 6, 8)
cas.mo_coeff = no_coeff

# 4. CASSCF计算
e_cas, ci_vec, mo, mo_occ = cas.kernel()

# 5. NEVPT2动态相关
e_nevpt2 = mrpt.NEVPT2(cas).kernel()

print(f'CASSCF能量: {e_cas:.6f} Hartree')
print(f'NEVPT2能量: {e_nevpt2:.6f} Hartree')
print(f'总能量: {e_cas + e_nevpt2:.6f} Hartree')
```

## 性能优化

### 1. 并行计算

```python
# 多核并行
import os
os.environ['OMP_NUM_THREADS'] = '8'

# MPI并行（通过mpi4py）
# mpirun -np 4 python script.py
```

### 2. 内存优化

```python
# 密度拟合减少内存
mf = mf.density_fit()

# 外存积分存储
mf.direct_scf = False  # 保存积分到磁盘
mf.chkfile = 'chkfile.h5'
```

### 3. 基组选择

```python
# 快速预计算
mol = gto.M(atom='...', basis='sto-3g')

# 高精度计算
mol = gto.M(atom='...', basis='def2-tzvp')

# 拟合基组
auxbasis = 'cc-pvdz-jkfit'
```

## 输入输出格式

### 常用基组

| 体系类型 | 基组 | 说明 |
|---------|------|------|
| 快速测试 | sto-3g, 3-21G | 最小基 |
| 常规计算 | cc-pVDZ, def2-SVP | 双ζ |
| 高精度 | cc-pVTZ, def2-TZVP | 三ζ |
| 关键精度 | cc-pVQZ, def2-QZVP | 四ζ |

### 输出解读

```python
# SCF收敛信息
print(mf.converged)  # True/False
print(mf.scf_summary)  # SCF总结

# 能量分解
print(f'总能量: {mf.e_tot:.6f}')
print(f'1电子能: {mf.e_1e:.6f}')
print(f'2电子能: {mf.e_2e:.6f}')
```

## 常见问题

### Q1: SCF不收敛
```python
# 尝试不同初始猜测
mf.init_guess = 'huckel'
mf.init_guess = '1e'

# 使用DIIS
mf.diis_start_cycle = 3

# 水平移动
mf.level_shift = 0.5

# 直接反演迭代子空间
mf.damp_factor = 0.2
```

### Q2: 内存不足
```python
# 密度拟合
mf = mf.density_fit()

# 外存计算
mf.direct_scf = False
```

### Q3: 激发态失败
```python
# 先用TDA
td_tda = tdscf.TDA(mf)

# 检查基态质量
print(f'能隙: {(mf.mo_energy[mf.mo_occ==0][0] - mf.mo_energy[mf.mo_occ>0][-1])*27.2114:.2f} eV')
```

## 扩展与集成

### GradDFT
PySCF可与GradDFT集成用于梯度计算和泛函优化：
```python
# 示例：计算泛函梯度
from grad_dft import functional_gradient
grad = functional_gradient(mf)
```

### 与机器学习结合
```python
# 使用PySCF生成训练数据
# 1. 构建分子数据库
# 2. 计算DFT能量/力
# 3. 训练ML模型（如GNN）
```

## References
详见 `references/` 目录：
- PySCF官方文档
- API参考
- 示例代码
- 性能优化指南

## Related Skills
- `dft-theory`: DFT理论基础
