# PySCF API详细参考

## 模块导入

```python
from pyscf import gto      # 基组和分子定义
from pyscf import scf      # SCF方法
from pyscf import dft      # DFT方法
from pyscf import mp       # 微扰理论（MP2等）
from pyscf import cc       # 耦合簇（CCSD等）
from pyscf import mcscf    # 多组态SCF
from pyscf import fci      # 全CI
from pyscf import tdscf    # 含时SCF（TDDFT）
from pyscf import ao2mo    # AO到MO积分变换
from pyscf import df       # 密度拟合
from pyscf import grad     # 能量梯度
from pyscf import geomopt  # 几何优化
from pyscf import solvent  # 溶剂效应
from pyscf.pbc import gto as pbcgto  # 周期体系GTO
from pyscf.pbc import scf as pbcscf  # 周期体系SCF
from pyscf.pbc import dft as pbcdft  # 周期体系DFT
```

## gto模块 - 分子和基组

### M类 - 分子/晶胞对象

#### 构造函数
```python
mol = gto.M(
    atom=None,              # 原子坐标
    basis=None,             # 基组
    spin=0,                 # 自旋多重度2S
    charge=0,                # 总电荷
    symmetry=False,         # 对称性
    unit='Ang',             # 坐标单位（'Ang'或'Bohr'）
    parse_arg=False,        # 是否解析参数
    verbose=4,              # 输出级别
)
```

#### 常用属性
```python
mol.atom          # 原子列表
mol.natm          # 原子数
mol.nel           # 电子数
mol.nao           # AO数
mol.nelectron     # 电子数
mol.spin          # 自旋2S
mol.charge        # 电荷
mol.basis         # 基组名
mol.symmetry      # 对称性
mol.irrep_id      # 不可约表示ID
mol.symm_orb      # 对称轨道
mol.cart          # 是否使用笛卡尔高斯
```

#### 常用方法

**分子几何**
```python
# 内坐标格式
mol = gto.M(
    atom='''
    O  0.0  0.0  0.0
    H  0.0  0.76 0.59
    H  0.0 -0.76 0.59
    ''',
    basis='cc-pvdz'
)

# Z矩阵格式
mol = gto.M(
    atom='''
    O
    H 1 0.96
    H 1 0.96 2 104.5
    ''',
    basis='cc-pvdz'
)

# 从XYZ文件读取
mol = gto.M(atom='h2o.xyz', basis='cc-pvdz')
```

**基组设置**
```python
# 使用内置基组
mol = gto.M(atom='...', basis='cc-pvdz')

# 自定义基组
mol = gto.M(
    atom='...',
    basis={
        'H': gto.parse('''
    H    S
      3.42525091              0.15432897
      0.62391373              0.53532814
      0.16885540              0.44463454
    '''),
        'O': gto.parse('...')
    }
)

# 混合基组
mol = gto.M(
    atom='O 0 0 0; H 0 1 0; H 0 0 1',
    basis={'O': 'cc-pvdz', 'H': 'cc-pvsz'}
)
```

**对称性**
```python
# 启用对称性
mol = gto.M(atom='...', symmetry=True)

# 指定点群
mol = gto.M(atom='...', symmetry='d2h')

# 可用点群：'d2h', 'c2v', 'c2h', 'd2', 'ci', 'cs', 'c1'
```

**积分计算**
```python
# 重叠矩阵
S = mol.intor('int1e_ovlp')

# 动能矩阵
T = mol.intor('int1e_kin')

# 核吸引矩阵
V = mol.intor('int1e_nuc')

# 核-核排斥能
Vnn = mol.energy_nuc()

# 双电子积分（8重对称）
eri = mol.intor('int2e')

# 双电子积分（4重对称）
eri_4 = mol.intor('int2e_sph', aosym=4)
```

**晶胞（PBC）**
```python
from pyscf.pbc import gto as pbcgto

cell = pbcgto.M(
    atom='...',
    basis='gth-szv',
    a=[[3.5, 0, 0], [0, 3.5, 0], [0, 0, 3.5]],  # 晶格向量
    pseudo='gth-pade',  # 赝势
    ke_cutoff=100,       # 动能截断
)

# k点采样
kpts = cell.make_kpts([2, 2, 2])
```

## scf模块 - SCF方法

### SCF基类

#### 常用类
```python
from pyscf import scf

# 闭壳层
mf = scf.RHF(mol)  # 限制性HF
mf = scf.UHF(mol)  # 非限制性HF
mf = scf.ROHF(mol) # 限制性开壳层HF

# 自由基
mf = scf.RKS(mol)  # 限制性KS-DFT
mf = scf.UKS(mol)  # 非限制性KS-DFT
```

#### 核心属性
```python
mf.mol           # 分子对象
mf.mo_coeff      # MO系数 (nao, nmo)
mf.mo_energy     # MO轨道能量
mf.mo_occ        # MO占据数
mf.e_tot         # 总能量
mf.converged     # 是否收敛
mf.scf_summary   # SCF总结字典
mf.chkfile       # 检查点文件
```

#### 核心方法

**运行SCF**
```python
# 标准运行
mf.kernel()
mf.kernel(dm0=None)  # 指定初始密度矩阵

# 返回能量和收敛标志
e_conv = mf.kernel()
print(f'能量: {e_conv:.6f}')
print(f'收敛: {mf.converged}')
```

**密度矩阵**
```python
# 总密度矩阵
dm = mf.make_rdm1()

# α和β密度矩阵（开壳层）
dma, dmb = mf.make_rdm1()
```

**Fock矩阵**
```python
# 有效Fock矩阵
fock = mf.get_fock()

# AO基Fock矩阵
fock_ao = mf.get_fock(dm=mf.make_rdm1())
```

**轨道分析**
```python
# HOMO能量
homo_idx = np.where(mf.mo_occ > 0)[0][-1]
homo_energy = mf.mo_energy[homo_idx]

# LUMO能量
lumo_idx = np.where(mf.mo_occ == 0)[0][0]
lumo_energy = mf.mo_energy[lumo_idx]

# 能隙
gap = lumo_energy - homo_energy
```

#### SCF控制

**初始猜测**
```python
# 核哈密顿量
mf.init_guess = '1e'

# Hückel方法
mf.init_guess = 'huckel'

# 从检查点读取
mf.init_guess = 'chkfile'
mf.chkfile = 'previous.chk'

# 从密度矩阵开始
dm_guess = mf.get_init_guess()
mf.kernel(dm0=dm_guess)
```

**收敛加速**
```python
# DIIS
mf.diis_start_cycle = 3  # DIIS开始周期
mf.diis = True

# 阻尼
mf.damp_factor = 0.2  # 阻尼因子

# 水平位移
mf.level_shift = 0.5  # 虚轨道位移

# 最大迭代
mf.max_cycle = 100

# 收敛阈值
mf.conv_tol = 1e-8  # 能量阈值
mf.conv_tol_grad = 1e-5  # 梯度阈值
```

**保存/读取**
```python
# 保存到检查点
mf.chkfile = 'calc.chk'
mf.dump_chk()  # 手动保存

# 从检查点读取
mf = scf.RHF(mol)
mf.chkfile = 'calc.chk'
mf.kernel()  # 自动读取

# 仅读取密度
mf_chk = scf.RHF(mol)
dm = mf_chk.from_chk('calc.chk')
```

## dft模块 - DFT方法

### RKS/UKS类

#### 基本用法
```python
from pyscf import dft

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.kernel()
```

#### 泛函设置

**内置泛函**
```python
# LDA
mf.xc = 'svwn'
mf.xc = 'vwn'

# GGA
mf.xc = 'pbe'
mf.xc = 'blyp'
mf.xc = 'bop'
mf.xc = 'bp86'

# meta-GGA
mf.xc = 'tpss'
mf.xc = 'm06-l'
mf.xc = 'scan'

# 杂化泛函
mf.xc = 'b3lyp'
mf.xc = 'pbe0'
mf.xc = 'wb97x-d'
mf.xc = 'm06-2x'

# 范围分离泛函
mf.xc = 'cam-b3lyp'
mf.xc = 'wb97x'
mf.xc = 'lc-wpbe'
```

**自定义泛函**
```python
# 方程式定义
mf.xc = '.2*HF + .08*LDA + .72*B88, .81*LYP + .19*VWN'

# 完全自定义函数
def my_xc(rho, spin=0, deriv=1, **kwargs):
    # rho: (N, 4) array [density, grad_x, grad_y, grad_z]
    # 返回: (ex+vc, vrho, vgamma)
    e_xc = ...  # 交换相关能密度
    vrho = ...  # 密度势
    vgamma = ...  # 梯度势
    return e_xc, vrho, vgamma

mf._numint._xc = my_xc
```

#### 网格设置

**原子网格**
```python
# (角度网格点数, 径向网格点数)
mf.grids.atom_grid = (75, 302)  # 标准
mf.grids.atom_grid = (99, 590)  # 精细
mf.grids.atom_grid = (250, 974) # 超精细

# 是否修剪
mf.grids.prune = None  # 不修剪
mf.grids.prune = dft.gen_grid.sg1_prune  # SG1修剪
```

**密度拟合网格**
```python
# 用于非局部色散
mf.nlc = 'vv10'  # 或 'rvv10'
mf.nlcgrids.atom_grid = (50, 194)
```

#### 分子积分

**交换相关势**
```python
# 交换相关势矩阵
vxc = mf.get_veff(mol, mf.make_rdm1())

# 分解
vx = mf.get_j(mol, mf.make_rdm1())  # 交换
vc = mf.get_k(mol, mf.make_rdm1())  # 相关
```

## tdscf模块 - 含时SCF

### TDDFT类

#### 基本用法
```python
from pyscf import tdscf

# TDDFT
td = tdscf.TDDFT(mf)
td.nstates = 6
td.kernel()

# TDA近似
td_tda = tdscf.TDA(mf)
td_tda.nstates = 4
td_tda.kernel()
```

#### 核心属性
```python
td.nstates         # 激发态数
td.e              # 激发能 (Hartree)
td.oscillator_strength  # 振子强度
td.vectors        # 跃迁向量
td.mf             # 基态SCF对象
```

#### 提取激发态信息

```python
# 激发能（eV）
for i, e in enumerate(td.e):
    print(f'态{i+1}: {e*27.2114:.2f} eV')

# 波长（nm）
for i, e in enumerate(td.e):
    wavelength = 1240/(e*27.2114)
    print(f'态{i+1}: {wavelength:.1f} nm')

# 振子强度
for i, f in enumerate(td.oscillator_strength):
    print(f'态{i+1}: f = {f:.3f}')

# 跃迁性质
for i in range(td.nstates):
    xy = td.transition_dipole[i]  # 跃迁偶极矩
    print(f'态{i+1}: μ = {xy}')
```

#### NTO分析

```python
# 自然跃迁轨道
weights, nto = td.get_nto(state=0)

# 空穴轨道（初始态）
hole = nto[0]  # (nao, nocc)

# 电子轨道（终态）
electron = nto[1]  # (nao, nvirt)

# 主导跃迁权重
dominant_weight = weights.max()
print(f'主导权重: {dominant_weight:.3f}')
```

#### 态平均

```python
# 态平均TDDFT
td = tdscf.TDDFT(mf)
td.state_average = True
td.nstates = 5
td.kernel()
```

## mcscf模块 - 多组态SCF

### CASSCF类

#### 基本用法
```python
from pyscf import mcscf

# CAS(nele, norb): 电子数, 轨道数
cas = mcscf.CASSCF(mf, 6, 8)
e_cas, ci_vec, mo, mo_occ = cas.kernel()
```

#### 核心属性
```python
cas.ncas           # 活性轨道数
cas.nelecas        # 活性电子数
cas.mo_coeff       # 优化后的MO系数
cas.ci             # CI向量
cas.e_cas          # CASSCF能量
cas.frozen         # 冻结轨道数
```

#### 初始轨道选择

```python
# 自动选择（HOMO/LUMO附近）
cas = mcscf.CASSCF(mf, 6, 8)

# 手动指定
cas = mcscf.CASSCF(mf, 6, 8)
cas.mo_coeff = mf.mo_coeff  # 从SCF开始

# 自然轨道（推荐）
no_coeff, no_occ = mcscf.cas_natorb(cas)
cas.mo_coeff = no_coeff
```

#### 态平均

```python
# 态平均CASSCF
cas = mcscf.CASSCF(mf, 6, 8)
cas.state_average([0.5, 0.3, 0.2])  # 3个态，权重0.5/0.3/0.2
cas.kernel()
```

#### 激发态

```python
# 态特定CASSCF
cas = mcscf.CASSCF(mf, 6, 8)
cas.state_specific(2)  # 第3激发态
cas.kernel()
```

#### 密度拟合CASSCF

```python
# 密度拟合加速
cas_df = mcscf.DFCASSCF(mf, 6, 8, auxbasis='cc-pvtz-jkfit')
cas_df.kernel()
```

## grad模块 - 能量梯度

### 梯度方法

#### HF/DFT梯度

```python
from pyscf import grad

# 获取梯度方法
g = mf.nuc_grad_method()
grad = g.kernel()  # (natm, 3)

# 打印梯度
print('能量梯度:')
for i in range(mol.natm):
    atom = mol.atom_symbol(i)
    print(f'{atom}: {grad[i]}')
```

#### CASSCF梯度

```python
# CASSCF梯度
g = cas.nuc_grad_method()
grad = g.kernel()
```

#### 梯度优化

```python
# 使用梯度进行几何优化
scanner = mf.as_scanner()
optimizer = scanner.nuc_grad_method().optimizer()
optimized_mol = optimizer.kernel()
```

## ao2mo模块 - 积分变换

### AO到MO变换

```python
from pyscf import ao2mo

# 完全变换
eri_mo = ao2mo.incore.full(eri_ao, mo_coeff)

# 部分变换（节省内存）
eri_mo = ao2mo.incore.general(eri_ao, [mo_occ, mo_occ, mo_occ, mo_occ])

# 外存变换（大分子）
ao2mo.outcore.full(mol, mo_coeff, 'eri_mo.h5')
```

## df模块 - 密度拟合

### 密度拟合SCF

```python
from pyscf import df

# 自动密度拟合
mf_df = mf.density_fit(auxbasis='def2-universal-jfit')

# 手动密度拟合
mf_df = df.density_fit(scf.RHF(mol), auxbasis='cc-pvtz-jkfit')

# 常用拟合基组
auxbasis = 'def2-universal-jfit'  # 通用拟合基
auxbasis = 'cc-pvtz-jkfit'       # JK拟合基
auxbasis = 'cc-pvtz-ri'          # RI拟合基
```

## 工具函数

### 坐标操作

```python
from pyscf.tools import cubegen, molden

# 立方体文件（用于可视化）
cubegen.orbital(mol, 'h2o_homo.cube',
                mo_coeff[:, mo_occ>0][:, -1])
cubegen.density(mol, 'h2o_density.cube',
                mf.make_rdm1())

# Molden文件（轨道可视化）
molden.from_mo(mol, 'h2o.molden', mf.mo_coeff)
```

### 布居分析

```python
from pyscf import lo

# Mulliken布居
pop = mf.mulliken_pop(mol, mf.make_rdm1())

# Lowdin布居
pop = lo.vec_lowdin(mf.mo_coeff, mf.get_ovlp())

# 自然布居分析（NPA）
# 需要外部工具或手动实现
```

## 并行计算

### OMP并行

```python
import os

# 设置线程数
os.environ['OMP_NUM_THREADS'] = '8'

# 自动并行
mf.kernel()
```

### MPI并行

```python
# 通过mpi4py运行
# mpirun -np 4 python script.py

from pyscf import lib

# 检查MPI
print(lib.num_threads())  # 线程数
print(lib.nproc)          # 进程数（如果使用MPI）
```

## 常见错误处理

### SCF不收敛

```python
# 检查收敛信息
print(mf.converged)
print(mf.scf_summary)

# 常见修复方法
mf.init_guess = 'huckel'    # 更好的初始猜测
mf.diis_start_cycle = 5    # 延迟DIIS
mf.level_shift = 0.5        # 水平位移
mf.damp_factor = 0.2       # 阻尼

# Newton-Raphson方法
mf_newton = mf.newton()
mf_newton.kernel()
```

### 内存不足

```python
# 密度拟合
mf = mf.density_fit()

# 外存计算
mf.direct_scf = False
mf.chkfile = 'calc.chk'

# 减小网格
mf.grids.atom_grid = (75, 302)
```

### 激发态问题

```python
# 先用TDA
td_tda = tdscf.TDA(mf)

# 检查基态质量
print(f'能隙: {gap*27.2114:.2f} eV')

# 增加网格精度
mf.grids.atom_grid = (99, 590)
```

## 检查点文件

### 保存所有数据

```python
# 自动保存
mf.chkfile = 'calculation.chk'
mf.kernel()

# 手动保存
mf.dump_chk(mf.chkfile, mf.mo_coeff, mf.mo_energy,
            mf.mo_occ, mf.e_tot)
```

### 从检查点读取

```python
# 读取分子
mol_chk = scf.chkfile.load('calculation.chk', 'mol')

# 读取MO
mo_coeff = scf.chkfile.load('calculation.chk', 'mo_coeff')
mo_energy = scf.chkfile.load('calculation.chk', 'mo_energy')
mo_occ = scf.chkfile.load('calculation.chk', 'mo_occ')

# 读取密度
dm = scf.chkfile.load('calculation.chk', 'dm')
```

## 性能优化建议

### 内存vs速度权衡

```python
# 大内存：存储所有积分
mf.direct_scf = False

# 小内存：即时计算积分
mf.direct_scf = True
```

### 基组选择

```python
# 快速预计算
basis = 'sto-3g'

# 标准计算
basis = 'def2-svp'

# 高精度
basis = 'def2-tzvp'

# 关键精度
basis = 'def2-qzvp'
```

### 并行效率

```python
# SCF并行：效率好
mf.kernel()

# 密度拟合：并行效率好
mf_df = mf.density_fit()

# MP2/CCSD：内存限制，并行效率中等
```

## PySCF扩展

### 自定义哈密顿量

```python
import numpy as np

# 1D Hubbard模型
h1 = np.zeros((N, N))
for i in range(N-1):
    h1[i, i+1] = h1[i+1, i] = -1.0

eri = np.zeros((N, N, N, N))
for i in range(N):
    eri[i, i, i, i] = U

# 构建SCF
mol = gto.M()
mol.nelectron = N
mf = scf.RHF(mol)
mf.get_hcore = lambda *args: h1
mf.get_ovlp = lambda *args: np.eye(N)
mf._eri = ao2mo.restore(8, eri, N)
mf.kernel()
```

### 与JAX集成

```python
import jax
import jax.numpy as jnp
from pyscf.jax import scf as jax_scf

# 转换SCF对象
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
mf = scf.RHF(mol)
mf_jax = jax_scf.RHF(mol)

# 自动微分梯度
grad_func = jax.grad(lambda coords: mf_jax.kernel()[0])
grad = grad_func(coords)
```

## 参考文献

1. PySCF文档: https://pyscf.org/
2. PySCF GitHub: https://github.com/pyscf/pyscf
3. PySCF示例: https://github.com/pyscf/pyscf/tree/master/examples
4. Sun, Q. et al. "Recent developments in the PySCF program package." WIREs Comput. Mol. Sci. 2020.
