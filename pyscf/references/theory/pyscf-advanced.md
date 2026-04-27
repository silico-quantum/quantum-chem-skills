# PySCF高级功能参考

## 自定义交换相关泛函

### 1. 泛函基本结构

#### 泛函接口
```python
from pyscf.dft import libxc, numint

# 泛函定义函数
def my_xc(rho, spin=0, deriv=1, **kwargs):
    """
    自定义交换相关泛函

    参数:
        rho: (N, 4) 或 (N, 2) 数组
              spin=0: [rho, grad_x, grad_y, grad_z]
              spin=1: [rho_a, rho_b, grad_ax, grad_ay, grad_az,
                       grad_bx, grad_by, grad_bz]
        spin: 自旋极化 (0: 限制性, 1: 非限制性)
        deriv: 导数阶数 (1: 能量和势, 2: 势的导数等)

    返回:
        deriv=1: (ex+vc, vrho, vgamma)
            ex+vc: (N,) 交换相关能密度
            vrho: (N, 2) 密度势 (dE/drho_a, dE/drho_b)
            vgamma: (N, 3) 梯度势 (dE/dgamma_aa, dE/dgamma_ab, dE/dgamma_bb)

            其中 gamma = grad_rho·grad_rho
    """
    # 提取密度
    if spin == 0:
        r = rho[:, 0]
        grho2 = rho[:, 1]**2 + rho[:, 2]**2 + rho[:, 3]**2
    else:
        ra = rho[:, 0]
        rb = rho[:, 1]
        grho2_a = rho[:, 2]**2 + rho[:, 3]**2 + rho[:, 4]**2
        grho2_b = rho[:, 5]**2 + rho[:, 6]**2 + rho[:, 7]**2

    # 你的泛函实现
    exc = ...  # 交换相关能密度
    vrho = ...  # 密度势
    vgamma = ...  # 梯度势

    return exc, vrho, vgamma
```

### 2. LDA泛函实现

#### Slater交换 (LDA)
```python
def lda_exchange(rho, spin=0, deriv=1):
    """
    Slater交换泛函
    Ex = -3/4 * (3/π)^(1/3) * ∫ ρ^(4/3) dr
    """
    if spin == 0:
        r = rho[:, 0]
        # 交换能密度
        exc = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        if deriv == 1:
            # Vx = dEx/drho = 4/3 * Ex/rho
            vrho = exc * (4/3) / (r + 1e-12)
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    else:
        ra = rho[:, 0]
        rb = rho[:, 1]
        # 自旋极化LDA交换
        exc = -0.75 * (3/np.pi)**(1/3) * (
            ra**(4/3) + rb**(4/3)
        )

        if deriv == 1:
            vrho_a = exc * (4/3) / (ra + 1e-12)
            vrho_b = exc * (4/3) / (rb + 1e-12)
            vrho = np.stack([vrho_a, vrho_b], axis=1)
            vgamma = np.zeros((len(rho), 3))
            return exc, vrho, vgamma

    return exc
```

#### VWN相关 (LDA)
```python
def vwn_correlation(rho, spin=0, deriv=1):
    """
    Vosko-Wilk-Nusair (VWN)相关泛函
    """
    # VWN参数 (参数来自原始论文)
    A = 0.0310907
    x0 = -0.10498
    b = 3.72744
    c = 12.9352

    def vwn_func(rs):
        """
        rs = (3/(4πρ))^(1/3)  Wigner-Seitz半径
        """
        x = np.sqrt(rs)
        X = x + b*x0 + c*x0**2
        Q = np.sqrt(4*c - b**2)

        fx = X**2 + Q**2
        fx0 = x0**2 + b*x0 + c*x0**2

        ec = A * (np.log(x**2/fx) +
                  2*b/Q * np.arctan((2*x + b)/Q) -
                  (b*x0)/Q * np.arctan((2*x0 + b)/Q) -
                  (x0*x)/fx)

        # 导数
        dec_drs = -A/3 * (
            1/x - (2*x + b)/X +
            2*b/Q * (1/(2*x + b + Q*x) + 1/(2*x + b - Q*x))
        )

        return ec, dec_drs

    if spin == 0:
        r = rho[:, 0]
        rs = (3/(4*np.pi*r))**(1/3)
        ec, dec_drs = vwn_func(rs)

        exc = ec

        if deriv == 1:
            # Vc = dEc/drho = dEc/drs * drs/drho
            # drs/drho = -rs/(3*rho)
            drs_drho = -rs/(3*r)
            vrho = dec_drs * drs_drho
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    else:
        # 自旋极化VWN (更复杂)
        ra = rho[:, 0]
        rb = rho[:, 1]
        r = ra + rb
        zeta = (ra - rb) / r  # 自旋极化参数

        rs = (3/(4*np.pi*r))**(1/3)
        ec0, dec_drs = vwn_func(rs)  # 无自旋极化

        # 自旋极化校正 (简化版)
        alpha = 0.001  # 参数
        exc = ec0 * (1 + alpha * f(zeta))

        if deriv == 1:
            # 完整的自旋极化导数更复杂
            vrho_a = ...  # dE/drho_a
            vrho_b = ...  # dE/drho_b
            vrho = np.stack([vrho_a, vrho_b], axis=1)
            vgamma = np.zeros((len(rho), 3))
            return exc, vrho, vgamma

    return exc
```

### 3. GGA泛函实现

#### Becke88交换 (B88)
```python
def becke88_exchange(rho, spin=0, deriv=1):
    """
    Becke 1988 GGA交换泛函
    Ex = ∫ εx^LDA * Fx(s) dr
    其中 s = |∇ρ| / (2kF ρ) 是约化密度梯度
    """
    beta = 0.0042  # Becke参数

    if spin == 0:
        r = rho[:, 0]
        grx = rho[:, 1]
        gry = rho[:, 2]
        grz = rho[:, 3]

        grho2 = grx**2 + gry**2 + grz**2

        # 约化密度梯度
        kF = (3*np.pi**2*r)**(1/3)
        s = np.sqrt(grho2) / (2*kF*r + 1e-12)

        # LDA交换
        ex_lda = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        # Becke增强因子
        s2 = s**2
        asinh_s = np.arcsinh(s)
        Fx = 1 + beta * s2 / (1 + 6*beta*s*asinh_s)

        # 交换能密度
        exc = ex_lda * Fx

        if deriv == 1:
            # 势计算需要dFx/ds
            # Vx = dEx/drho - ∇·(dEx/d∇ρ)
            # 这是一个复杂的泛函导数计算

            # 简化版：仅返回LDA势
            vrho = ex_lda * (4/3) / (r + 1e-12) * Fx
            vgamma = np.zeros_like(rho[:, 0])

            # 完整的GGA势需要tau项
            return exc, vrho, vgamma

    return exc
```

#### PBE交换
```python
def pbe_exchange(rho, spin=0, deriv=1):
    """
    Perdew-Burke-Ernzerhof交换泛函
    """
    kappa = 0.804
    mu = 0.21951

    if spin == 0:
        r = rho[:, 0]
        grho2 = rho[:, 1]**2 + rho[:, 2]**2 + rho[:, 3]**2

        kF = (3*np.pi**2*r)**(1/3)
        s = np.sqrt(grho2) / (2*kF*r + 1e-12)

        # PBE增强因子
        Fx = 1 + kappa - kappa / (1 + mu*s**2/kappa)

        # LDA交换
        ex_lda = -0.75 * (3/np.pi)**(1/3) * r**(4/3)

        exc = ex_lda * Fx

        if deriv == 1:
            # PBE势 (简化版)
            vrho = ex_lda * (4/3) / (r + 1e-12) * Fx
            vgamma = np.zeros_like(rho[:, 0])
            return exc, vrho, vgamma

    return exc
```

### 4. 杂化泛函

#### B3LYP手动实现
```python
def b3lyp_functional(mol):
    """
    手动实现B3LYP泛函
    Ex = 0.2*Ex^HF + 0.08*Ex^LDA + 0.72*Ex^B88
    Ec = 0.19*Ec^VWN + 0.81*Ec^LYP
    """
    from pyscf import dft

    mf = dft.RKS(mol)

    # 设置杂化参数
    mf.xc = {
        'hybrid': {
            'hf_fraction': 0.2,  # 20% HF交换
        },
        'RXC': {
            'type': 'GGA',
            'exch': {
                'LDA': 0.08,
                'B88': 0.72,
            },
            'corr': {
                'VWN': 0.19,
                'LYP': 0.81,
            }
        }
    }

    # 或者更简单的方式
    mf.xc = 'b3lyp'

    return mf
```

#### 自定义杂化泛函
```python
def custom_hybrid(mol, hf_fraction=0.25, xc_dft='pbe'):
    """
    自定义杂化泛函
    Ex = α*Ex^HF + (1-α)*Ex^DFT
    """
    mf = dft.RKS(mol)

    # 通过字符串指定
    mf.xc = f'{hf_fraction}*HF + {1-hf_fraction}*{xc_dft}'

    # 或者通过字典
    mf.xc = f'{hf_fraction}*HF,{xc_dft}'

    return mf
```

### 5. 应用自定义泛函

```python
from pyscf import gto, dft, lib

# 定义分子
mol = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# 方法1: 直接赋值泛函函数
mf = dft.RKS(mol)
mf._numint._xc = my_xc  # 替换泛函计算函数
e = mf.kernel()

# 方法2: 通过libxc注册
# 注册新泛函
libxc.define_xc_ = libxc.define_xc_.copy()
libxc.define_xc_('MY_XC', my_xc)

# 使用
mf = dft.RKS(mol)
mf.xc = 'MY_XC'
e = mf.kernel()

# 方法3: 使用NumInt对象
from pyscf.dft import numint

ni = numint.NumInt()
ni._xc = my_xc

mf = dft.RKS(mol)
mf._numint = ni
e = mf.kernel()
```

## 积分操作

### 1. 双电子积分

#### 原始双电子积分 (8重对称)
```python
from pyscf import gto

mol = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)

# 完全双电子积分 (8重对称)
# (ij|kl) ~ (ji|kl) ~ (ij|lk) ~ ...
eri = mol.intor('int2e')
print(f"积分形状: {eri.shape}")  # (nao, nao, nao, nao)

# 特定积分
i, j, k, l = 0, 0, 0, 0
value = eri[i, j, k, l]
print(f"({i}{j}|{k}{l}) = {value:.6f}")
```

#### 4重对称积分
```python
# 4重对称: (ij|kl) = (kl|ij)
eri_4 = mol.intor('int2e_sph', aosym=4)
print(f"4重对称形状: {eri_4.shape}")  # (nao*(nao+1)/2, nao*(nao+1)/2)
```

#### 密度拟合积分
```python
from pyscf import df

# 创建DF对象
dfobj = df.DF(mol)
dfobj.auxbasis = 'cc-pvdz-ri'

# 计算3中心积分 (ij|A)
ijA = dfobj.get_2c2e()  # (naux, nao, nao)

# 计算2中心积分 (A|B)
AB = dfobj.get_2c2e_aux()  # (naux, naux)

# 拟合双电子积分
eri_df = dfobj.get_j()  # 近似双电子积分
```

### 2. AO到MO积分变换

#### 完全变换 (内存密集)
```python
from pyscf import ao2mo

# AO基双电子积分
eri_ao = mol.intor('int2e')

# MO系数
C = mf.mo_coeff  # (nao, nmo)

# 完全变换到MO基
eri_mo = ao2mo.incore.full(eri_ao, C)
print(f"MO积分形状: {eri_mo.shape}")  # (nmo, nmo, nmo, nmo)

# 特定MO积分
i, j, k, l = 0, 1, 2, 3  # 轨道索引
value = eri_mo[i, j, k, l]
print(f"({i}{j}|{k}{l})_MO = {value:.6f}")
```

#### 部分变换 (占据-占据)
```python
# 占据轨道
nocc = mol.nelectron // 2
C_occ = C[:, :nocc]  # (nao, nocc)

# 仅计算占据-占据积分
eri_oooo = ao2mo.incore.general(eri_ao,
                                [C_occ, C_occ, C_occ, C_occ])
print(f"(oo|oo)形状: {eri_oooo.shape}")  # (nocc*nocc, nocc*nocc)
```

#### 外存变换 (大分子)
```python
# 保存到磁盘
ao2mo.outcore.full(mol, C, 'eri_mo.h5')

# 读取特定块
from pyscf import lib
with lib.H5File('eri_mo.h5', 'r') as f:
    eri_block = f['eri_mo'][:10, :10, :10, :10]

# 或者使用内存映射
eri_mo = ao2mo.load('eri_mo.h5')
```

### 3. 积分导数

#### 基函数导数积分
```python
# 基函数梯度: d/dx χ_a(r)
# dS_ij/dx_A = ∫ [dχ_i/dx_A χ_j + χ_i dχ_j/dx_A] dr

# 重叠矩阵导数
dS_dx = mol.intor('int1e_ovlp_ip1', comp=3)  # (3, nao, nao)
# dS_dx[0] = d/dx, dS_dx[1] = d/dy, dS_dx[2] = d/dz

# 动能矩阵导数
dT_dx = mol.intor('int1e_kin_ip1', comp=3)

# 核吸引矩阵导数
dV_dx = mol.intor('int1e_nuc_ip1', comp=3)
```

#### 双电子积分导数
```python
# 双电子积分导数 (3中心)
# d(ij|kl)/dA = d/dx_A (ij|kl)
eri_ip1 = mol.intor('int2e_ip1', comp=3)
print(f"积分导数形状: {eri_ip1.shape}")  # (3, nao, nao, nao, nao)
```

### 4. 相对论积分

#### DKH (Douglas-Kroll-Hess)
```python
# Douglas-Kroll-Hess哈密顿量
from pyscf.dkh import dkh

# DKH2阶
dkh2 = dkh.make_dkh2_hamiltonian(mol)

# DKH修正的动能和势能
T_dkh, V_dkh = dkh2

# 在SCF中使用
mf = scf.RHF(mol)
mf.get_hcore = lambda *args: T_dkh + V_dkh
mf.kernel()
```

#### ZORA (Zeroth-Order Regular Approximation)
```python
from pyscf.sfx2c1e import sfx2c1e

# Spin-free X2C1e
mf = scf.RHF(mol)
mf = sfx2c1e.get_x2c(mf)  # 转换为相对论SCF
mf.kernel()
```

## SCF高级控制

### 1. DIIS加速

#### 自定义DIIS
```python
from pyscf.scf import diis

# DIIS类
my_diis = diis.DIIS()

# SCF循环
for cycle in range(max_cycle):
    # 构建Fock矩阵
    fock = build_fock(dm)

    # DIIS外推
    if cycle >= diis_start:
        fock = my_diis.update(fock, x=dm)

    # 对角化
    mo_energy, mo_coeff = eig(fock, S)

    # 新密度
    dm_new = build_density(mo_coeff)

    # 收敛检查
    if converged(dm, dm_new):
        break

    dm = dm_new
```

#### DIIS误差向量
```python
# 自定义误差向量
def error_vector(fock, dm, S):
    """
    DIIS误差: [F,D]S = FDS - SDF
    """
    return fock @ dm @ S - S @ dm @ fock

# 在DIIS中使用
my_diis = diis.DIIS()
error = error_vector(fock, dm, S)
fock = my_diis.update(fock, x=dm, errvec=error)
```

### 2. 直接SCF

#### 控制直接计算
```python
mf = scf.RHF(mol)

# 完全直接 (每次重新计算积分)
mf.direct_scf = True

# 半直接 (缓存一些积分)
mf.direct_scf = False
mf._eri = mol.intor('int2e')  # 预计算所有积分

# 混合模式
mf.direct_scf = True
mf.direct_scf_tol = 1e-12  # 积分阈值
```

#### 内存管理
```python
# 设置最大内存
lib.num_threads(4)  # 4线程
lib.chkfile.save_temp(mol.chkfile, mf.get_fock())  # 保存中间结果

# 外存SCF
mf.chkfile = 'calc.chk'
mf.kernel()
# 中间结果会保存到检查点
```

### 3. SCF稳定化

#### 稳定性分析
```python
# 检查SCF解的稳定性
mf = scf.RHF(mol)
mf.kernel()

# 稳定性分析
from pyscf.scf import stability

# 内部稳定性 (UHF不稳定?)
stable = stability.stability_internal(mf)
print(f"内部稳定: {stable}")

# 外部稳定性 (存在更低能解?)
stable2, e_new, mf_new = stability.stability_external(mf)
print(f"外部稳定: {stable2}")

# 如果不稳定，使用新的SCF
if not stable2:
    mf_new.kernel()
```

#### 阻尼和水平位移
```python
# 阻尼 (振荡问题)
mf = scf.RHF(mol)
mf.damp = 0.2  # 20%阻尼
mf.kernel()

# 水平位移 (收敛慢)
mf = scf.RHF(mol)
mf.level_shift = 0.5  # 虚轨道位移0.5 Ha
mf.kernel()

# 组合使用
mf.damp = 0.2
mf.level_shift = 0.5
mf.diis_start_cycle = 5  # 延迟DIIS
```

### 4. Newton-Raphson SCF

#### 二阶SCF
```python
# Newton-Raphson方法 (快速收敛，但不稳定)
mf = scf.RHF(mol)
mf_nr = mf.newton()  # 创建NR-SCF
mf_nr.kernel()

# 带约束的NR-SCF
mf_nr.constrain_occ = True  # 保持占据
mf_nr.kernel()
```

#### 奇数电子系统
```python
# ROHF稳定性
mol = gto.M(atom='C 0 0 0', charge=0, spin=2, basis='sto-3g')
mf = scf.ROHF(mol)
mf.kernel()

# 转换为UHF并稳定化
mf_uhf = scf.UHF(mol)
mf_uhf.kernel()

# 稳定性检查
mf_uhf_stable = stability.stability_internal(mf_uhf)
```

## TDDFT高级功能

### 1. 自然跃迁轨道 (NTO)

#### 深度NTO分析
```python
from pyscf import tdscf

td = tdscf.TDDFT(mf)
td.nstates = 5
td.kernel()

# 所有态的NTO
for n in range(td.nstates):
    weights, nto = td.get_nto(state=n)

    print(f"\n态 {n+1}:")
    print(f"  激发能: {td.e[n]*27.2114:.2f} eV")
    print(f"  主导权重: {weights.max():.4f}")

    # 空穴和电子轨道
    hole = nto[0]  # (nao, nocc)
    electron = nto[1]  # (nao, nvirt)

    # 主导NTO对
    dominant_pair = np.argmax(weights)
    print(f"  主导对权重: {weights[dominant_pair]:.4f}")

    # 分析主导NTO的轨道成分
    hole_orb = hole[:, dominant_pair]  # (nao,)
    electron_orb = electron[:, dominant_pair]

    # 轨道成分分析
    from pyscf.lo import iao
    hole_iao = iao.mulliken_pop(mol, hole_orb)
    electron_iao = iao.mulliken_pop(mol, electron_orb)
```

#### NTO可视化
```python
# 保存NTO到Molden文件
from pyscf.tools import molden

for n in range(td.nstates):
    weights, nto = td.get_nto(state=n)

    # 主导空穴轨道
    hole = nto[0][:, np.argmax(weights)]
    # 主导电子轨道
    electron = nto[1][:, np.argmax(weights)]

    # 保存
    molden.from_mo(mol, f'hole_state{n+1}.molden', hole)
    molden.from_mo(mol, f'electron_state{n+1}.molden', electron)
```

### 2. TDDFT密度分析

#### 激发态密度差
```python
# 基态密度
rho0 = mf.get_rho()

# 激发态密度 (态1)
dm1 = td.generate_density(state=0)  # 第1激发态
rho1 = td.get_rho(dm1)

# 密度差
delta_rho = rho1 - rho0

# 保存密度差
from pyscf.tools import cubegen
cubegen.density(mol, 'delta_rho.cube', dm1 - mf.make_rdm1())

# 分析
print(f"密度差积分 (应接近0): {delta_rho.sum():.6f}")
print(f"密度差极值: {delta_rho.min():.6f}, {delta_rho.max():.6f}")
```

#### 电荷转移分析
```python
def charge_transfer_analysis(mol, dm_excited):
    """
    分析激发态的电荷转移
    """
    from pyscf import lo

    # Mulliken布居
    pop_excited = mf.mulliken_pop(mol, dm_excited)
    pop_ground = mf.mulliken_pop(mol, mf.make_rdm1())

    # 电荷变化
    delta_charge = pop_excited[0] - pop_ground[0]

    print("电荷变化 (e):")
    for i in range(mol.natm):
        atom = mol.atom_pure_symbol(i)
        print(f"  {atom}{i}: {delta_charge[i]:+.4f}")

    # 电荷转移距离
    # d_CT = |Σ q_i * r_i| / Σ |q_i|
    coords = np.array([mol.atom_coord(i) for i in range(mol.natm)])
    ct_distance = np.abs(np.sum(delta_charge[:, None] * coords, axis=0))
    ct_distance /= np.sum(np.abs(delta_charge))

    print(f"\n电荷转移距离: {ct_distance:.3f} Å")

    return delta_charge

# 分析第一激发态
dm1 = td.generate_density(state=0)
charge_transfer_analysis(mol, dm1)
```

### 3. 自旋轨道耦合 (SOC)

#### 单组态相互作用
```python
from pyscf import lib

# 自旋轨道耦合矩阵
# SOC = <Ψ_S,M_S|H_SO|Ψ'_S',M'_S'>

# 对于TDDFT态，需要计算SOC矩阵元
# 这是一个高级功能，需要扩展PySCF

# 简化方案: 使用PySCF的spin-orbit模块
from pyscf.spinorbit import soc

# SOC需要4分量 relativistic SCF
from pyscf.x2c import X2C
mol_r = X2C(mol)
mf_r = scf.RHF(mol_r)
mf_r.kernel()

# SOC计算
soc_mat = soc.spin_orbit_coulomb(mf_r, td)
print(f"SOC矩阵形状: {soc_mat.shape}")

# 耦合强度
print(f"SOC强度: {np.linalg.norm(soc_mat):.4f} cm^-1")
```

## MP2和后-HF方法

### 1. MP2

#### RMP2计算
```python
from pyscf import mp

# MP2对象
mp2 = mp.MP2(mf)
e_mp2, t2 = mp2.kernel()

print(f"MP2相关能: {e_mp2:.6f} Ha")
print(f"总能量 (HF+MP2): {mf.e_tot + e_mp2:.6f} Ha")

# MP2能量分解
# E(MP2) = E(HF) + E_2
# E_2 = Σ_{ijab} |<ij||ab>|^2 / (ε_i+ε_j-ε_a-ε_b)

# 轨道能量
eps_i = mf.mo_energy[:nocc]  # 占据
eps_a = mf.mo_energy[nocc:]  # 虚轨道

# 顶点: <ij||ab> = <ij|ab> - <ij|ba>
# 存储在t2中

# 诊断参数
# D1 = Σ_{ijab} |t2_{ijab}|^2
# T1 = Σ_{ia} |t1_{ia}|  (for MP3+)

d1 = mp2.get_d1()
t1 = mp2.get_t1()

print(f"D1诊断: {d1:.3f} (理想 < 0.02)")
print(f"T1诊断: {t1:.3f} (理想 < 0.02)")

# 多参考警告
if d1 > 0.02 or t1 > 0.02:
    print("警告: 多参考效应显著，考虑CASSCF")
```

#### 局域MP2 (LMP2)
```python
# 密度拟合MP2
mp2_df = mp.DFMP2(mf)
mp2_df.auxbasis = 'cc-pvtz-ri'
e_mp2_df = mp2_df.kernel()[0]

print(f"DF-MP2相关能: {e_mp2_df:.6f} Ha")
print(f"与标准MP2差: {abs(e_mp2 - e_mp2_df):.6f} Ha")
```

#### RI-MP2 (密度拟合)
```python
# 使用Resolution-of-Identity
mp2_ri = mp.MP2(mf)
mp2_ri.with_df = df.DF(mol)
mp2_ri.with_df.auxbasis = 'cc-pvdz-ri'
e_mp2_ri = mp2_ri.kernel()[0]

print(f"RI-MP2相关能: {e_mp2_ri:.6f} Ha")
```

### 2. 耦合簇

#### CCSD
```python
from pyscf import cc

# CCSD对象
ccsd = cc.CCSD(mf)
e_ccsd, t1, t2 = ccsd.kernel()

print(f"CCSD相关能: {e_ccsd:.6f} Ha")
print(f"总能量: {mf.e_tot + e_ccsd:.6f} Ha")

# T1诊断
t1_norm = np.linalg.norm(t1)
print(f"T1诊断: {t1_norm:.3f}")

# CCSD(T) - 黄金标准
e_ccsdt = ccsd.ccsd_t()
print(f"(T)校正: {e_ccsdt:.6f} Ha")
print(f"CCSD(T)总能量: {mf.e_tot + e_ccsd + e_ccsdt:.6f} Ha")
```

#### 分解
```python
# E(CCSD) = Σ_{ia} t_{ia} <ia||ia> + Σ_{ijab} t_{ijab} <ij||ab>

# 部分能量
e_1 = ccsd.energy1()  # 单激发
e_2 = ccsd.energy2()  # 双激发

print(f"E(单): {e_1:.6f} Ha")
print(f"E(双): {e_2:.6f} Ha")
```

#### Lambda方程 (激发态)
```python
# CCSD Lambda方程 (用于激发态)
ccsd_lambda = cc.ccsd_lambda.CCSD_Lambda(ccsd)

# 基态Lambda
l1, l2 = ccsd_lambda.kernel()

# EOM-CCSD (方程运动耦合簇)
from pyscf import cc

eom = cc.eom_rccsd.EOMIP(ccsd)  # 电离
# eom = cc.eom_rccsd.EOMEA(ccsd)  # 电子亲和
# eom = cc.eom_rccsd.EOMEES(ccsd)  # 激发态

ip_e, ip_vec = eom.ipccsd(nroots=3)  # 3个电离态
print(f"垂直电离能 (eV):")
for i, e in enumerate(ip_e):
    print(f"  态{i+1}: {e*27.2114:.2f} eV")
```

### 3. CI方法

#### CIS (构型相互作用单激发)
```python
from pyscf import ci

# CIS对象
cis = ci.CIS(mf)
cis.nstates = 5
e_cis, civec = cis.kernel()

print(f"CIS激发态:")
for i in range(cis.nstates):
    print(f"  态{i+1}: {e_cis[i]*27.2114:.2f} eV")
```

#### CISD
```python
# CISD (CIS + 双激发)
cisd = ci.CISD(mf)
e_cisd, ci_vec = cisd.kernel()

print(f"CISD相关能: {e_cisd - mf.e_tot:.6f} Ha")

# Davidson算法
cisd.max_space = 100  # Davidson子空间大小
e_cisd, ci_vec = cisd.kernel()
```

### 4. FCI (全CI)

#### 小系统FCI
```python
from pyscf import fci

# 小分子
mol_small = gto.M(
    atom='H 0 0 0; H 0 0 0.74',
    basis='sto-3g'
)
mf_small = scf.RHF(mol_small)
mf_small.kernel()

# FCI计算
cisolver = fci.FCI(mf_small)
e_fci, ci_vec = cisolver.kernel()

print(f"FCI能量: {e_fci:.6f} Ha")
print(f"FCI vs HF: {(e_fci - mf_small.e_tot)*27.2114:.2f} eV")

# 精确波函数分析
print(f"CI系数数: {len(ci_vec)}")
print(f"主导构型权重: {abs(ci_vec).max():.4f}")
```

#### 选定CI
```python
# 有限空间选定的CI
fcisolver = fci.FCI(mf_small)
fcisolver.nroots = 3  # 多个态

# 选定空间
norb = mol_small.nao
nelec = mol_small.nelectron

# 手动选择轨道空间
# 这里使用完整FCI
e_fci, civec = fcisolver.kernel()
```

## 密度拟合 (DF)

### 1. 基础DF

#### DF-SCF
```python
from pyscf import df

# 自动DF
mf = scf.RHF(mol)
mf_df = mf.density_fit()

# 或手动指定拟合基组
mf_df = df.density_fit(mf, auxbasis='def2-universal-jfit')

e_df = mf_df.kernel()
print(f"DF-RHF能量: {e_df:.6f} Ha")

# 与标准RHF比较
print(f"DF vs RHF能量差: {abs(e_df - mf.kernel()):.8f} Ha")
```

#### DF对象
```python
# DF对象
dfobj = df.DF(mol)
dfobj.auxbasis = 'cc-pvdz-ri'

# 3中心积分 (ij|P)
ijP = dfobj.get_2c2e()
print(f"(ij|P)形状: {ijP.shape}")  # (naux, nao, nao)

# 2中心积分 (P|Q)
PQ = dfobj.get_2c2e_aux()
print(f"(P|Q)形状: {PQ.shape}")  # (naux, naux)

# 求解 (P|Q) * J^P = (ij|P)
# J = Σ_PQ (ij|P) * (P|Q)^(-1)
```

### 2. MP2-DF

#### DF-MP2
```python
from pyscf import mp

# MP2 + DF
mp2 = mp.MP2(mf)
mp2_df = mp.density_fit(mp2, auxbasis='cc-pvdz-ri')

e_mp2_df = mp2_df.kernel()[0]
print(f"DF-MP2相关能: {e_mp2_df:.6f} Ha")

# 检查与标准MP2的差别
e_mp2_std = mp.MP2(mf).kernel()[0]
print(f"DF vs 标准MP2: {abs(e_mp2_df - e_mp2_std):.6f} Ha")
```

### 3. CCSD-DF

#### DF-CCSD
```python
# CCSD + DF
ccsd = cc.CCSD(mf)
ccsd_df = ccsd.density_fit()

e_ccsd_df = ccsd_df.kernel()[0]
print(f"DF-CCSD相关能: {e_ccsd_df:.6f} Ha")
```

## 周期体系

### 1. 周期性分子定义

#### 1D周期
```python
from pyscf.pbc import gto as pbcgto

# 1D晶胞 (线性链)
cell = pbcgto.M(
    atom='H 0 0 0',
    basis='sto-3g',
    a=[[2.0, 0, 0], [0, 20, 0], [0, 0, 20]],  # 晶格向量
    unit='Ang',
)

# k点采样
kpts = cell.make_kpts([2, 1, 1])  # 2个k点
```

#### 2D周期
```python
# 2D晶胞 (石墨烯层)
cell = pbcgto.M(
    atom='''
    C  0.0000  0.0000  0.0000
    C  0.0000  1.4200  0.0000
    C  1.2308  0.7100  0.0000
    C  1.2308  2.1300  0.0000
    ''',
    basis='gth-szv',
    a=[[2.4616, 0, 0], [0, 4.2600, 0], [0, 0, 20]],  # 平面周期
    pseudo='gth-pade',
    ke_cutoff=100,  # 动能截断 (eV)
)

# k点
kpts = cell.make_kpts([4, 4, 1])
```

#### 3D周期
```python
# 3D晶胞 (硅)
cell = pbcgto.M(
    atom='Si 0 0 0; Si 0.25 0.25 0.25',
    basis='gth-szv',
    a=[[5.43, 0, 0], [0, 5.43, 0], [0, 0, 5.43]],
    pseudo='gth-pade',
    ke_cutoff=100,
)

# k点网格
kpts = cell.make_kpts([4, 4, 4])
```

### 2. 周期SCF

#### Gamma-only
```python
from pyscf.pbc import scf as pbcscf

# 仅Gamma点
mf = pbcscf.RHF(cell)
e_gam = mf.kernel()
print(f"Gamma-only能量: {e_gam:.6f} Ha")
```

#### 多k点SCF
```python
# 多k点
mf = pbcscf.RHF(cell, kpts)
e_kpts = mf.kernel()
print(f"多k点能量: {e_kpts:.6f} Ha")

# k点权重
print(f"k点数: {len(kpts)}")
print(f"总权重: {mf.kpts.sum()}")
```

### 3. 周期DFT

#### PBC-DFT
```python
from pyscf.pbc import dft as pbcdft

# 周期DFT
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'pbe'
mf.grids.level = 3  # 网格精度
e_pbe = mf.kernel()
print(f"PBE能量: {e_pbe:.6f} Ha")
```

#### 杂化泛函
```python
# 周期杂化泛函 (计算昂贵)
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'hse06'  # 屏蔽杂化
e_hse = mf.kernel()
print(f"HSE06能量: {e_hse:.6f} Ha")

# 或使用精确交换的密度拟合
mf = pbcdft.RKS(cell, kpts)
mf.xc = 'pbe0'
mf = mf.density_fit()  # 加速HF交换
e_pbe0 = mf.kernel()
```

### 4. 能带结构

#### 能带计算
```python
# 计算能带路径
from pyscf.pbc import tools

# 高对称点路径
Gamma = [0, 0, 0]
X = [0.5, 0, 0]
M = [0.5, 0.5, 0]

# 生成k路径
kpath = tools.get_bandpath([Gamma, X, M], cell, 20)

# 在每个k点计算能带
band_energies = []
for k in kpath:
    mf_k = pbcdft.RKS(cell, k.reshape(1, 3))
    mf_k.xc = 'pbe'
    mf_k.kernel()
    band_energies.append(mf_k.mo_energy[0])

band_energies = np.array(band_energies)  # (nk, nmo)

# 绘制能带
import matplotlib.pyplot as plt
plt.plot(band_energies)
plt.xlabel('k-path')
plt.ylabel('Energy (Ha)')
plt.title('Band Structure')
plt.show()
```

## 溶剂效应

### 1. PCM (极化连续介质)

#### PCM-SCF
```python
from pyscf import solvent

# PCM模型
pcm = solvent.PCM(mol)
pcm.eps = 78.4  # 水的介电常数
pcm.method = 'IEFPCM'  # 积分方程形式PCM

# SCF with PCM
mf = scf.RHF(mol)
mf = pcm.run(mf)
e_pcm = mf.e_tot
print(f"PCM能量: {e_pcm:.6f} Ha")

# 溶剂化能
e_gas = scf.RHF(mol).kernel()
G_sol = e_pcm - e_gas
print(f"溶剂化能: {G_sol*27.2114:.2f} eV")
```

#### 不同溶剂
```python
# 不同溶剂
solvents = {
    'water': 78.4,
    'ethanol': 24.3,
    'methanol': 32.6,
    'acetone': 20.7,
    'dichloromethane': 8.93,
}

for name, eps in solvents.items():
    pcm = solvent.PCM(mol)
    pcm.eps = eps
    mf_pcm = scf.RHF(mol)
    mf_pcm = pcm.run(mf_pcm)
    print(f"{name:15s}: {mf_pcm.e_tot:.6f} Ha")
```

### 2. SMD模型

#### SMD-SCF
```python
# SMD (Solvation Model based on Density)
pcm = solvent.PCM(mol)
pcm.method = 'SMD'
pcm.solvent = 'water'  # 溶剂名称

mf = scf.RHF(mol)
mf = pcm.run(mf)
print(f"SMD能量: {mf.e_tot:.6f} Ha")
```

### 3. 显式溶剂

#### 微团簇模型
```python
# 水+6个水分子的微团簇
cluster = gto.M(
    atom='''
    O  0.000000  0.000000  0.000000
    H  0.000000  0.758602 -0.504284
    H  0.000000 -0.758602 -0.504284
    O  2.800000  0.000000  1.500000
    H  3.558602  0.000000  1.995716
    H  2.041398  0.000000  1.995716
    ... (更多水分子)
    ''',
    basis='cc-pvdz'
)

mf_cluster = scf.RHF(cluster)
e_cluster = mf_cluster.kernel()
print(f"团簇能量: {e_cluster:.6f} Ha")
```

## 并行计算

### 1. OpenMP并行

#### 设置线程数
```python
import os

# 设置OpenMP线程数
os.environ['OMP_NUM_THREADS'] = '8'

# 或使用PySCF设置
lib.num_threads(8)

# SCF自动并行
mf = scf.RHF(mol)
e = mf.kernel()
print(f"使用的线程数: {lib.num_threads()}")
```

#### 性能分析
```python
# 启用计时
lib.logger.TIMER_LEVEL = 3

mf = scf.RHF(mol)
e = mf.kernel()

# 查看计时信息
# 计时信息会输出到日志
```

### 2. MPI并行

#### MPI-SCF
```python
# 需要使用mpi4py启动
# mpirun -np 4 python script.py

from pyscf import lib

# 检查MPI
if hasattr(lib, 'mpi'):
    print(f"MPI进程数: {lib.mpi.rank} of {lib.mpi.size}")
else:
    print("MPI未启用")

# SCF自动并行
mf = scf.RHF(mol)
e = mf.kernel()
```

## 参考文献和资源

### PySCF文档
1. PySCF官方文档: https://pyscf.org/
2. PySCF示例: https://github.com/pyscf/pyscf/tree/master/examples
3. PySCF API文档: https://pyscf.org/api.html

### 论文
4. Sun et al., "Recent developments in the PySCF program package", WIREs Comput Mol Sci, 2020
5. Sun et al., "PySCF: the Python-based simulations of chemistry framework", WIREs Comput Mol Sci, 2018

### 教程和课程
6. PySCF workshop: https://pyscf.org/workshop/
7. Quantum chemistry with PySCF: https://pyscf.org/tutorial.html
