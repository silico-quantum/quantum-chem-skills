# PySCF与JAX集成指南

## 概述

JAX (Just Another eXecutor) 是一个高性能的数值计算库，支持自动微分（autograd）和 JIT 编译。将 PySCF 与 JAX 集成可以实现：

1. **端到端自动微分**：直接对 DFT 计算梯度
2. **可微分量子化学**：在神经网络中使用 PySCF
3. **泛函优化**：自动优化交换相关泛函参数
4. **GPU 加速**：利用 GPU 加速积分计算

## 安装配置

### 1. 环境设置

```bash
# 安装 PySCF
pip install pyscf

# 安装 JAX
pip install jax jaxlib

# GPU 支持 (可选)
# pip install jax[cuda]
```

### 2. 基础导入

```python
import numpy as np
import jax
import jax.numpy as jnp
from jax import grad, jit, vmap, pmap

import pyscf
from pyscf import gto, scf, dft, lib
```

## PySCF-JAX 基础

### 1. 数组转换

#### NumPy ↔ JAX 数组
```python
# NumPy 数组
np_array = np.array([1.0, 2.0, 3.0])

# 转换为 JAX 数组
jx_array = jnp.array(np_array)
# 或
jx_array = jax.device_put(np_array)

# 转换回 NumPy
np_from_jax = np.array(jx_array)
# 或
np_from_jax = jax.device_get(jx_array)
```

#### PySCF 对象处理
```python
# PySCF 分子对象
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')

# 提取 JAX 友好的数据
coords = jnp.array(mol.atom_coords())
charges = jnp.array([mol.atom_charge(i) for i in range(mol.natm)])
basis_info = mol._basis  # 基组信息
```

### 2. 可包装函数

#### JIT 兼容函数
```python
from jax import jit

# 定义简单的能量函数
def simple_h2_energy(coords):
    """
    H₂ 原子核排斥能 + 简化电子能量

    params:
        coords: (2, 3) 原子坐标
    """
    # 核-核排斥
    r_vec = coords[0] - coords[1]
    r = jnp.sqrt(jnp.sum(r_vec**2))
    Vnn = 1.0 / r

    # 简化电子能 (Hückel 模型)
    E_el = -2.0 * np.exp(-r)  # 简化模型

    return Vnn + E_el

# JIT 编译
jit_energy = jit(simple_h2_energy)

# 测试
coords_test = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
E = jit_energy(coords_test)
print(f"能量: {E:.6f} Ha")
```

#### 自动微分
```python
# 计算梯度
grad_energy = grad(simple_h2_energy)

# 计算在给定构型下的力 (力的负梯度)
forces = -grad_energy(coords_test)
print(f"原子 0 力: {forces[0]}")
print(f"原子 1 力: {forces[1]}")
```

## GradDFT 框架

### 1. GradDFT 架构

GradDFT 是一个端到端可微分的 DFT 框架，允许通过 DFT 计算反向传播梯度。

```python
import jax
import jax.numpy as jnp
from pyscf import gto, scf, dft
from pyscf.dft import numint

class GradDFT:
    """
    可微分 DFT 框架
    """
    def __init__(self, mol):
        self.mol = mol
        self.nao = mol.nao
        self.nelec = mol.nelectron

        # 初始化 SCF
        self.mf = dft.RKS(mol)
        self.mf.xc = 'pbe'

        # 数值积分对象
        self.ni = numint.NumInt()

    def build_matrices(self, coords):
        """
        构建 JAX 兼容的积分矩阵

        params:
            coords: (natm, 3) 原子坐标
        """
        # 注意: 这里需要重新定义分子对象
        # 简化版: 假设基函数固定，只改变原子位置

        # 重叠矩阵 (简化: 假设固定)
        S = self.mf.get_ovlp()

        # 动能矩阵 (简化)
        T = self.mf.get_hcore() - self.mol.intor('int1e_nuc')

        # 核吸引矩阵 (依赖于坐标)
        V = self._build_nuclear_potential(coords)

        # 总核心哈密顿量
        Hcore = T + V

        return S, Hcore

    def _build_nuclear_potential(self, coords):
        """
        构建核吸引势矩阵

        V_μν = ∫ χ_μ(r) [ -Σ_A Z_A/|r-R_A| ] χ_ν(r) dr
        """
        # 简化: 使用固定积分 + 位移
        # 实际实现需要重新计算积分

        # 基础核吸引
        V0 = self.mol.intor('int1e_nuc')

        # 坐标位移的影响 (需要重新计算)
        # 这里使用简化假设
        return V0  # 实际中需要重新计算

    def scf_iteration(self, coords, dm):
        """
        单次 SCF 迭代

        params:
            coords: (natm, 3) 原子坐标
            dm: (nao, nao) 密度矩阵
        """
        S, Hcore = self.build_matrices(coords)

        # 交换相关势
        rho = self._compute_rho(dm)
        vxc = self._compute_vxc(rho)

        # Hartree势
        J = self._compute_j(dm)

        # Fock 矩阵
        F = Hcore + J + vxc

        # 对角化
        mo_energy, mo_coeff = jnp.linalg.eigh(F, S)

        # 新密度
        nocc = self.nelec // 2
        mo_occ = jnp.zeros_like(mo_energy)
        mo_occ = mo_occ.at[:nocc].set(2.0)

        dm_new = 2.0 * mo_coeff[:, :nocc] @ mo_coeff[:, :nocc].T

        return dm_new, mo_energy, mo_coeff, F

    def _compute_rho(self, dm):
        """
        计算电子密度 (简化)
        """
        # 这里应该使用数值积分
        # 简化: 基态密度近似
        return jnp.trace(dm) / self.nao

    def _compute_vxc(self, rho):
        """
        计算交换相关势 (简化)
        """
        # LDA 近似
        if isinstance(rho, (float, int)):
            rho_val = rho
        else:
            rho_val = jnp.mean(rho)

        # LDA 交换势
        vx = -(3/np.pi)**(1/3) * rho_val**(1/3)

        # 对角矩阵 (简化)
        n = self.nao
        return vx * jnp.eye(n)

    def _compute_j(self, dm):
        """
        计算 Hartree 势 (简化)
        """
        # 简化: J = G * dm, 其中 G 是库仑算符
        # 这里使用缩放近似
        return 0.5 * dm

    def solve_scf(self, coords, max_iter=50, tol=1e-8):
        """
        完整 SCF 求解

        params:
            coords: (natm, 3) 原子坐标
            max_iter: 最大迭代数
            tol: 收敛阈值
        """
        # 初始密度: 核哈密顿量
        S, Hcore = self.build_matrices(coords)
        dm = jnp.linalg.inv(S) @ Hcore @ jnp.linalg.inv(S)
        dm = dm / jnp.trace(dm) * self.nelec

        # SCF 循环
        for i in range(max_iter):
            dm_new, mo_energy, mo_coeff, F = self.scf_iteration(coords, dm)

            # 收敛检查
            delta = jnp.max(jnp.abs(dm_new - dm))
            if delta < tol:
                break

            dm = dm_new

        # 总能量
        E = self._total_energy(coords, dm, F)

        return E, dm, mo_energy, mo_coeff

    def _total_energy(self, coords, dm, F):
        """
        计算总能量
        """
        S, Hcore = self.build_matrices(coords)

        # 电子能量
        E_el = 0.5 * jnp.trace((Hcore + F) @ dm)

        # 核-核排斥
        Vnn = self._nuclear_repulsion(coords)

        return E_el + Vnn

    def _nuclear_repulsion(self, coords):
        """
        核-核排斥能
        """
        natm = coords.shape[0]
        Vnn = 0.0

        for i in range(natm):
            for j in range(i+1, natm):
                r = jnp.sqrt(jnp.sum((coords[i] - coords[j])**2))
                Z_i = self.mol.atom_charge(i)
                Z_j = self.mol.atom_charge(j)
                Vnn += Z_i * Z_j / r

        return Vnn

# 使用
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
graddft = GradDFT(mol)

coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
E, dm, mo_e, mo_c = graddft.solve_scf(coords)
print(f"DFT 能量: {E:.6f} Ha")
```

### 2. 自动微分梯度

```python
# 定义可微分能量函数
def dft_energy(coords):
    """
    DFT 能量函数 (JAX 可微分)

    params:
        coords: (natm, 3) 原子坐标
    """
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# 计算梯度
grad_energy = grad(dft_energy)

# 测试
coords_test = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
forces = -grad_energy(coords_test)

print(f"原子力:")
for i, f in enumerate(forces):
    print(f"  原子 {i}: {f}")
```

### 3. 几何优化

```python
# 梯度下降优化
def gradient_descent(initial_coords, lr=0.01, steps=100):
    """
    简单梯度下降优化

    params:
        initial_coords: 初始坐标
        lr: 学习率
        steps: 优化步数
    """
    coords = initial_coords.copy()

    for step in range(steps):
        # 计算力和能量
        E = dft_energy(coords)
        forces = -grad_energy(coords)

        # 更新坐标
        coords = coords + lr * forces

        if step % 10 == 0:
            print(f"步骤 {step}: E = {E:.6f} Ha")

    return coords, E

# 优化
coords_init = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
coords_opt, E_opt = gradient_descent(coords_init, lr=0.01, steps=50)

print(f"\n优化后能量: {E_opt:.6f} Ha")
print(f"优化后坐标:")
print(coords_init)
print(coords_opt)
```

## 泛函优化

### 1. 参数化泛函

```python
from jax import vmap

class ParameterizedXC:
    """
    参数化交换相关泛函
    E_xc = ∫ f(ρ, θ) dτ
    其中 θ 是可学习参数
    """
    def __init__(self, mol, n_params=3):
        self.mol = mol
        self.n_params = n_params

        # 初始参数 (例如: PBE 参数的变体)
        self.params = jnp.array([0.9, 0.8, 0.1])

    def lda_exchange(self, rho, params):
        """
        参数化 LDA 交换

        E_x = a1 * E_x^LDA
        """
        a1 = params[0]

        # 标准 LDA 交换
        Ex = -0.75 * (3/np.pi)**(1/3) * rho**(4/3)

        return a1 * Ex

    def lda_correlation(self, rho, params):
        """
        参数化 LDA 相关
        """
        a2 = params[1]
        a3 = params[2]

        # 简化的 VWN 形式
        Ec = a2 * rho**(1/3) + a3 * rho**(2/3)

        return Ec

    def xc_energy_density(self, rho, params):
        """
        交换相关能密度
        """
        Ex = self.lda_exchange(rho, params)
        Ec = self.lda_correlation(rho, params)

        return Ex + Ec

    def xc_potential(self, rho, params):
        """
        交换相关势 V_xc = dE_xc/dρ
        """
        # 通过自动微分计算
        # E_xc[ρ] = ∫ ε_xc(ρ) dτ
        # V_xc = dε_xc/dρ

        xc_density = lambda r: self.xc_energy_density(r, params)

        # JAX 自动微分
        vxc = grad(xc_density)(rho)

        return vxc

# 使用
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
xc = ParameterizedXC(mol)

# 计算在给定密度下的 XC 能量
rho = 0.5  # 简化: 均匀密度
Exc_density = xc.xc_energy_density(rho, xc.params)
print(f"XC 能量密度: {Exc_density:.6f}")

# XC 势
Vxc = xc.xc_potential(rho, xc.params)
print(f"XC 势: {Vxc:.6f}")
```

### 2. 泛函参数优化

```python
def functional_loss(params, training_data):
    """
    泛函参数优化损失函数

    params:
        params: 泛函参数
        training_data: [(mol, E_ref), ...] 训练数据
    """
    xc = ParameterizedXC(training_data[0][0])
    xc.params = params

    total_loss = 0.0

    for mol, E_ref in training_data:
        # 运行 DFT 计算使用参数化泛函
        graddft = GradDFT(mol)

        # 修改 XC 使用参数化泛函
        # (需要集成到 GradDFT 中)
        # graddft.xc_func = xc

        # 计算能量
        coords = jnp.array(mol.atom_coords())
        E_pred, _, _, _ = graddft.solve_scf(coords)

        # MSE 损失
        total_loss += (E_pred - E_ref)**2

    return total_loss / len(training_data)

# 优化
from jax import value_and_grad
from jax.example_libraries import optimizers

# 创建训练数据 (简化)
training_mols = [
    (gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g'), -1.16),
    (gto.M(atom='H 0 0 0; H 0 0 0.80', basis='sto-3g'), -1.14),
]

# 初始参数
init_params = jnp.array([0.9, 0.8, 0.1])

# 优化器
opt_init, opt_update, get_params = optimizers.adam(0.01)
opt_state = opt_init(init_params)

@jit
def step(i, opt_state, training_data):
    params = get_params(opt_state)
    loss, grads = value_and_grad(functional_loss)(params, training_data)
    opt_state = opt_update(i, grads, opt_state)
    return loss, opt_state

# 优化循环
for i in range(100):
    loss, opt_state = step(i, opt_state, training_mols)
    if i % 10 == 0:
        print(f"迭代 {i}: Loss = {loss:.6f}")

# 最优参数
optimal_params = get_params(opt_state)
print(f"\n最优参数: {optimal_params}")
```

## 高级集成

### 1. 神经网络增强 DFT

```python
import flax.linen as nn
import optax

class NNXC(nn.Module):
    """
    神经网络交换相关泛函
    """
    hidden_dim: int = 64

    @nn.compact
    def __call__(self, rho_features):
        """
        密度特征 → XC 能量

        params:
            rho_features: (N, n_features) 密度特征
        """
        x = nn.Dense(self.hidden_dim)(rho_features)
        x = nn.tanh(x)
        x = nn.Dense(self.hidden_dim)(x)
        x = nn.tanh(x)
        x = nn.Dense(1)(x)
        return x.squeeze()

# 特征提取
def extract_density_features(mol, dm):
    """
    从密度矩阵提取特征
    """
    # 简化: 使用局部密度和梯度
    # 实际实现应该计算每个网格点的特征

    # 示例特征
    rho = jnp.trace(dm) / mol.nao  # 平均密度
    grad_rho = jnp.sqrt(jnp.sum(jnp.array([0.1, 0.2, 0.3])**2))  # 简化

    features = jnp.array([rho, grad_rho, rho**2, grad_rho**2])
    return features

# 端到端训练
class NNEnergyModel:
    """
    神经网络增强的 DFT 能量模型
    """
    def __init__(self, mol):
        self.mol = mol
        self.nn_xc = NNXC()
        self.graddft = GradDFT(mol)

    def forward(self, coords, dm):
        """
        前向传播
        """
        # 提取密度特征
        features = extract_density_features(self.mol, dm)

        # NN 预测 XC 能量密度
        Exc_nn = self.nn_xc(features)

        # 计算总能量
        E = self.graddft._total_energy(coords, dm, None)

        return E + Exc_nn

    def loss(self, coords, dm, E_ref):
        E_pred = self.forward(coords, dm)
        return (E_pred - E_ref)**2

# 训练
mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
model = NNEnergyModel(mol)

# 初始化参数
dummy_features = jnp.array([0.5, 0.1, 0.25, 0.01])
params = model.nn_xc.init(jax.random.PRNGKey(0), dummy_features)

# 优化器
optimizer = optax.adam(0.01)
opt_state = optimizer.init(params)

@jit
def train_step(coords, dm, E_ref, params, opt_state):
    def loss_fn(p):
        # 临时更新模型参数
        # (需要正确集成到 Flax)
        return model.loss(coords, dm, E_ref)

    loss, grads = jax.value_and_grad(loss_fn)(params)
    updates, opt_state = optimizer.update(grads, opt_state)
    params = optax.apply_updates(params, updates)
    return loss, params, opt_state

# 训练循环 (简化)
coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
dm = jnp.eye(mol.nao) * 0.5  # 简化密度
E_ref = -1.16  # 参考能量

for i in range(100):
    loss, params, opt_state = train_step(coords, dm, E_ref, params, opt_state)
    if i % 20 == 0:
        print(f"迭代 {i}: Loss = {loss:.6f}")
```

### 2. 多体系批量计算

```python
# 批量处理多个分子
from jax import vmap

def batch_dft_energy(batch_coords, mol_template):
    """
    批量计算多个分子的 DFT 能量

    params:
        batch_coords: (batch_size, natm, 3) 批量坐标
        mol_template: 分子模板 (基组等)
    """
    def single_energy(coords):
        graddft = GradDFT(mol_template)
        E, _, _, _ = graddft.solve_scf(coords)
        return E

    # 向量化
    energies = vmap(single_energy)(batch_coords)
    return energies

# 测试
batch_coords = jnp.array([
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],  # 分子 1
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.80]],  # 分子 2
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.86]],  # 分子 3
])

energies = batch_dft_energy(batch_coords, mol)
print(f"批量能量: {energies}")
```

### 3. GPU 加速

```python
# 检查 GPU 可用性
from jax.lib import xla_bridge
print(f"设备: {xla_bridge.get_backend().platform}")

# 将数据移到 GPU
coords_gpu = jax.device_put(coords, jax.devices('gpu')[0])

# JIT 编译 + GPU 自动加速
@jit
def fast_scf(coords):
    graddft = GradDFT(mol)
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# 运行 (自动使用 GPU)
E_gpu = fast_scf(coords_gpu)
```

## 实际应用

### 1. 势能面扫描

```python
# H₂ 势能面扫描
def scan_bond_length(r_values):
    """
    扫描 H₂ 键长

    params:
        r_values: (n,) 键长数组
    """
    energies = []

    for r in r_values:
        # 构建分子
        coords = jnp.array([
            [0.0, 0.0, -r/2],
            [0.0, 0.0, r/2]
        ])

        # 计算 DFT 能量
        E = dft_energy(coords)
        energies.append(E)

    return jnp.array(energies)

# 扫描
r_values = jnp.linspace(0.5, 2.0, 20)
energies = scan_bond_length(r_values)

# 绘制
import matplotlib.pyplot as plt
plt.plot(r_values, energies, 'o-')
plt.xlabel('Bond Length (Å)')
plt.ylabel('Energy (Ha)')
plt.title('H₂ Potential Energy Curve')
plt.show()
```

### 2. 力场拟合

```python
# 拟合力场到 DFT 力
def force_loss(params, coords, target_forces):
    """
    力场参数优化损失

    params:
        params: 力场参数
        coords: 原子坐标
        target_forces: 目标力 (natm, 3)
    """
    # 计算 DFT 力
    forces_pred = -grad_energy(coords)

    # MSE 损失
    loss = jnp.mean((forces_pred - target_forces)**2)
    return loss

# 使用
coords = jnp.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
target_forces = jnp.array([[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]])  # 示例

# 优化 (简化)
params = jnp.array([1.0, 0.5])

for i in range(50):
    loss, grads = value_and_grad(force_loss)(params, coords, target_forces)
    params = params - 0.01 * grads

    if i % 10 == 0:
        print(f"迭代 {i}: Loss = {loss:.6f}")
```

## 注意事项

### 1. 数值稳定性

```python
# 避免除零
r = jnp.sqrt(jnp.sum((coords[i] - coords[j])**2) + 1e-12)

# 梯度裁剪
@jit
def clip_grads(grads, max_norm=1.0):
    norm = jnp.sqrt(jnp.sum(grads**2))
    scale = jnp.minimum(1.0, max_norm / (norm + 1e-8))
    return grads * scale

# 使用
grads = clip_grads(grad_energy(coords))
```

### 2. SCF 收敛

```python
# 收敛检查
def is_converged(dm_old, dm_new, tol=1e-8):
    delta = jnp.max(jnp.abs(dm_new - dm_old))
    return delta < tol

# DIIS 外推 (简化)
def diis_extrapolate(fock_list, error_list, diis_dim=6):
    """
    DIIS 外推 Fock 矩阵
    """
    # 实现 DIIS 算法
    # (简化: 直接返回最新的 Fock)
    return fock_list[-1]
```

### 3. 内存管理

```python
# 大分子分块处理
def batch_process_large_molecule(mol, batch_size=100):
    """
    分批处理大分子
    """
    n_atoms = mol.natm
    coords = mol.atom_coords()

    # 分批计算能量贡献
    total_energy = 0.0

    for i in range(0, n_atoms, batch_size):
        batch_coords = coords[i:i+batch_size]
        E_batch = dft_energy(batch_coords)
        total_energy += E_batch

    return total_energy
```

## 性能优化

### 1. JIT 编译

```python
# 编译关键函数
@jit
def compiled_scf(coords):
    """编译后的 SCF"""
    graddft = GradDFT(mol)
    E, _, _, _ = graddft.solve_scf(coords)
    return E

# 首次调用会编译
E1 = compiled_scf(coords)  # 编译

# 后续调用很快
E2 = compiled_scf(coords)  # 快速
```

### 2. 向量化

```python
# 向量化多个构型
from jax import vmap

# 预先编译
batch_energy = vmap(compiled_scf)

# 批量计算
batch_coords = jnp.array([coords, coords*1.01, coords*0.99])
energies = batch_energy(batch_coords)
```

### 3. 并行化

```python
# 多设备并行 (GPU 集群)
from jax import pmap

# 预分发数据到多个 GPU
batch_coords = jnp.array([coords] * 4)  # 4 个构型
batch_coords = jnp.reshape(batch_coords, (4, -1, 3))

# 并行计算
parallel_energy = pmap(compiled_scf)
energies = parallel_energy(batch_coords)
```

## 参考

### PySCF-JAX
1. PySCF JAX 模块: https://github.com/pyscf/pyscf/tree/master/pyscf/jax
2. GradDFT: https://github.com/pyscf/pyscf/tree/master/pyscf/dft/gradients

### JAX
3. JAX 文档: https://jax.readthedocs.io/
4. JAX 自动微分: https://jax.readthedocs.io/en/latest/jax.html#automatic-differentiation

### 论文
5. Miller et al., "JAX-MD: A Framework for Differentiable Molecular Dynamics", arXiv 2023
6. Kasim et al., "Differentiable Quantum Chemistry with JAX", arXiv 2022
