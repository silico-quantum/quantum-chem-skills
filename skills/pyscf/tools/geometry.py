#!/usr/bin/env python3
"""
Geometry Optimization - 几何优化

支持:
    - 能量最小化
    - 过渡态搜索
    - 频率计算（验证极小点/过渡态）

用法:
    python geometry.py <xyz_file> <method> [basis] [functional]
    python geometry.py water_guess.xyz RKS cc-pvdz b3lyp
    python geometry.py ts_guess.xyz Saddle cc-pvdz
"""

import sys
import numpy as np

from pyscf import gto, dft, scf, geomopt, hessian


def read_xyz(xyz_file):
    """从 XYZ 文件读取分子"""
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    return symbols, np.array(coords)


def write_xyz(symbols, coords, filename, energy=None, comment=None):
    """写入 XYZ 文件"""
    natoms = len(symbols)
    with open(filename, 'w') as f:
        if comment:
            f.write(f"{natoms}\n{comment}\n")
        elif energy:
            f.write(f"{natoms}\nenergy = {energy:.10f}\n")
        else:
            f.write(f"{natoms}\noptimized structure\n")
        for s, c in zip(symbols, coords):
            f.write(f"{s:2s} {c[0]:15.10f} {c[1]:15.10f} {c[2]:15.10f}\n")


def run_geomopt(xyz_file, method='RKS', basis='cc-pvdz', functional='pbe',
                conv_grad=1e-4, conv_energy=1e-6, max_steps=100):
    """
    几何优化（能量最小化）
    
    Args:
        xyz_file: 初始 XYZ 文件
        method: 'RKS', 'UKS', 'RHF', 'UHF'
        basis: 基组
        functional: DFT 泛函（用于 RKS/UKS）
        conv_grad: 力收敛标准 (Hartree/Angstrom)
        conv_energy: 能量收敛标准 (Hartree)
        max_steps: 最大优化步数
    
    Returns:
        dict: 优化后的分子对象和能量
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    # 构建 SCF 对象
    if method == 'RKS':
        mf = dft.RKS(mol, xc=functional)
    elif method == 'UKS':
        mf = dft.UKS(mol, xc=functional)
    elif method == 'RHF':
        mf = scf.RHF(mol)
    elif method == 'UHF':
        mf = scf.UHF(mol)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # 几何优化
    opt = geomopt.GeometryOptimizer(mf)
    opt.conv_grad = conv_grad
    opt.conv_energy = conv_energy
    opt.max_steps = max_steps
    
    print("开始几何优化...")
    mol_opt = opt.run()
    
    print("\n" + "=" * 50)
    print("优化完成")
    print("=" * 50)
    print(f"总能量: {mol_opt.energy:.10f} Hartree")
    print("\n优化后坐标 (Angstrom):")
    for i in range(mol_opt.natm):
        print(f"  {mol_opt.atom_symbol(i):2s}: {mol_opt.atom_coord(i)}")
    
    # 保存
    output_file = xyz_file.replace('.xyz', '_opt.xyz')
    write_xyz(
        [mol_opt.atom_symbol(i) for i in range(mol_opt.natm)],
        mol_opt.atom_coords(),
        output_file,
        energy=mol_opt.energy
    )
    print(f"\n已保存到: {output_file}")
    
    return {'mol_opt': mol_opt, 'output_file': output_file}


def run_saddle(xyz_file, method='RKS', basis='cc-pvdz', functional='pbe',
                conv_grad=1e-4):
    """
    过渡态优化（寻找一阶鞍点）
    
    Uses dimer method 或 Berny algorithm
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    if method == 'RKS':
        mf = dft.RKS(mol, xc=functional)
    elif method == 'RHF':
        mf = scf.RHF(mol)
    else:
        raise ValueError("过渡态目前仅支持 RKS/RHF")
    
    # 过渡态优化
    opt = geomopt.GeometryOptimizer(mf, 'saddle')
    opt.conv_grad = conv_grad
    
    print("开始过渡态搜索...")
    mol_ts = opt.run()
    
    print(f"\n过渡态能量: {mol_ts.energy:.10f} Hartree")
    
    output_file = xyz_file.replace('.xyz', '_ts.xyz')
    write_xyz(
        [mol_ts.atom_symbol(i) for i in range(mol_ts.natm)],
        mol_ts.atom_coords(),
        output_file,
        energy=mol_ts.energy
    )
    
    return {'mol_opt': mol_ts, 'output_file': output_file}


def run_freq(xyz_file, method='RKS', basis='cc-pvdz', functional='pbe'):
    """
    频率计算（验证优化结构）
    
    返回:
        - 所有频率 (cm⁻¹)
        - 虚频数量（0 = 极小点, 1 = 过渡态）
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    if method == 'RKS':
        mf = dft.RKS(mol, xc=functional).run()
    elif method == 'RHF':
        mf = scf.RHF(mol).run()
    else:
        raise ValueError("频率计算仅支持 RKS/RHF")
    
    # 计算 Hessian
    if method == 'RKS':
        hess = hessian.RKS(mf).run()
    else:
        hess = hessian.RHF(mf).run()
    
    # 频率
    freq = hess.kernel()
    freq_cm = np.sqrt(np.abs(freq)) * 5140.0  # atomic unit -> cm⁻¹
    
    # 统计虚频
    n_imag = np.sum(freq_cm < 0)
    
    print("\n" + "=" * 50)
    print("频率分析")
    print("=" * 50)
    print(f"虚频数量: {n_imag}")
    if n_imag == 0:
        print("结论: 极小点（稳定的几何构型）")
    else:
        print("结论: 鞍点（过渡态）")
    
    print(f"\n前 10 个频率 (cm⁻¹):")
    sorted_freq = sorted(freq_cm)
    for i, f in enumerate(sorted_freq[:10]):
        tag = " (虚频)" if f < 0 else ""
        print(f"  {i+1:3d}: {f:10.4f}{tag}")
    
    return {'freq': freq_cm, 'n_imag': n_imag}


def scan_bond_length(xyz_file, atom_pair, r_range, n_steps=10,
                     method='RKS', functional='pbe', basis='cc-pvdz'):
    """
    扫描指定键长（1D PES）
    
    Args:
        xyz_file: 初始结构
        atom_pair: tuple, e.g., (0, 1) 表示第 0 和第 1 个原子之间的键
        r_range: (r_min, r_max) 扫描范围 (Angstrom)
        n_steps: 扫描点数
    
    Returns:
        dict: {'r': 键长数组, 'E': 能量数组}
    """
    symbols, coords = read_xyz(xyz_file)
    i, j = atom_pair
    
    # 计算初始键长
    r_init = np.linalg.norm(coords[i] - coords[j])
    
    print(f"键长扫描: 原子 {i}-{j}")
    print(f"初始键长: {r_init:.4f} Angstrom")
    print(f"范围: {r_range[0]:.2f} - {r_range[1]:.2f} Angstrom, {n_steps} 步")
    
    r_values = np.linspace(r_range[0], r_range[1], n_steps)
    energies = []
    
    for idx, r in enumerate(r_values):
        # 修改坐标：保持其他原子不动，只改变键长
        new_coords = coords.copy()
        direction = (coords[j] - coords[i]) / np.linalg.norm(coords[j] - coords[i])
        new_coords[j] = new_coords[i] + r * direction
        
        atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                              for s, c in zip(symbols, new_coords)])
        
        mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0)
        
        if method == 'RKS':
            mf = dft.RKS(mol, xc=functional).run()
        else:
            mf = scf.RHF(mol).run()
        
        energies.append(mf.e_tot)
        print(f"  Step {idx+1}/{n_steps}: r={r:.4f} E={mf.e_tot:.10f}")
    
    result = {'r': r_values, 'E': np.array(energies)}
    
    # 保存扫描数据
    np.savetxt(xyz_file.replace('.xyz', '_scan.dat'), 
               np.column_stack([result['r'], result['E']]),
               header="r(Angstrom) E(Hartree)")
    
    return result


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\n示例:")
        print("  python geometry.py water_guess.xyz RKS cc-pvdz b3lyp")
        print("  python geometry.py ts_guess.xyz Saddle cc-pvdz")
        print("  python geometry.py opt.xyz Freq RKS cc-pvdz b3lyp")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    task = sys.argv[2].lower()
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    functional = sys.argv[4] if len(sys.argv) > 4 else 'pbe'
    
    if task in ['opt', 'min', 'minimize']:
        run_geomopt(xyz_file, 'RKS', basis, functional)
    elif task in ['ts', 'saddle', 'transition']:
        run_saddle(xyz_file, 'RKS', basis, functional)
    elif task in ['freq', 'frequency']:
        run_freq(xyz_file, 'RKS', basis, functional)
    else:
        print(f"Unknown task: {task}")
        print("Available: opt, ts, freq")
        sys.exit(1)


if __name__ == '__main__':
    main()
