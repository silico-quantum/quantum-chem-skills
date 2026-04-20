#!/usr/bin/env python3
"""
MP2 - Moller-Plesset 二阶微扰理论

支持:
    - RMP2 (闭壳层)
    - UMP2 (开壳层)
    - 密度拟合加速

用法:
    python mp2.py <xyz_file> [basis] [frozen]
    python mp2.py water.xyz cc-pvdz
"""

import sys
import numpy as np

from pyscf import gto, scf, mp


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


def run_rmp2(xyz_file, basis='cc-pvdz', frozen=0, density_fit=False):
    """
    闭壳层 RMP2 计算
    
    Args:
        xyz_file: XYZ 文件路径
        basis: 基组名
        frozen: frozen core 电子数（默认 0 = 全部电子参与）
        density_fit: 是否使用密度拟合加速
    
    Returns:
        dict: {'mp': MP2对象, 'mf': SCF对象, 'mol': 分子对象}
    """
    symbols, coords = read_xyz(xyz_file)
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    # HF 计算
    mf = scf.RHF(mol).run()
    
    # MP2 计算
    if density_fit:
        mmp = mp.MP2(mf).density_fit(auxbasis='cc-pvdz-ri').run()
    else:
        if frozen > 0:
            mmp = mp.MP2(mf, frozen=frozen).run()
        else:
            mmp = mp.MP2(mf).run()
    
    print("\n" + "=" * 50)
    print("MP2 结果")
    print("=" * 50)
    print(f"E(MP2)     = {mmp.e_tot:.10f} Hartree")
    print(f"E(corr)    = {mmp.e_corr:.10f} Hartree")
    print(f"E(SCF)     = {mmp.e_mp2 + mmp.e_corr - mmp.e_mp2 + mf.e_tot:.10f}")
    
    print(f"\nT2 amplitudes shape: {mmp.t2.shape}")
    
    return {'mp': mmp, 'mf': mf, 'mol': mol}


def run_ump2(xyz_file, basis='cc-pvdz', frozen=0):
    """
    开壳层 UMP2 计算
    
    Args:
        xyz_file: XYZ 文件
        basis: 基组
        frozen: frozen core 电子数
    """
    symbols, coords = read_xyz(xyz_file)
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    mf = scf.UHF(mol).run()
    
    if frozen > 0:
        mmp = mp.UMP2(mf, frozen=frozen).run()
    else:
        mmp = mp.UMP2(mf).run()
    
    print("\nUMP2 结果:")
    print(f"E(corr) = {mmp.e_corr:.10f} Hartree")
    print(f"T2 shape (aa): {mmp.t2[0].shape}")
    print(f"T2 shape (bb): {mmp.t2[1].shape}")
    
    return {'mp': mmp, 'mf': mf, 'mol': mol}


def run_direct_mp2(xyz_file, basis='cc-pvdz'):
    """
    直接 MP2（基于 RI-MP2，省内存）
    使用 density fitting 近似双电子积分
    """
    symbols, coords = read_xyz(xyz_file)
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin)
    
    mf = scf.RHF(mol).run()
    
    # 使用 density fitting
    from pyscf import df
    mmp = mp.MP2(mf).density_fit(auxbasis='cc-pvdz-ri').run()
    
    return {'mp': mmp, 'mf': mf}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\n示例:")
        print("  python mp2.py water.xyz cc-pvdz")
        print("  python mp2.py molecule.xyz 6-31G*")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    basis = sys.argv[2] if len(sys.argv) > 2 else 'cc-pvdz'
    frozen = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    
    symbols, _ = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    print(f"分子: {' '.join(symbols)}, 基组: {basis}, 自旋: {spin}")
    
    if spin == 0:
        result = run_rmp2(xyz_file, basis, frozen)
    else:
        result = run_ump2(xyz_file, basis, frozen)


if __name__ == '__main__':
    main()
