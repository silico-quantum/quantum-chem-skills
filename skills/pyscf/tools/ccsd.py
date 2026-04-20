#!/usr/bin/env python3
"""
CCSD / CCSD(T) - 耦合簇单双（含三次项）

支持:
    - RCCSD (闭壳层)
    - UCCSD (开壳层)
    - CCSD(T) 修正
    - 密度拟合加速

用法:
    python ccsd.py <xyz_file> [basis] [nroots]
    python ccsd.py water.xyz cc-pvdz
"""

import sys
import numpy as np

from pyscf import gto, scf, cc


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


def run_ccsd(xyz_file, basis='cc-pvdz', frozen=0, density_fit=False, ccdo=0):
    """
    CCSD 计算（通用入口，自动判断 R/U）
    
    Args:
        xyz_file: XYZ 文件
        basis: 基组
        frozen: frozen 电子数
        density_fit: 使用 DF 加速
        ccdo: CCSD 迭代次数 (0 = 默认收敛)
    
    Returns:
        dict: {'cc': CCSD对象, 'mf': SCF对象, 'e_tot': CCSD总能量}
    """
    symbols, coords = read_xyz(xyz_file)
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    # HF 计算
    if spin == 0:
        mf = scf.RHF(mol).run()
    else:
        mf = scf.UHF(mol).run()
    
    print(f"HF 完成: E = {mf.e_tot:.10f}")
    
    # CCSD
    if spin == 0:
        mycc = cc.CCSD(mf, frozen=frozen)
    else:
        mycc = cc.UCCSD(mf, frozen=frozen)
    
    if density_fit:
        print("使用密度拟合加速...")
        # DF-CCSD 需要特殊设置
        from pyscf import df
        mycc = mycc.density_fit(auxbasis='cc-pvdz-ri')
    
    mycc.max_cycle = 100 if ccdo == 0 else ccdo
    mycc.verbose = 5
    mycc.kernel()
    
    print("\n" + "=" * 50)
    print("CCSD 结果")
    print("=" * 50)
    print(f"E(CCSD) = {mycc.e_tot:.10f} Hartree")
    print(f"E(corr) = {mycc.e_tot - mf.e_tot:.10f} Hartree")
    print(f"T1 shape: {mycc.t1.shape}")
    print(f"T2 shape: {mycc.t2.shape}")
    
    return {'cc': mycc, 'mf': mf, 'mol': mol}


def run_ccsd_t(xyz_file, basis='cc-pvdz', frozen=0):
    """
    CCSD(T) - CCSD + perturbative triple correction
    CCSD(T) 的计算成本约为 CCSD 的 3-4 倍
    """
    result = run_ccsd(xyz_file, basis, frozen)
    mycc = result['cc']
    mf = result['mf']
    
    # CCSD(T) correction
    eccsd_t = mycc.ccsd_t()
    
    e_total = mycc.e_tot + eccsd_t
    
    print("\n" + "=" * 50)
    print("CCSD(T) 结果")
    print("=" * 50)
    print(f"E(CCSD(T)) = {e_total:.10f} Hartree")
    print(f"ΔE(T)      = {eccsd_t:.10f} Hartree")
    
    result['e_tot'] = e_total
    result['eccsd_t'] = eccsd_t
    
    return result


def run_rccsd(xyz_file, basis='cc-pvdz', frozen=0):
    """
    限制性闭壳层 RCCSD
    """
    symbols, coords = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=0, verbose=4)
    
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf, frozen=frozen).run()
    
    return {'cc': mycc, 'mf': mf, 'mol': mol}


def run_uccsd(xyz_file, basis='cc-pvdz', frozen=0):
    """
    非限制性开壳层 UCCSD
    """
    symbols, coords = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    mf = scf.UHF(mol).run()
    mycc = cc.UCCSD(mf, frozen=frozen).run()
    
    print(f"E(UCCSD) = {mycc.e_tot:.10f}")
    
    return {'cc': mycc, 'mf': mf, 'mol': mol}


def get_one_particle_density_matrix(result):
    """
    计算单粒子密度矩阵 (1PDM)
    用于后续分析
    """
    cc = result['cc']
    dm1 = cc.make_rdm1()
    print(f"1PDM shape: {dm1.shape}")
    return dm1


def get_two_particle_density_matrix(result):
    """
    计算双粒子密度矩阵 (2PDM)
    用于后 HF 相关分析
    """
    cc = result['cc']
    dm2 = cc.make_rdm2()
    print(f"2PDM shape: {dm2.shape}")
    return dm2


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\n示例:")
        print("  python ccsd.py water.xyz cc-pvdz")
        print("  python ccsd.py mol.xyz 6-31G*")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    basis = sys.argv[2] if len(sys.argv) > 2 else 'cc-pvdz'
    include_t = '-t' in sys.argv or '--ccsdt' in sys.argv
    
    symbols, _ = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    print(f"分子: {' '.join(symbols)}, 基组: {basis}, 自旋: {spin}")
    
    if include_t:
        result = run_ccsd_t(xyz_file, basis)
    else:
        result = run_ccsd(xyz_file, basis)


if __name__ == '__main__':
    main()
