#!/usr/bin/env python3
"""
CASSCF / CASCI / CASPT2 / NEVPT2 - 多参考态方法

支持:
    - CASSCF (完全活化空间 SCF)
    - CASCI
    - CASPT2 (二级微扰修正)
    - NEVPT2 (N-Electron Valence State Perturbation Theory)

用法:
    python cascf.py <xyz_file> <ncas> <nelec> [basis] [functional]
    python cascf.py mol.xyz 8 10 cc-pvdz pbe
    python cascf.py o2_triplet.xyz 8 12 sto-3g
"""

import sys
import numpy as np

from pyscf import gto, scf, dft, mcscf


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


def auto_nelec_from_atoms(symbols, charge=0):
    """
    自动计算活化电子数
    
    常见活性空间:
    - 1,4-丁二烯 (C4H6): 4电子/4轨道
    - 苯 (C6H6): 6电子/6轨道
    - 甲醛 (CH2O): 6电子/4轨道 (C=O pi, O lone pair)
    - 萘 (C10H8): 10电子/10轨道
    
    Returns:
        tuple: (ncas, nelec)
    """
    # 简单估算：每个原子 1 电子参与活化
    # 这只是一个估计，实际需要化学直觉
    nelec = sum(gto.charge(s) for s in symbols) - charge
    ncas = nelec  # 保守估计
    return ncas, nelec


def guess_active_space(mf, ncas, nelec, loc_method='minao'):
    """
    猜测活化空间轨道
    
    Args:
        mf: SCF 对象
        ncas: 活化轨道数
        nelec: 活化电子数
        loc_method: 'minao', 'meta_lowdin', 'ER', ' Boys'
    
    Returns:
        array: 排序后的分子轨道系数
    """
    from pyscf import lo
    
    # 本地化轨道
    mo = lo.orth_ao(mf.mol, loc_method)
    mo_cas = mo[:, :ncas]  # 选择前 ncas 个
    
    return mo_cas


def run_casscf(xyz_file, ncas, nelec, basis='cc-pvdz', functional='pbe', 
               nstates=1, freeze_core=True, localization='minao'):
    """
    运行 CASSCF 计算
    
    Args:
        xyz_file: XYZ 文件
        ncas: 活化轨道数
        nelec: 活化电子数
        basis: 基组
        functional: DFT 泛函（用于初始轨道）
        nstates: 多态 CASSCF 的态数量
        freeze_core: 是否冻结内层轨道
        localization: 轨道本地化方法
    
    Returns:
        dict: {'mc': CASSCF对象, 'mf': SCF对象, 'mo': 轨道系数}
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    # 确定自旋
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    # DFT 计算（提供初始轨道）
    if functional:
        if spin == 0:
            mf = dft.RKS(mol, xc=functional).run()
        else:
            mf = dft.UKS(mol, xc=functional).run()
    else:
        if spin == 0:
            mf = scf.RHF(mol).run()
        else:
            mf = scf.UHF(mol).run()
    
    print(f"\n活化空间: {nelec}电子 / {ncas}轨道")
    print(f"DFT 完成: E = {mf.e_tot:.10f}")
    
    # CASSCF
    mc = mcscf.CASSCF(mf, ncas, nelec)
    
    if nstates > 1:
        print(f"多态 CASSCF: {nstates} 态")
        mc = mcscf.CASSCF(mf, ncas, nelec, nstates)
    
    # 冻结内层轨道
    if freeze_core:
        mc.frozen = 1  # 冻结 1s 轨道（用于重原子）
    
    # 运行
    mc.verbose = 5
    e_casscf = mc.run()
    
    print("\n" + "=" * 50)
    print(f"CASSCF({ncas},{nelec}) 结果")
    print("=" * 50)
    print(f"E(CASSCF) = {e_casscf:.10f} Hartree")
    print(f"E(corr)   = {e_casscf - mf.e_tot:.10f} Hartree")
    
    return {'mc': mc, 'mf': mf, 'mol': mol, 'e_tot': e_casscf}


def run_casci(xyz_file, ncas, nelec, basis='cc-pvdz', functional='pbe'):
    """
    CASCI（固定轨道，只优化 CI 系数，比 CASSCF 快）
    """
    symbols, coords = read_xyz(xyz_file)
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    if functional:
        if spin == 0:
            mf = dft.RKS(mol, xc=functional).run()
        else:
            mf = dft.UKS(mol, xc=functional).run()
    else:
        mf = scf.RHF(mol).run()
    
    mc = mcscf.CASCI(mf, ncas, nelec).run()
    
    print(f"\nCASCI({ncas},{nelec}) E = {mc.e_tot:.10f}")
    
    return {'mc': mc, 'mf': mf}


def run_caspt2(xyz_file, ncas, nelec, basis='cc-pvdz', functional='pbe',
               root=0, shift=0.0, imacro=0):
    """
    CASPT2（二级微扰理论，对 CASSCF 能量修正）
    
    Args:
        xyz_file: XYZ 文件
        ncas, nelec: 活化空间
        root: 修正的态（0 = 基态）
        shift: level shift（防 intruder state）
        imacro: 冻结 macro-core
    """
    # 先跑 CASSCF
    result = run_casscf(xyz_file, ncas, nelec, basis, functional)
    mc = result['mc']
    
    # CASPT2
    mcpt2 = mcscf.CASPT2(mc)
    mcpt2.shift = shift  # intruder state avoidance
    mcpt2.verbose = 5
    
    if root > 0:
        mcpt2.frozen = 1
    
    e_pt2 = mcpt2.kernel()
    
    print("\n" + "=" * 50)
    print("CASPT2 结果")
    print("=" * 50)
    print(f"E(CASPT2) = {e_pt2:.10f} Hartree")
    print(f"ΔE(PT2)   = {e_pt2 - mc.e_tot[0] if hasattr(mc, 'e_tot') and isinstance(mc.e_tot, list) else e_pt2 - mc.e_tot:.10f} Hartree")
    
    result['mcpt2'] = mcpt2
    result['e_tot'] = e_pt2
    
    return result


def run_nevpt2(xyz_file, ncas, nelec, basis='cc-pvdz', functional='pbe'):
    """
    NEVPT2（N-electron valence state perturbation theory）
    比 CASPT2 更稳健，对 intruder state 不敏感
    """
    result = run_casscf(xyz_file, ncas, nelec, basis, functional)
    mc = result['mc']
    
    mnevpt = mcscf.NEVPT2(mc).run()
    
    print(f"\nNEVPT2 E = {mnevpt.e_tot:.10f}")
    
    result['mnevpt'] = mnevpt
    return result


def run_ms-caspt2(xyz_file, ncas, nelec, nstates, basis='cc-pvdz', functional='pbe',
                   shifts=None):
    """
    ms-CASPT2（多态 CASPT2，考虑态相互作用）
    
    Args:
        nstates: CASSCF 态数量
        shifts: dict, 每个态的 level shift
    """
    result = run_casscf(xyz_file, ncas, nelec, basis, functional, nstates=nstates)
    mc = result['mc']
    
    from pyscf import lib
    mcpt2 = mcscf.CASPT2(mc, energy=True).mspt2()
    
    if shifts:
        for state, s in shifts.items():
            mcpt2.level_shift[state] = s
    
    e_ms = mcpt2.kernel()
    
    print(f"\nms-CASPT2: {nstates} states")
    for i, e in enumerate(e_ms):
        print(f"  State {i}: {e:.10f}")
    
    return result


def analyze_mo_contributions(mc, top_n=10):
    """
    分析活化轨道的组成
    """
    mo_coeff = mc.mo_coeff
    
    print("\n活化空间轨道分析:")
    print("（需要人工判断哪些是 pi, sigma, n 等）")
    
    # 简单的轨道能量排序
    mo_energy = mc.mo_energy * 27.2114
    cas_space = mc.ncas
    
    print(f"\n前 {min(top_n, cas_space)} 个活化轨道能量 (eV):")
    for i in range(min(top_n, cas_space)):
        print(f"  CAS orbital {i}: {mo_energy[i]:.4f} eV")


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        print("\n示例:")
        print("  python cascf.py ethene.xyz 2 4 cc-pvdz")
        print("  python cascf.py o2_triplet.xyz 8 12 sto-3g")
        print("  python cascf.py formaldehyde.xyz 4 6 6-31g* pbe")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    ncas = int(sys.argv[2])
    nelec = int(sys.argv[3])
    basis = sys.argv[4] if len(sys.argv) > 4 else 'cc-pvdz'
    functional = sys.argv[5] if len(sys.argv) > 5 else 'pbe'
    
    print(f"CASSCF({ncas},{nelec})")
    print(f"基组: {basis}, 泛函: {functional}")
    
    result = run_casscf(xyz_file, ncas, nelec, basis, functional)


if __name__ == '__main__':
    main()
