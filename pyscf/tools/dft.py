#!/usr/bin/env python3
"""
DFT - 密度泛函理论专用工具

泛函分类:
    LDA: SVWN5
    GGA: PBE, BLYP, BP86
    meta-GGA: TPSS, SCAN
    hybrid: B3LYP, PBE0, M06, M06-2X
    range-separated: WB97X, WB97X-D, CAM-B3LYP

用法:
    python dft.py <xyz_file> <functional> [basis]
    python dft.py water.csv pbe cc-pvdz
    python dft.py mol.xyz wb97x-d3 def2-tzvp
"""

import sys
import numpy as np

from pyscf import gto, dft


# 泛函速查表
FUNCTIONAL_CATEGORIES = {
    'lda': ['svwn5', 'lda'],
    'gga': ['pbe', 'blyp', 'bp86', 'pbesol'],
    'meta_gga': ['tpss', 'scan', 'm06-l', 'tpss0'],
    'hybrid': ['b3lyp', 'pbe0', 'm06', 'm06-2x', 'b3pw91', 'wb97x', 'scan0'],
    'range_separated': ['wb97x', 'wb97x-d', 'wb97x-d3bj', 'camb3lyp', 'lc-wpbe'],
    'dispersion': ['pbe-d3', 'b3lyp-d3', 'wb97x-d3bj'],
}

# 推荐泛函速查
RECOMMENDED = {
    '常规有机分子 (平衡)': 'b3lyp',
    '常规有机分子 (快速)': 'pbe',
    '电荷转移 / CT 态': 'wb97x-d3bj',
    '弱相互作用 / 范德华': 'wb97x-d3bj',
    '激发态 (TDDFT)': 'b3lyp',
    '溶剂化': 'pbe0',
    '过渡金属': 'b3lyp-d3',
    '重原子 (含重元素)': 'wpbesol',
}


def list_functionals(category=None):
    """列出所有可用泛函"""
    if category:
        return FUNCTIONAL_CATEGORIES.get(category, [])
    else:
        all_funcs = set()
        for funcs in FUNCTIONAL_CATEGORIES.values():
            all_funcs.update(funcs)
        return sorted(all_funcs)


def recommend_functional(use_case):
    """根据用途推荐泛函"""
    return RECOMMENDED.get(use_case, 'b3lyp')


def run_dft(xyz_file, functional='pbe', basis='cc-pvdz', 
            restricted=True, spin=0, charge=0):
    """
    运行 DFT 计算
    
    Args:
        xyz_file: XYZ 文件
        functional: 泛函名
        basis: 基组
        restricted: True = RKS, False = UKS
        spin: 自旋多重度
        charge: 电荷
    
    Returns:
        dict: {'mf': SCF对象, 'mol': 分子对象}
    """
    # 读取 XYZ
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    coords = np.array(coords)
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=4)
    
    if restricted:
        mf = dft.RKS(mol, xc=functional).run()
    else:
        mf = dft.UKS(mol, xc=functional).run()
    
    print("\n" + "=" * 50)
    print(f"DFT 结果 ({functional}/{basis})")
    print("=" * 50)
    print(f"总能量: {mf.e_tot:.10f} Hartree")
    print(f"总电子: {mol.nelectron}")
    
    # 轨道能级
    mo_e = mf.mo_energy * 27.2114
    occ = mf.mo_occ
    
    homo_idx = np.where(occ > 0)[0].max()
    lumo_idx = np.where(occ == 0)[0].min()
    
    print(f"\nHOMO: {mo_e[homo_idx]:.4f} eV")
    print(f"LUMO: {mo_e[lumo_idx]:.4f} eV")
    print(f"Gap:  {mo_e[lumo_idx] - mo_e[homo_idx]:.4f} eV")
    
    return {'mf': mf, 'mol': mol, 'symbols': symbols}


def compare_functionals(xyz_file, functionals, basis='cc-pvdz'):
    """
    比较多个泛函的能量和 HOMO-LUMO gap
    
    Args:
        functionals: 泛函名列表
    """
    print("\n泛函比较:")
    print(f"{'Functional':<20} {'E(tot)':>15} {'HOMO':>10} {'LUMO':>10} {'Gap':>10}")
    print("-" * 70)
    
    results = []
    for xc in functionals:
        try:
            result = run_dft(xyz_file, xc, basis, output_file=None)
            mf = result['mf']
            mo_e = mf.mo_energy * 27.2114
            occ = mf.mo_occ
            
            homo = mo_e[occ > 0].max()
            lumo = mo_e[occ == 0].min()
            gap = lumo - homo
            
            print(f"{xc:<20} {mf.e_tot:>15.8f} {homo:>10.4f} {lumo:>10.4f} {gap:>10.4f}")
            results.append({'xc': xc, 'e_tot': mf.e_tot, 'homo': homo, 
                          'lumo': lumo, 'gap': gap})
        except Exception as e:
            print(f"{xc:<20} ERROR: {e}")
    
    return results


def grid_convergence_test(xyz_file, functional='pbe', 
                         grids=[(50, 110), (75, 195), (100, 290)]):
    """
    测试 SCF 收敛性 vs 原子轨道网格大小
    
    grids: list of (atomic_radii, lebedev_n) tuples
    """
    print("\n网格收敛测试:")
    
    for rad, leb in grids:
        try:
            mol = gto.M(atom=open(xyz_file).read(), basis='cc-pvdz')
            mf = dft.RKS(mol, xc=functional)
            mf.grids.atom_radius = rad
            mf.grids.lebedev_grids = {'': (leb, 194)}
            mf.run()
            
            print(f"  Grid ({rad}, {leb}): E = {mf.e_tot:.10f}")
        except Exception as e:
            print(f"  Grid ({rad}, {leb}): FAILED - {e}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\n推荐泛函:")
        for use, func in RECOMMENDED.items():
            print(f"  {use}: {func}")
        print("\n可用泛函类别:")
        for cat in FUNCTIONAL_CATEGORIES:
            funcs = FUNCTIONAL_CATEGORIES[cat]
            print(f"  {cat}: {', '.join(funcs[:5])}...")
        print("\n示例:")
        print("  python dft.py water.csv b3lyp cc-pvdz")
        print("  python dft.py mol.csv wb97x-d3bj def2-tzvp")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    run_dft(xyz_file, functional, basis)


if __name__ == '__main__':
    main()
