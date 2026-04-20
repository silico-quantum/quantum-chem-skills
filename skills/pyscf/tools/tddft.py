#!/usr/bin/env python3
"""
TDDFT / TDA - 激发态计算

支持:
    - TDDFT (线性响应)
    - TDA (Tamm-Dancoff 近似)
    - CIS (Configuration Interaction Singles)
    - SF-DFT (态选择 DFT)

用法:
    python tddft.py <xyz_file> <functional> [nstates] [basis]
    python tddft.py mol.xyz b3lyp 10 cc-pvdz
    python tddft.py mol.xyz pbe 20 6-31g*
"""

import sys
import numpy as np

from pyscf import gto, dft, tddft, scf


def read_xyz(xyz_file):
    """从 XYZ 文件读取分子"""
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols = []
    coords = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def run_tddft(xyz_file, functional='pbe', nstates=10, basis='cc-pvdz', 
              restricted=True, nstates_tda=None):
    """
    运行 TDDFT 计算
    
    Args:
        xyz_file: XYZ 坐标文件
        functional: DFT 泛函名 ('pbe', 'b3lyp', 'wb97x', etc.)
        nstates: 激发态数量
        basis: 基组名
        restricted: True = RKS, False = UKS
        nstates_tda: TDA 单独指定（用于加速）
    
    Returns:
        dict: {
            'td': TDDFT 对象,
            'mf': SCF 对象,
            'mol': 分子对象,
            'e_ev': 激发能 (eV),
            'f': 振子强度,
            'mo_E': 轨道能量 (eV)
        }
    """
    # 读取分子
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    coords = np.array(coords)
    
    # 确定电荷和自旋
    nelec_total = sum(gto.charge(s) for s in symbols)
    spin = nelec_total % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=4)
    
    # DFT 计算
    if restricted:
        mf = dft.RKS(mol, xc=functional).run()
    else:
        mf = dft.UKS(mol, xc=functional).run()
    
    print(f"DFT 基态计算完成: E = {mf.e_tot:.10f} Hartree")
    
    # TDDFT 计算
    td = tddft.TDDFT(mf).run(nstates=nstates)
    
    # TDA 计算（可选，更快但对 double excitation 不准）
    if nstates_tda:
        td_tda = tddft.TDA(mf).run(nstates=nstates_tda)
    
    # 提取结果
    e_ev = td.e * 27.2114  # Hartree -> eV
    f = td.f
    
    print("\n" + "=" * 60)
    print(f"TDDFT 激发态 (functional={functional})")
    print("=" * 60)
    print(f"{'态':>4} {'能量(eV)':>12} {'波长(nm)':>12} {'振子强度':>12}")
    print("-" * 60)
    for i in range(min(10, len(e_ev))):
        wavelength = 1240.0 / e_ev[i] if e_ev[i] > 0 else float('inf')
        print(f"{i+1:>4} {e_ev[i]:>12.4f} {wavelength:>12.2f} {f[i]:>12.6f}")
    
    return {
        'td': td,
        'mf': mf,
        'mol': mol,
        'e_ev': e_ev,
        'f': f,
        'symbols': symbols,
        'coords': coords
    }


def run_tda(xyz_file, functional='pbe', nstates=20, basis='cc-pvdz'):
    """
    运行 TDA (Tamm-Dancoff) 计算
    比 TDDFT 快约 2 倍，适合快速筛选
    """
    result = run_tddft(xyz_file, functional, nstates, basis)
    td_tda = tddft.TDA(result['mf']).run(nstates=nstates)
    result['td_tda'] = td_tda
    result['e_ev_tda'] = td_tda.e * 27.2114
    result['f_tda'] = td_tda.f
    return result


def run_cis(xyz_file, basis='cc-pvdz', nstates=10):
    """
    CIS (Configuration Interaction Singles)
    基于 HF 而非 DFT 的激发态计算
    """
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    symbols, coords = [], []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[j]) for j in range(1, 4)])
    coords = np.array(coords)
    
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                          for s, c in zip(symbols, coords)])
    mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin)
    
    mf = scf.RHF(mol).run()
    
    from pyscf import cis
    mycis = cis.CIS(mf).run(nstates=nstates)
    
    print("\nCIS 激发态:")
    for i in range(min(10, len(mycis.e))):
        print(f"  态 {i+1}: {mycis.e[i]*27.2114:.4f} eV, f={mycis.f[i]:.6f}")
    
    return {'cis': mycis, 'mf': mf}


def analyze_excited_states(result, top_n=5):
    """
    分析激发态：打印主要跃迁成分
    """
    td = result['td']
    print("\n" + "=" * 60)
    print("主要跃迁分析 (TOP %d 激发态)" % top_n)
    print("=" * 60)
    
    td.analyze()
    
    return td.osc_filter


def filter_by_emission(result, emin=400, emax=800):
    """
    筛选指定波长范围的激发态
    
    Args:
        result: run_tddft() 返回的结果
        emin, emax: 波长范围 (nm)
    
    Returns:
        list: 符合条件的激发态 index
    """
    e_ev = result['e_ev']
    f = result['f']
    
    wavelengths = 1240.0 / e_ev
    valid = []
    
    print(f"\n波长筛选 ({emin}-{emax} nm):")
    print(f"{'态':>4} {'波长(nm)':>12} {'振子强度':>12}")
    print("-" * 40)
    
    for i in range(len(e_ev)):
        wl = wavelengths[i]
        if emin <= wl <= emax and f[i] > 0:
            valid.append(i)
            print(f"{i+1:>4} {wl:>12.2f} {f[i]:>12.6f}")
    
    print(f"\n共 {len(valid)} 个态符合条件")
    return valid


def generate_spectrum(result, nstates=None, sigma=0.05, emin=0, emax=15):
    """
    生成 UV-Vis 吸收光谱
    
    Args:
        result: run_tddft() 返回结果
        nstates: 使用的激发态数量（None = 全部）
        sigma: 高斯展宽 (eV)
        emin, emax: 能量范围 (eV)
    
    Returns:
        tuple: (energies, spectrum) 数组
    """
    if nstates is None:
        nstates = len(result['e_ev'])
    
    e = result['e_ev'][:nstates]
    f = result['f'][:nstates]
    
    energies = np.linspace(emin, emax, 2000)
    spectrum = np.zeros_like(energies)
    
    for ei, fi in zip(e, f):
        if fi > 0:
            spectrum += fi * np.exp(-(energies - ei)**2 / (2 * sigma**2))
    
    return energies, spectrum


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\n示例:")
        print("  python tddft.py mol.xyz b3lyp 20 cc-pvdz")
        print("  python tddft.py mol.xyz pbe 10 6-31g*")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    functional = sys.argv[2]
    nstates = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    basis = sys.argv[4] if len(sys.argv) > 4 else 'cc-pvdz'
    
    print(f"计算: {xyz_file}")
    print(f"泛函: {functional}, 基组: {basis}, nstates={nstates}")
    
    result = run_tddft(xyz_file, functional, nstates, basis)
    
    # 分析
    analyze_excited_states(result)
    
    # 光谱筛选
    valid = filter_by_emission(result, emin=400, emax=800)


if __name__ == '__main__':
    main()
