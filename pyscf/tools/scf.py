#!/usr/bin/env python3
"""
SCF - Self-Consistent Field calculations
Hartree-Fock (RHF, UHF, ROHF) 基础框架

用法:
    python scf.py <xyz_file> <method> [basis]
    python scf.py water.xyz RHF cc-pvdz
    python scf.py molecule.xyz UHF sto-3g
"""

import sys
import numpy as np

from pyscf import gto, scf


def read_xyz(xyz_file):
    """从 XYZ 文件读取分子结构"""
    with open(xyz_file) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    coords = []
    symbols = []
    for i in range(2, 2 + natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords)


def build_mol(symbols, coords, basis='cc-pvdz', charge=0, spin=0, unit='Angstrom'):
    """构建分子对象"""
    atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) for s, c in zip(symbols, coords)])
    mol = gto.M(
        atom=atom_str,
        basis=basis,
        charge=charge,
        spin=spin,
        unit=unit,
        verbose=4
    )
    return mol


def run_rhf(mol, diis=True, direct_scf=False, chkfile=None):
    """
    闭壳层 RHF 计算
    
    Args:
        mol: PySCF 分子对象
        diis: 是否使用 DIIS 加速收敛
        direct_scf: 是否使用 direct SCF（省内存）
        chkfile: checkpoint 文件路径
    """
    mf = scf.RHF(mol)
    
    if diis:
        mf = scf.RHF(mol).diis()
    if direct_scf:
        mf.direct_scf = True
    if chkfile:
        mf.chkfile = chkfile
    
    return mf.run()


def run_uhf(mol, diis=True, chkfile=None):
    """
    开壳层 UHF 计算
    
    Args:
        mol: PySCF 分子对象（spin > 0）
        diis: 是否使用 DIIS 加速
    """
    mf = scf.UHF(mol)
    if diis:
        mf = mf.diis()
    return mf.run()


def run_rohf(mol, diis=True):
    """
    限制性开壳层 ROHF
    """
    mf = scf.ROHF(mol)
    if diis:
        mf = mf.diis()
    return mf.run()


def run_rhf_with_guess(mol, guess='atom', method='RHF', chkfile=None):
    """
    带指定初始猜测的 SCF 计算
    
    Args:
        mol: PySCF 分子对象
        guess: 'atom', 'hcore', '1e', 'chk'
        method: 'RHF', 'UHF', 'ROHF', 'RKS', 'UKS'
    """
    if method == 'RHF':
        mf_class = scf.RHF
    elif method == 'UHF':
        mf_class = scf.UHF
    elif method == 'ROHF':
        mf_class = scf.ROHF
    elif method == 'RKS':
        mf_class = scf.RKS
    elif method == 'UKS':
        mf_class = scf.UKS
    else:
        raise ValueError(f"Unknown method: {method}")
    
    mf = mf_class(mol)
    mf.init_guess = guess
    if chkfile:
        mf.chkfile = chkfile
    return mf.run()


def newton_rhf(mol, xc=None):
    """
    Newton-Raphson 加速的 SCF（适合难以收敛的体系）
    
    Args:
        mol: PySCF 分子对象
        xc: 泛函名（如 'b3lyp'），None 则为 HF
    """
    if xc:
        from pyscf import dft
        if mol.spin == 0:
            mf_class = dft.RKS
        else:
            mf_class = dft.UKS
        mf = mf_class(mol, xc=xc)
    else:
        if mol.spin == 0:
            mf_class = scf.RHF
        else:
            mf_class = scf.UHF
        mf = mf_class(mol)
    
    return mf.newton().run()


def density_fit_scf(mol, method='RHF', auxbasis=None):
    """
    密度拟合加速的 SCF（省内存，适合大体系）
    
    Args:
        mol: PySCF 分子对象
        method: 'RHF', 'UHF', 'RKS', 'UKS'
        auxbasis: 辅助基组，默认自动选择
    """
    if method == 'RHF':
        mf = scf.RHF(mol).density_fit(auxbasis=auxbasis).run()
    elif method == 'UHF':
        mf = scf.UHF(mol).density_fit(auxbasis=auxbasis).run()
    elif method == 'RKS':
        from pyscf import dft
        mf = dft.RKS(mol).density_fit(auxbasis=auxbasis).run()
    elif method == 'UKS':
        from pyscf import dft
        mf = dft.UKS(mol).density_fit(auxbasis=auxbasis).run()
    return mf


def get_orbital_info(mf):
    """
    打印轨道信息
    
    Returns:
        dict: HOMO, LUMO, HOMO-LUMO gap (eV)
    """
    occ = mf.mo_occ
    energy = mf.mo_energy * 27.2114  # Hartree -> eV
    
    homo_idx = np.where(occ > 0)[0].max()
    lumo_idx = np.where(occ == 0)[0].min()
    
    homo = energy[homo_idx]
    lumo = energy[lumo_idx]
    gap = lumo - homo
    
    print("=" * 50)
    print("轨道信息 (eV)")
    print("=" * 50)
    print(f"HOMO:  {homo:.4f} eV")
    print(f"LUMO:  {lumo:.4f} eV")
    print(f"Gap:   {gap:.4f} eV")
    print(f"Gap:   {gap*8065.5:.2f} cm⁻¹")
    
    return {'homo': homo, 'lumo': lumo, 'gap': gap, 'homo_idx': homo_idx, 'lumo_idx': lumo_idx}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\n示例:")
        print("  python scf.py water.xyz RHF cc-pvdz")
        print("  python scf.py o2.xyz UHF sto-3g")
        print("  python scf.py molecule.xyz RKS b3lyp")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    method = sys.argv[2].upper()  # RHF, UHF, RKS, UKS
    basis = sys.argv[3] if len(sys.argv) > 3 else 'cc-pvdz'
    
    # 读取分子
    symbols, coords = read_xyz(xyz_file)
    
    # 确定自旋
    nelec_total = sum(gto.charge(s) for s in symbols)
    nelectron = nelec_total - 0  # charge = 0
    spin = nelectron % 2  # 0 = closed shell, 1 = open shell
    
    print(f"分子: {' '.join(symbols)}")
    print(f"基组: {basis}")
    print(f"电荷: 0, 自旋: {spin}")
    
    # 构建分子
    mol = build_mol(symbols, coords, basis=basis, charge=0, spin=spin)
    
    # 运行 SCF
    if method == 'RHF':
        mf = run_rhf(mol)
    elif method == 'UHF':
        mf = run_uhf(mol)
    elif method == 'RKS':
        xc = sys.argv[4] if len(sys.argv) > 4 else 'pbe'
        from pyscf import dft
        mf = dft.RKS(mol, xc=xc).run()
    elif method == 'UKS':
        xc = sys.argv[4] if len(sys.argv) > 4 else 'pbe'
        from pyscf import dft
        mf = dft.UKS(mol, xc=xc).run()
    elif method == 'ROHF':
        mf = run_rohf(mol)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    print(f"\n总能量: {mf.e_tot:.10f} Hartree")
    
    # 轨道信息
    get_orbital_info(mf)


if __name__ == '__main__':
    main()
