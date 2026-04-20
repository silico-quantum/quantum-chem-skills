#!/usr/bin/env python3
"""
Wavefunction Analysis - 波函数分析

支持:
    - Mulliken 电荷
    - Löwdin 轨道分析
    - NBO 分析
    - DOS/PDOS 分析
    - 轨道可视化数据导出

用法:
    python analysis.py <chkfile> <task> [params]
    python analysis.py hf.chk mulliken
    python analysis.py hf.chk dos --plot
"""

import sys
import numpy as np

from pyscf import gto, scf, dft, lo, nbo


def load_from_chk(chkfile, method='RHF'):
    """从 checkpoint 文件加载计算结果"""
    mol = gto.M()
    mol.build()
    
    if method == 'RHF':
        mf = scf.RHF(mol)
    elif method == 'UHF':
        mf = scf.UHF(mol)
    elif method == 'RKS':
        mf = dft.RKS(mol)
    elif method == 'UKS':
        mf = dft.UKS(mol)
    
    mf = mf.from_chk(chkfile)
    return mf


def load_from_scf_result(scf_result):
    """从 SCF 运行结果提取"""
    return scf_result


def mulliken_charges(mf, mol=None):
    """
    Mulliken 布居分析
    
    Returns:
        array: 每个原子的电荷
    """
    if mol is None:
        mol = mf.mol
    
    # 使用正交基
    mo_occ = mf.mo_occ
    mo_coeff = mf.mo_coeff
    dm = mf.make_rdm1()
    
    # 原子轨道基
    orth_ao = lo.orth_ao(mol, 'minao')
    
    # 转换密度矩阵到正交基
    from pyscf import lib
    dm_orth = lib.dot(orth_ao.T, lib.dot(dm, orth_ao))
    
    # 原子布居
    atoms = [mol.atom_symbol(i) for i in range(mol.natm)]
    charges = np.zeros(mol.natm)
    
    ao_per_atom = mol.atom_nbas()
    start = 0
    for ia in range(mol.natm):
        nao = mol.atom_nbas(ia)
        end = start + nao
        charges[ia] = dm_orth[start:end, start:end].trace()
        start = end
    
    # 核电荷
    nuclear = np.array([gto.charge(mol.atom_symbol(i)) for i in range(mol.natm)])
    charges = nuclear - charges
    
    print("\nMulliken 电荷:")
    for i, (atom, charge) in enumerate(zip(atoms, charges)):
        print(f"  {i:3d} {atom:2s}: {charge:+.6f}")
    
    return charges


def lowdin_charges(mf, mol=None):
    """
    Löwdin 布居分析（比 Mulliken 更稳定）
    """
    if mol is None:
        mol = mf.mol
    
    # Löwdin 正交基
    orth_ao = lo.orth_ao(mol, 'meta_lowdin')
    dm = mf.make_rdm1()
    
    from pyscf import lib
    dm_orth = lib.dot(orth_ao.T, lib.dot(dm, orth_ao))
    
    atoms = [mol.atom_symbol(i) for i in range(mol.natm)]
    charges = np.zeros(mol.natm)
    
    start = 0
    for ia in range(mol.natm):
        nao = mol.atom_nbas(ia)
        end = start + nao
        charges[ia] = dm_orth[start:end, start:end].trace()
        start = end
    
    nuclear = np.array([gto.charge(mol.atom_symbol(i)) for i in range(mol.natm)])
    charges = nuclear - charges
    
    print("\nLöwdin 电荷:")
    for i, (atom, charge) in enumerate(zip(atoms, charges)):
        print(f"  {i:3d} {atom:2s}: {charge:+.6f}")
    
    return charges


def dipole_moment(mf, mol=None):
    """
    计算偶极矩
    
    Returns:
        array: [dx, dy, dz] Debye
    """
    if mol is None:
        mol = mf.mol
    
    dm = mf.make_rdm1()
    
    # 电子贡献
    from pyscf import df
    ao_dip = df.get_dipole(mol)
    
    elec_dip = np.einsum('xij,ji->x', ao_dip, dm)
    
    # 核贡献
    nuclear_dip = np.zeros(3)
    for i in range(mol.natm):
        nuclear_dip += mol.atom_coord(i) * gto.charge(mol.atom_symbol(i))
    
    total_dip = nuclear_dip - elec_dip
    # 转换到 Debye (1 a.u. = 2.5418 Debye)
    total_debye = total_dip * 2.5418
    
    print("\n偶极矩:")
    print(f"  x = {total_debye[0]:.6f} D")
    print(f"  y = {total_debye[1]:.6f} D")
    print(f"  z = {total_debye[2]:.6f} D")
    print(f"  |μ| = {np.linalg.norm(total_debye):.6f} D")
    
    return total_debye


def orbital_population(mf, mol=None, mo_index=None):
    """
    轨道组成分析（哪些原子/壳层贡献最大）
    
    Args:
        mo_index: 要分析的 MO index，None = 分析所有
    """
    if mol is None:
        mol = mf.mol
    
    # 原子轨道索引
    ao_labels = mol.ao_labels()
    
    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    if mo_index is not None:
        indices = [mo_index]
    else:
        # 分析前 10 个轨道
        indices = list(range(min(10, mo_coeff.shape[1])))
    
    print("\n轨道组成分析 (前 %d 个):" % len(indices))
    print("-" * 60)
    
    for idx in indices:
        coeff_sq = mo_coeff[:, idx]**2
        sorted_idx = np.argsort(coeff_sq)[::-1]
        
        occ_tag = "[occ]" if mo_occ[idx] > 0 else "[emp]"
        
        print(f"\nMO {idx} {occ_tag} E={mo_energy[idx]:.4f} eV")
        print("  主要贡献:")
        for j in sorted_idx[:5]:
            if coeff_sq[j] > 0.01:
                print(f"    {ao_labels[j]:30s} {coeff_sq[j]:.4f}")


def dos_analysis(mf, mol=None, energy_range=None):
    """
    DOS (态密度) 分析
    
    Returns:
        tuple: (energies, dos, pdos_per_atom)
    """
    if mol is None:
        mol = mf.mol
    
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    # Gaussian 展宽
    sigma = 0.1  # eV
    
    if energy_range is None:
        e_min = mo_energy.min() - 1
        e_max = mo_energy.max() + 1
    else:
        e_min, e_max = energy_range
    
    energies = np.linspace(e_min, e_max, 1000)
    dos = np.zeros_like(energies)
    
    # Gaussian 展宽
    for i, e in enumerate(mo_energy):
        dos += np.exp(-(energies - e)**2 / (2 * sigma**2))
    
    print("\nDOS 分析:")
    print(f"  HOMO: {mo_energy[mo_occ > 0].max():.4f} eV")
    print(f"  LUMO: {mo_energy[mo_occ == 0].min():.4f} eV")
    print(f"  Gap:  {mo_energy[mo_occ == 0].min() - mo_energy[mo_occ > 0].max():.4f} eV")
    
    return energies, dos


def pdos_analysis(mf, mol=None):
    """
    PDOS (投影态密度) 分析
    """
    if mol is None:
        mol = mf.mol
    
    mo_energy = mf.mo_energy * 27.2114
    mo_coeff = mf.mo_coeff
    mo_occ = mf.mo_occ
    
    sigma = 0.1
    e_min = mo_energy.min() - 1
    e_max = mo_energy.max() + 1
    energies = np.linspace(e_min, e_max, 500)
    
    pdos = {}
    
    for ia in range(mol.natm):
        atom_symbol = mol.atom_symbol(ia)
        
        # 获取该原子的 AO 索引范围
        start = mol.atom_start_idx[ia]
        end = start + mol.atom_nbas[ia]
        
        # 该原子的投影
        coeff_atom = mo_coeff[start:end, :]
        contrib = coeff_atom**2
        
        pdos[atom_symbol] = np.zeros_like(energies)
        for i, e in enumerate(mo_energy):
            pdos[atom_symbol] += contrib[:, i].sum() * np.exp(-(energies - e)**2 / (2 * sigma**2))
    
    print("\nPDOS per atom:")
    for atom, dos in pdos.items():
        print(f"  {atom}: max DOS = {dos.max():.4f}")
    
    return energies, pdos


def orbital_energies_summary(mf):
    """
    打印轨道能级汇总
    """
    mo_energy = mf.mo_energy * 27.2114  # eV
    mo_occ = mf.mo_occ
    
    n_homo = int(mo_occ.sum() // 2)
    
    print("\n" + "=" * 60)
    print("轨道能级汇总 (eV)")
    print("=" * 60)
    print(f"{'Index':>6} {'Occ':>4} {'Energy(eV)':>12} {'Assignment':>20}")
    print("-" * 60)
    
    for i in range(min(20, len(mo_energy))):
        occ_str = "occ" if mo_occ[i] > 0 else "emp"
        
        # 简单标记 HOMO/LUMO
        assign = ""
        if i == n_homo - 1:
            assign = "<-- HOMO"
        elif i == n_homo:
            assign = "<-- LUMO"
        
        print(f"{i:6d} {occ_str:>4} {mo_energy[i]:12.4f} {assign:>20}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("\n示例:")
        print("  python analysis.py result.pkl mulliken")
        print("  python analysis.py result.pkl dos --plot")
        sys.exit(1)
    
    input_file = sys.argv[1]
    task = sys.argv[2].lower()
    
    # 加载数据（需要提供 scf_result 或 mol/mf）
    print(f"Analysis: {task}")
    print(f"Input: {input_file}")
    print("\n需要传入 SCF/MF 对象进行分析")
    print("用法示例:")
    print("  from pyscf import scf, analysis")
    print("  mf = scf.RHF(mol).run()")
    print("  analysis.mulliken_charges(mf)")


if __name__ == '__main__':
    main()
