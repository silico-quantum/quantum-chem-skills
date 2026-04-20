"""
PySCF DFT计算完整示例
演示水分子DFT计算、激发态分析和性质计算
"""

from pyscf import gto, dft, tdscf, lib
import numpy as np

def water_dft():
    """水分子DFT计算"""

    print("="*60)
    print("1. 分子定义")
    print("="*60)

    # 定义水分子
    mol = gto.M(
        atom='''
        O  0.000000  0.000000  0.117790
        H  0.000000  0.755453 -0.471161
        H  0.000000 -0.755453 -0.471161
        ''',
        basis='cc-pvdz',
        charge=0,
        spin=0,
        unit='Ang',
        symmetry='c2v'
    )

    print(f"分子: {mol.atom[0][0]}{'H'*2}")
    print(f"电子数: {mol.nelectron}")
    print(f"AO数: {mol.nao}")
    print(f"对称性: {mol.symmetry}")
    print()

    print("="*60)
    print("2. DFT计算 (B3LYP)")
    print("="*60)

    # RKS-DFT
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.grids.atom_grid = (75, 302)

    # 运行SCF
    e_tot = mf.kernel()

    print(f"\n基态能量: {e_tot:.10f} Hartree")
    print(f"能量: {e_tot*27.2114:.6f} eV")
    print(f"收敛: {mf.converged}")
    print()

    # 分子轨道分析
    print("="*60)
    print("3. 分子轨道分析")
    print("="*60)

    # 轨道能量和占据数
    mo_energy = mf.mo_energy
    mo_occ = mf.mo_occ

    # HOMO
    homo_idx = np.where(mo_occ > 0)[0][-1]
    homo_energy = mo_energy[homo_idx]

    # LUMO
    lumo_idx = np.where(mo_occ == 0)[0][0]
    lumo_energy = mo_energy[lumo_idx]

    # 能隙
    gap = lumo_energy - homo_energy

    print(f"HOMO能量: {homo_energy:.6f} Hartree ({homo_energy*27.2114:.3f} eV)")
    print(f"LUMO能量: {lumo_energy:.6f} Hartree ({lumo_energy*27.2114:.3f} eV)")
    print(f"能隙: {gap:.6f} Hartree ({gap*27.2114:.3f} eV)")
    print()

    # 占据轨道
    print("占据轨道:")
    for i in range(mol.nao):
        if mo_occ[i] > 0:
            print(f"  MO {i+1:2d}: E = {mo_energy[i]:10.6f} Ha,  occ = {mo_occ[i]:.2f}")
    print()

    # 虚轨道
    print("虚轨道 (前5个):")
    for i in range(lumo_idx, min(lumo_idx+5, mol.nao)):
        print(f"  MO {i+1:2d}: E = {mo_energy[i]:10.6f} Ha,  occ = {mo_occ[i]:.2f}")
    print()

    print("="*60)
    print("4. LR-TDDFT激发态")
    print("="*60)

    # TDDFT计算
    td = tdscf.TDDFT(mf)
    td.nstates = 6
    td.kernel()

    # 分析激发态
    print("\n激发态分析:")
    print(f"{'态':<4} {'能量(eV)':<10} {'波长(nm)':<12} {'振子强度':<12}")
    print("-"*60)

    for i in range(td.nstates):
        e_ev = td.e[i] * 27.2114
        wavelength = 1240.0 / e_ev
        f = td.oscillator_strength[i]
        print(f"{i+1:<4} {e_ev:<10.3f} {wavelength:<12.2f} {f:<12.3f}")
    print()

    # 最强激发态
    strongest = np.argmax(td.oscillator_strength)
    print(f"最强激发态: 第{strongest+1}激发态, 振子强度 = {td.oscillator_strength[strongest]:.3f}")
    print()

    print("="*60)
    print("5. NTO分析 (第1激发态)")
    print("="*60)

    # 自然跃迁轨道
    weights, nto = td.get_nto(state=0)
    hole = nto[0]  # 空穴轨道
    electron = nto[1]  # 电子轨道

    print(f"主导跃迁权重: {weights.max():.6f}")
    print(f"空穴轨道数: {hole.shape[1]}")
    print(f"电子轨道数: {electron.shape[1]}")
    print()

    print("="*60)
    print("6. 偶极矩")
    print("="*60)

    # 电偶极矩
    dip = mf.dip_moment(unit='Debye')
    print(f"电偶极矩 (Debye):")
    print(f"  μ_x = {dip[0]:.6f}")
    print(f"  μ_y = {dip[1]:.6f}")
    print(f"  μ_z = {dip[2]:.6f}")
    print(f"  |μ| = {np.linalg.norm(dip):.6f} Debye")
    print()

    print("="*60)
    print("7. 电荷布居分析")
    print("="*60)

    # Mulliken布居
    pop = mf.mulliken_pop(mol, mf.make_rdm1())

    print("Mulliken布居:")
    for i in range(mol.natm):
        atom = mol.atom_pure_symbol(i)
        charge = pop[0][i]
        print(f"  {atom}: {charge:.6f} e")
    print()

    print("="*60)
    print("8. 不同泛函比较")
    print("="*60)

    functionals = ['pbe', 'b3lyp', 'wb97x-d', 'cam-b3lyp']
    results = {}

    for xc in functionals:
        mf_xc = dft.RKS(mol)
        mf_xc.xc = xc
        e_xc = mf_xc.kernel()
        results[xc] = e_xc
        print(f"{xc:12s}: {e_xc:.10f} Hartree")

    # 能量差
    print("\n相对能量 (相对于PBE):")
    e_ref = results['pbe']
    for xc, e in results.items():
        delta = (e - e_ref) * 27.2114 * 627.509  # 转换为kcal/mol
        print(f"{xc:12s}: {delta:10.3f} kcal/mol")
    print()

    print("="*60)
    print("9. 保存结果")
    print("="*60)

    # 保存到检查点
    mf.chkfile = 'water_dft.chk'
    mf.dump_chk()

    print(f"检查点文件: {mf.chkfile}")
    print(f"包含: MO系数、轨道能量、密度矩阵等")
    print()

    return mf, td

if __name__ == '__main__':
    # 设置输出级别
    lib.logger.TIMER_LEVEL = 3
    lib.logger.INFO = True

    # 运行计算
    mf, td = water_dft()

    print("="*60)
    print("计算完成!")
    print("="*60)
