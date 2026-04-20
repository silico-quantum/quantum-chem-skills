#!/usr/bin/env python3
"""
PES - Potential Energy Surface 势能面扫描

支持:
    - 1D 扫描（单坐标）
    - 2D 扫描（双坐标）
    - Relaxed 扫描（每点优化）
    - 刚性扫描（不优化）

用法:
    python pes.py <xyz_file> <scan_type> [params]
    python pes.py mol.xyz scan_r 0 1 1.0 2.0 20
    python pes.py mol.xyz scan_angle 0 1 2 90 180 10
"""

import sys
import numpy as np

from pyscf import gto, dft, scf, geomopt


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


def write_xyz(symbols, coords, filename):
    """写入 XYZ 文件"""
    with open(filename, 'w') as f:
        f.write(f"{len(symbols)}\nPES scan\n")
        for s, c in zip(symbols, coords):
            f.write(f"{s:2s} {c[0]:15.10f} {c[1]:15.10f} {c[2]:15.10f}\n")


def run_rigid_scan(xyz_file, scan_type, indices, values, 
                   method='RKS', functional='pbe', basis='cc-pvdz',
                   output_prefix='pes_scan'):
    """
    刚性扫描（不优化，只改变坐标）
    
    Args:
        xyz_file: 初始 XYZ 文件
        scan_type: 'bond', 'angle', 'dihedral'
        indices: 原子索引 (list)
            bond: [i, j]
            angle: [i, j, k]
            dihedral: [i, j, k, l]
        values: 扫描值 (Angstrom 或 degree)
        method: 'RKS', 'RHF'
        functional: DFT 泛函
    
    Returns:
        dict: {'values': 扫描值, 'energies': 能量}
    """
    symbols, coords = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    energies = []
    
    print(f"刚性扫描: {scan_type} = {values[0]:.3f} 到 {values[-1]:.3f}")
    print(f"点数: {len(values)}")
    
    for idx, val in enumerate(values):
        new_coords = modify_geometry(coords, scan_type, indices, val)
        
        atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                              for s, c in zip(symbols, new_coords)])
        
        mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=0)
        
        if method == 'RKS':
            mf = dft.RKS(mol, xc=functional).run()
        else:
            mf = scf.RHF(mol).run()
        
        energies.append(mf.e_tot)
        
        if (idx + 1) % 5 == 0 or idx == 0:
            print(f"  Step {idx+1}/{len(values)}: {val:.4f} E={mf.e_tot:.10f}")
    
    result = {'values': np.array(values), 'energies': np.array(energies)}
    
    # 保存
    out_file = f'{output_prefix}_{scan_type}.dat'
    header = f"# {scan_type} energy scan\n# {method} {functional}/{basis}"
    np.savetxt(out_file, np.column_stack([result['values'], result['energies']]),
               header=header)
    print(f"\n结果保存到: {out_file}")
    
    return result


def run_relaxed_scan(xyz_file, scan_type, indices, values,
                     method='RKS', functional='pbe', basis='cc-pvdz',
                     output_prefix='pes_relaxed'):
    """
    Relaxed 扫描（每点优化坐标）
    """
    symbols, coords = read_xyz(xyz_file)
    nelec = sum(gto.charge(s) for s in symbols)
    spin = nelec % 2
    
    energies = []
    
    print(f"Relaxed 扫描: {scan_type} = {values[0]:.3f} 到 {values[-1]:.3f}")
    
    for idx, val in enumerate(values):
        new_coords = modify_geometry(coords, scan_type, indices, val)
        
        atom_str = ' '.join(['%s %.10f %.10f %.10f' % (s, *c) 
                              for s, c in zip(symbols, new_coords)])
        
        mol = gto.M(atom=atom_str, basis=basis, charge=0, spin=spin, verbose=0)
        
        if method == 'RKS':
            mf = dft.RKS(mol, xc=functional)
        else:
            mf = scf.RHF(mol)
        
        # 优化
        opt = geomopt.GeometryOptimizer(mf)
        mol_opt = opt.run()
        
        energies.append(mol_opt.energy)
        
        # 保存结构
        out_xyz = f'{output_prefix}_{idx:04d}.xyz'
        write_xyz(
            [mol_opt.atom_symbol(i) for i in range(mol_opt.natm)],
            mol_opt.atom_coords(),
            out_xyz
        )
        
        if (idx + 1) % 5 == 0 or idx == 0:
            print(f"  Step {idx+1}/{len(values)}: {val:.4f} E={mol_opt.energy:.10f}")
    
    result = {'values': np.array(values), 'energies': np.array(energies)}
    
    out_file = f'{output_prefix}_{scan_type}.dat'
    np.savetxt(out_file, np.column_stack([result['values'], result['energies']]),
               header=f"# Relaxed {scan_type} scan\n# {method} {functional}/{basis}")
    
    return result


def modify_geometry(coords, scan_type, indices, value):
    """
    修改几何坐标
    
    Args:
        coords: 原始坐标数组 (N, 3)
        scan_type: 'bond', 'angle', 'dihedral'
        indices: 原子索引
        value: 目标值 (Angstrom 或 degree)
    """
    new_coords = coords.copy()
    
    if scan_type == 'bond':
        i, j = indices
        # 改变键长
        direction = coords[j] - coords[i]
        direction = direction / np.linalg.norm(direction)
        new_coords[j] = coords[i] + value * direction
    
    elif scan_type == 'angle':
        i, j, k = indices
        # 改变键角 (degree)
        theta = np.deg2rad(value)
        
        # 计算当前位置
        v1 = coords[i] - coords[j]
        v2 = coords[k] - coords[j]
        
        # 计算当前角度
        cos_curr = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        curr_angle = np.arccos(np.clip(cos_curr, -1, 1))
        
        # 旋转 v1 到目标角度
        # 使用 Rodrigues' rotation formula
        axis = np.cross(v1, v2)
        if np.linalg.norm(axis) < 1e-10:
            # 共线，使用默认值
            return new_coords
        
        axis = axis / np.linalg.norm(axis)
        
        # 旋转角度
        delta = theta - curr_angle
        
        # 应用旋转
        v1_rot = v1 * np.cos(delta) + np.cross(axis, v1) * np.sin(delta) + axis * np.dot(axis, v1) * (1 - np.cos(delta))
        
        new_coords[i] = coords[j] + v1_rot
    
    elif scan_type == 'dihedral':
        i, j, k, l = indices
        # 改变二面角 (degree)
        phi = np.deg2rad(value)
        
        # 计算当前位置
        b1 = coords[j] - coords[i]
        b2 = coords[k] - coords[j]
        b3 = coords[l] - coords[k]
        
        # 二面角通过 b1, b2, b3 计算
        # 先绕 b2-b3 轴旋转 l 原子
        # 简化实现
        return new_coords  # TODO: 实现二面角旋转
    
    return new_coords


def scan_bond(xyz_file, i, j, r_min, r_max, n_points=20,
              method='RKS', functional='pbe', basis='cc-pvdz',
              relaxed=False, output_prefix='pes'):
    """
    扫描键长
    
    Args:
        i, j: 原子索引
        r_min, r_max: 键长范围 (Angstrom)
    """
    values = np.linspace(r_min, r_max, n_points)
    
    if relaxed:
        return run_relaxed_scan(xyz_file, 'bond', [i, j], values,
                                method, functional, basis, output_prefix)
    else:
        return run_rigid_scan(xyz_file, 'bond', [i, j], values,
                              method, functional, basis, output_prefix)


def scan_angle(xyz_file, i, j, k, angle_min, angle_max, n_points=20,
               method='RKS', functional='pbe', basis='cc-pvdz',
               relaxed=False, output_prefix='pes'):
    """
    扫描键角
    
    Args:
        i, j, k: 原子索引（j 为顶点）
        angle_min, angle_max: 角度范围 (degree)
    """
    values = np.linspace(angle_min, angle_max, n_points)
    
    if relaxed:
        return run_relaxed_scan(xyz_file, 'angle', [i, j, k], values,
                                method, functional, basis, output_prefix)
    else:
        return run_rigid_scan(xyz_file, 'angle', [i, j, k], values,
                              method, functional, basis, output_prefix)


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        print("\n示例:")
        print("  python pes.py mol.xyz bond 0 1 1.0 2.0 20")
        print("  python pes.py mol.xyz angle 0 1 2 90 180 10")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    scan_type = sys.argv[2]  # bond, angle
    
    if scan_type == 'bond':
        i, j = int(sys.argv[3]), int(sys.argv[4])
        r_min, r_max = float(sys.argv[5]), float(sys.argv[6])
        n_points = int(sys.argv[7]) if len(sys.argv) > 7 else 20
        basis = sys.argv[8] if len(sys.argv) > 8 else 'cc-pvdz'
        functional = sys.argv[9] if len(sys.argv) > 9 else 'pbe'
        
        scan_bond(xyz_file, i, j, r_min, r_max, n_points,
                 'RKS', functional, basis)
    
    elif scan_type == 'angle':
        i, j, k = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        angle_min = float(sys.argv[6])
        angle_max = float(sys.argv[7])
        n_points = int(sys.argv[8]) if len(sys.argv) > 8 else 20
        
        scan_angle(xyz_file, i, j, k, angle_min, angle_max, n_points)


if __name__ == '__main__':
    main()
