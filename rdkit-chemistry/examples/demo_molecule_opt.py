#!/usr/bin/env python3
"""
RDKit 分子优化和可视化演示
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import os

# 选择一个典型的 TADF 分子: DMAC-TRZ (给体: DMAC, 受体: TRZ)
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== 分子优化演示: {NAME} ===\n")

# Step 1: 从 SMILES 构建 3D 分子
print("Step 1: 从 SMILES 构建 3D 构象...")
mol = Chem.MolFromSmiles(SMILES)
if mol is None:
    print("❌ SMILES 解析失败")
    exit(1)

mol = Chem.AddHs(mol)
print(f"✓ 分子构建完成: {mol.GetNumAtoms()} 原子")

# Step 2: 生成 3D 坐标
print("\nStep 2: 生成 3D 坐标...")
res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
if res == -1:
    print("⚠️  默认嵌入失败，尝试随机坐标...")
    res = AllChem.EmbedMolecule(mol, AllChem.ETKDG(), useRandomCoords=True)

if res == -1:
    print("❌ 3D 坐标生成失败")
    exit(1)
print("✓ 3D 坐标生成成功")

# Step 3: MMFF94 力场优化
print("\nStep 3: MMFF94 力场优化...")
try:
    # 构建力场
    mmff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    if mmff is None:
        raise Exception("MMFF 力场初始化失败")

    # 优化
    initial_energy = mmff.CalcEnergy()
    print(f"  初始能量: {initial_energy:.2f} kcal/mol")

    converged = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    if converged == 0:
        final_energy = mmff.CalcEnergy()
        print(f"  最终能量: {final_energy:.2f} kcal/mol")
        print(f"  能量降低: {initial_energy - final_energy:.2f} kcal/mol")
        print("✓ MMFF 优化收敛")
    else:
        print("⚠️  MMFF 优化未完全收敛")

except Exception as e:
    print(f"⚠️  MMFF 失败: {e}")
    print("  尝试 UFF 力场...")
    try:
        AllChem.UFFOptimizeMolecule(mol)
        print("✓ UFF 优化完成")
    except Exception as e2:
        print(f"❌ UFF 也失败: {e2}")
        exit(1)

# Step 4: 导出 XYZ 和 SDF 文件
print("\nStep 4: 导出文件...")
xyz_file = f"{NAME}.xyz"
sdf_file = f"{NAME}.sdf"

# XYZ 文件
Chem.MolToXYZFile(mol, xyz_file)
print(f"✓ 已保存: {xyz_file}")

# SDF 文件（包含键信息）
writer = Chem.SDWriter(sdf_file)
writer.write(mol)
writer.close()
print(f"✓ 已保存: {sdf_file}")

# 读取 XYZ 文件用于后续渲染
with open(xyz_file, 'r') as f:
    xyz_content = f.read()
    print(f"\nXYZ 文件内容预览:")
    print("```")
    lines = xyz_content.split('\n')
    print(lines[0])  # 原子数
    print(lines[1])  # 注释
    for i in range(2, min(7, len(lines))):
        print(lines[i])
    if len(lines) > 7:
        print("...")
    print("```")

# Step 5: 使用 xyzrender 渲染
print("\nStep 5: xyzrender 渲染...")
png_file = f"{NAME}.png"

try:
    # 使用 SDF 文件渲染（保留键信息）
    result = subprocess.run(
        ['xyzrender', sdf_file, '-o', png_file, '--transparent', '--bo'],
        capture_output=True,
        text=True,
        timeout=30
    )

    if result.returncode == 0 and os.path.exists(png_file):
        file_size = os.path.getsize(png_file)
        print(f"✓ 渲染成功: {png_file} ({file_size} bytes)")
    else:
        print(f"❌ 渲染失败:")
        print(result.stderr)
        exit(1)

except FileNotFoundError:
    print("❌ xyzrender 未安装")
    print("  安装: pip install xyzrender")
    exit(1)
except Exception as e:
    print(f"❌ 渲染错误: {e}")
    exit(1)

print("\n" + "="*50)
print("✅ 完成！")
print(f"   分子: {NAME}")
print(f"   原子数: {mol.GetNumAtoms()}")
print(f"   XYZ: {xyz_file}")
print(f"   PNG: {png_file}")
print("="*50)
