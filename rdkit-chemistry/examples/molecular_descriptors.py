#!/usr/bin/env python3
"""
RDKit 分子特征计算和可视化
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
import json

# 使用 DMAC-TRZ 作为示例
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== 分子特征计算: {NAME} ===\n")

# 构建分子
mol = Chem.MolFromSmiles(SMILES)
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)

print("=" * 60)
print("1. 基本分子信息")
print("=" * 60)
print(f"分子式: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")
print(f"分子量: {Descriptors.ExactMolWt(mol):.2f} Da")
print(f"总原子数: {mol_3d.GetNumAtoms()} (含 H)")
print(f"重原子数: {mol.GetNumAtoms()}")
print(f"总键数: {mol.GetNumBonds()}")

print("\n" + "=" * 60)
print("2. 物理化学性质")
print("=" * 60)
print(f"LogP (脂溶性): {Descriptors.MolLogP(mol):.2f}")
print(f"TPSA (拓扑极性表面积): {Descriptors.TPSA(mol):.2f} Å²")
print(f"可旋转键数: {Descriptors.NumRotatableBonds(mol)}")
print(f"氢键供体: {Descriptors.NumHDonors(mol)}")
print(f"氢键受体: {Descriptors.NumHAcceptors(mol)}")

print("\n" + "=" * 60)
print("3. 电子性质")
print("=" * 60)

# Gasteiger 电荷
AllChem.ComputeGasteigerCharges(mol_3d)
charges = [float(atom.GetProp('_GasteigerCharge')) for atom in mol_3d.GetAtoms()]
print(f"原子电荷范围: [{min(charges):.3f}, {max(charges):.3f}]")
print(f"平均电荷: {np.mean(charges):.3f}")

# 识别正负电荷中心
most_positive_idx = int(np.argmax(charges))
most_negative_idx = int(np.argmin(charges))
most_positive_atom = mol_3d.GetAtomWithIdx(most_positive_idx)
most_negative_atom = mol_3d.GetAtomWithIdx(most_negative_idx)

print(f"\n最正电荷原子: {most_positive_atom.GetSymbol()} (idx {most_positive_idx}, charge={charges[most_positive_idx]:.3f})")
print(f"最负电荷原子: {most_negative_atom.GetSymbol()} (idx {most_negative_idx}, charge={charges[most_negative_idx]:.3f})")

print("\n" + "=" * 60)
print("4. 3D 描述符")
print("=" * 60)

# 分子体积和表面积
volume = AllChem.ComputeMolVolume(mol_3d)
print(f"分子体积: {volume:.2f} Å³")

# 计算惯性张量和分子形状
conf = mol_3d.GetConformer()
coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
center = coords.mean(axis=0)

# 计算主轴
coords_centered = coords - center
cov_matrix = np.cov(coords_centered.T)
eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
eigenvalues = np.real(eigenvalues)
eigenvalues = np.sort(eigenvalues)[::-1]

print(f"主轴长度 (σ): {np.sqrt(eigenvalues[0]):.2f}, {np.sqrt(eigenvalues[1]):.2f}, {np.sqrt(eigenvalues[2]):.2f} Å")
print(f"扁平度 (c/a): {np.sqrt(eigenvalues[2]/eigenvalues[0]):.3f}")

print("\n" + "=" * 60)
print("5. 共轭体系分析")
print("=" * 60)

# 芳香环
num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
print(f"芳香环数量: {num_aromatic_rings}")

# 环信息
ring_info = mol.GetRingInfo()
print(f"总环数: {ring_info.NumRings()}")

# π 电子数（近似）
aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
print(f"芳香原子数: {aromatic_atoms}")

print("\n" + "=" * 60)
print("6. 给体-受体特征 (TADF 相关)")
print("=" * 60)

# 识别给体和受体基团
donor_pattern = Chem.MolFromSmarts('[N,n,O,o,S,s]')  # N, O, S
acceptor_pattern = Chem.MolFromSmarts('[N+!$(*[O-]),n+!$(*[O-]),O,o]')  # N+ 或 O

donor_matches = mol.GetSubstructMatches(donor_pattern)
acceptor_matches = mol.GetSubstructMatches(acceptor_pattern)

print(f"潜在给体位点: {len(donor_matches)}")
print(f"潜在受体位点: {len(acceptor_matches)}")

# 特征基团识别
functional_groups = {
    '叔胺': Chem.MolFromSmarts('[N;R0](C)(C)C'),
    '三嗪': Chem.MolFromSmarts('c1ncncn1'),
    '咔唑': Chem.MolFromSmarts('c1ccc2c(c1)[nH]c1ccccc12'),
}

print("\n功能基团:")
for name, pattern in functional_groups.items():
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        print(f"  {name}: {len(matches)} 个")

print("\n" + "=" * 60)
print("7. 生成可视化")
print("=" * 60)

# 7.1 2D 结构图（带原子索引）
print("生成 2D 结构...")
drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
drawer.SetFontSize(0.8)
opts = drawer.drawOptions()
opts.addAtomIndices = True
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
drawer.WriteDrawingText(f"{NAME}_2d.png")
print(f"✓ 保存: {NAME}_2d.png")

# 7.2 2D 结构图（带 Gasteiger 电荷）
print("生成电荷分布图...")
for atom in mol_3d.GetAtoms():
    charge = float(atom.GetProp('_GasteigerCharge'))
    # 归一化到 [-1, 1]
    norm_charge = np.clip(charge / 0.5, -1, 1)
    atom.SetProp('Note', f"{charge:.2f}")

drawer2 = rdMolDraw2D.MolDraw2DCairo(800, 600)
drawer2.DrawMolecule(mol_3d)
drawer2.FinishDrawing()
drawer2.WriteDrawingText(f"{NAME}_charges.png")
print(f"✓ 保存: {NAME}_charges.png")

# 7.3 3D 构象（SDF 格式）
sdf_file = f"{NAME}_3d.sdf"
writer = Chem.SDWriter(sdf_file)
writer.write(mol_3d)
writer.close()
print(f"✓ 保存: {sdf_file}")

# 7.4 导出特征到 JSON
features = {
    "分子信息": {
        "分子式": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "分子量": round(Descriptors.ExactMolWt(mol), 2),
        "重原子数": mol.GetNumAtoms(),
    },
    "物理化学性质": {
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "可旋转键数": Descriptors.NumRotatableBonds(mol),
        "氢键供体": Descriptors.NumHDonors(mol),
        "氢键受体": Descriptors.NumHAcceptors(mol),
    },
    "电子性质": {
        "电荷范围": [round(min(charges), 3), round(max(charges), 3)],
        "最正电荷原子": f"{most_positive_atom.GetSymbol()} (idx {most_positive_idx})",
        "最负电荷原子": f"{most_negative_atom.GetSymbol()} (idx {most_negative_idx})",
    },
    "3D 描述符": {
        "分子体积_Å3": round(volume, 2),
        "主轴长度_Å": [round(np.sqrt(e), 2) for e in eigenvalues],
        "扁平度": round(np.sqrt(eigenvalues[2]/eigenvalues[0]), 3),
    },
    "共轭体系": {
        "芳香环数": num_aromatic_rings,
        "芳香原子数": aromatic_atoms,
    },
    "TADF 特征": {
        "给体位点数": len(donor_matches),
        "受体位点数": len(acceptor_matches),
    }
}

json_file = f"{NAME}_features.json"
with open(json_file, 'w', encoding='utf-8') as f:
    json.dump(features, f, ensure_ascii=False, indent=2)
print(f"✓ 保存: {json_file}")

print("\n" + "=" * 60)
print("✅ 完成！")
print("=" * 60)
print(f"生成的文件:")
print(f"  - {NAME}_2d.png (2D 结构图)")
print(f"  - {NAME}_charges.png (电荷分布)")
print(f"  - {NAME}_3d.sdf (3D 构象)")
print(f"  - {NAME}_features.json (特征数据)")
print("=" * 60)
