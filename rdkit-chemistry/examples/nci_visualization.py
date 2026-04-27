#!/usr/bin/env python3
"""
电荷和非共价相互作用可视化（简化版，仅使用 RDKit）
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np

# 使用 DMAC-TRZ
SMILES = "Cc1c2c(cc3ccccc13)N(c1ccccc1)C2C1=CC=C2C=CN=C(c3ccccc3)N=C(c3ccccc3)N=C21"
NAME = "DMAC-TRZ"

print(f"=== 电荷和非共价相互作用分析: {NAME} ===\n")

# 构建分子
mol = Chem.MolFromSmiles(SMILES)
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
AllChem.MMFFOptimizeMolecule(mol_3d)

# 计算电荷
AllChem.ComputeGasteigerCharges(mol_3d)

def getCharge(atom):
    """获取原子电荷"""
    try:
        return float(atom.GetProp('_GasteigerCharge'))
    except:
        return 0.0

print("=" * 60)
print("1. 非共价相互作用位点识别")
print("=" * 60)

# 识别可能的非共价相互作用位点
nci_sites = {
    "氢键受体": [],
    "氢键供体": [],
    "π体系": [],
    "孤对电子": []
}

# 分析原子
for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    symbol = atom.GetSymbol()
    charge = getCharge(atom)
    
    # 氢键受体 (电负性原子，负电荷)
    if symbol in ['N', 'O', 'F'] and charge < -0.05:
        nci_sites["氢键受体"].append((idx, symbol, charge))
    
    # 氢键供体 (连接 H 的 N, O)
    if symbol in ['N', 'O'] and any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors()):
        nci_sites["氢键供体"].append((idx, symbol, charge))
    
    # 孤对电子 (N, O, S 的非键电子)
    if symbol in ['N', 'O', 'S']:
        nci_sites["孤对电子"].append((idx, symbol, charge))
    
    # π体系 (芳香碳)
    if atom.GetIsAromatic():
        nci_sites["π体系"].append((idx, symbol, charge))

# 打印结果
for site_type, sites in nci_sites.items():
    if sites:
        print(f"\n{site_type}:")
        for idx, symbol, charge in sites[:8]:  # 显示前8个
            print(f"  {symbol} (idx {idx}): 电荷 {charge:.3f}")
        if len(sites) > 8:
            print(f"  ... 共 {len(sites)} 个")

print("\n" + "=" * 60)
print("2. π-π 堆积倾向分析")
print("=" * 60)

# 计算芳香环的几何特征
ring_info = mol_3d.GetRingInfo()
aromatic_rings = [ring for ring in ring_info.AtomRings() if all(mol_3d.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

print(f"芳香环数量: {len(aromatic_rings)}")

# 计算每个环的中心和法向量
conf = mol_3d.GetConformer()
ring_centers = []
ring_normals = []

for ring in aromatic_rings:
    coords = np.array([list(conf.GetAtomPosition(i)) for i in ring])
    center = coords.mean(axis=0)
    ring_centers.append(center)
    
    # 计算法向量（使用前3个原子）
    if len(ring) >= 3:
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        normal = np.cross(v1, v2)
        normal = normal / np.linalg.norm(normal)
        ring_normals.append(normal)

# 分析环之间的夹角（π-π 堆积倾向）
print("\n芳香环相对取向:")
for i in range(len(ring_normals)):
    for j in range(i+1, len(ring_normals)):
        angle = np.degrees(np.arccos(np.abs(np.dot(ring_normals[i], ring_normals[j]))))
        distance = np.linalg.norm(ring_centers[i] - ring_centers[j])
        
        if distance < 6.0:  # 6 Å 以内的环
            print(f"  环 {i} - 环 {j}: 夹角 {angle:.1f}°, 距离 {distance:.2f} Å")
            if angle > 160 or angle < 20:
                print(f"    → 可能的 π-π 堆积!")

print("\n" + "=" * 60)
print("3. 静电势表面可视化 (原子电荷着色)")
print("=" * 60)

# 生成原子颜色映射（基于电荷）
charges = [getCharge(atom) for atom in mol_3d.GetAtoms()]
charge_min, charge_max = min(charges), max(charges)

print(f"电荷范围: [{charge_min:.3f}, {charge_max:.3f}]")
print("颜色映射: 红色(负) → 白色(中性) → 蓝色(正)")

# 创建原子颜色
atom_colors = {}
for i, atom in enumerate(mol_3d.GetAtoms()):
    charge = charges[i]
    # 归一化到 [0, 1]
    if charge_max != charge_min:
        norm_charge = (charge - charge_min) / (charge_max - charge_min)
    else:
        norm_charge = 0.5
    
    # RGB: 红(负) -> 白(中性) -> 蓝(正)
    if norm_charge < 0.5:
        # 负电荷：红色
        r = 1.0
        g = norm_charge * 2
        b = norm_charge * 2
    else:
        # 正电荷：蓝色
        r = (1 - norm_charge) * 2
        g = (1 - norm_charge) * 2
        b = 1.0
    
    atom_colors[i] = (r, g, b, 1.0)

# 生成可视化
drawer = rdMolDraw2D.MolDraw2DCairo(1000, 800)
drawer.SetFontSize(0.9)
drawer.DrawMolecule(mol_3d, highlightAtoms=list(atom_colors.keys()), highlightAtomColors=atom_colors)
drawer.FinishDrawing()
drawer.WriteDrawingText(f"{NAME}_charge_surface.png")
print(f"✓ 保存: {NAME}_charge_surface.png")

print("\n" + "=" * 60)
print("4. 非共价相互作用位点可视化")
print("=" * 60)

# 标记关键位点
hb_acceptors = [site[0] for site in nci_sites["氢键受体"]]
pi_centers = []

# 选择芳香环中心原子
for i, ring in enumerate(aromatic_rings):
    if len(ring) > 0:
        # 选择环的中心原子
        center_atom = ring[len(ring)//2]
        pi_centers.append(center_atom)

# 红色标记氢键受体，蓝色标记π体系
atom_colors2 = {}
for idx in hb_acceptors:
    atom_colors2[idx] = (1.0, 0.2, 0.2, 0.8)  # 红色

for idx in pi_centers:
    if idx not in atom_colors2:
        atom_colors2[idx] = (0.2, 0.4, 1.0, 0.8)  # 蓝色

# 生成可视化
drawer2 = rdMolDraw2D.MolDraw2DCairo(1000, 800)
drawer2.SetFontSize(0.9)
drawer2.DrawMolecule(mol_3d, highlightAtoms=list(atom_colors2.keys()), highlightAtomColors=atom_colors2)
drawer2.FinishDrawing()
drawer2.WriteDrawingText(f"{NAME}_nci_sites.png")
print(f"✓ 保存: {NAME}_nci_sites.png")
print(f"  红色: 氢键受体 ({len(hb_acceptors)} 个)")
print(f"  蓝色: π体系中心 ({len(pi_centers)} 个)")

print("\n" + "=" * 60)
print("5. 分子堆积倾向预测")
print("=" * 60)

# 计算分子形状因子
coords_3d = np.array([list(conf.GetAtomPosition(i)) for i in range(mol_3d.GetNumAtoms())])
cov_matrix = np.cov((coords_3d - coords_3d.mean(axis=0)).T)
eigenvalues = np.real(np.sort(np.linalg.eigvals(cov_matrix))[::-1])

L1, L2, L3 = np.sqrt(eigenvalues)
aspect_ratios = (L1/L3, L2/L3)

print(f"分子形状因子:")
print(f"  长宽比 (L1/L3): {aspect_ratios[0]:.2f}")
print(f"  宽厚比 (L2/L3): {aspect_ratios[1]:.2f}")

if aspect_ratios[0] > 2 and aspect_ratios[1] > 1.5:
    print(f"  → 扁平形状，有利于面-面堆积")
elif aspect_ratios[0] > 3:
    print(f"  → 长条形状，可能形成边-面堆积")
else:
    print(f"  → 较为球形，堆积倾向较弱")

# 计算偶极矩（近似）
print(f"\n偶极矩估算:")
dipole = np.zeros(3)
for i, atom in enumerate(mol_3d.GetAtoms()):
    charge = getCharge(atom)
    pos = np.array(list(conf.GetAtomPosition(i)))
    dipole += charge * pos

dipole_magnitude = np.linalg.norm(dipole)
print(f"  |μ| ≈ {dipole_magnitude:.2f} Debye (Gasteiger 电荷估算)")
if dipole_magnitude > 3:
    print(f"  → 较大偶极矩，可能影响固态堆积")

print("\n" + "=" * 60)
print("✅ 完成！")
print("=" * 60)
print(f"生成的可视化:")
print(f"  - {NAME}_charge_surface.png (原子电荷着色)")
print(f"  - {NAME}_nci_sites.png (非共价相互作用位点)")
print("\n建议: 使用 Gaussian + Multiwfn 进行精确的 ESP 和 RDG 分析")
print("  Gaussian: #P B3LYP/6-31G* SCF=Tight")
print("  Multiwfn: 功能 5 (ESP), 功能 20 (RDG)")
print("=" * 60)
