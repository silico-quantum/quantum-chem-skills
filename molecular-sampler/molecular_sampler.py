#!/usr/bin/env python3
"""
Molecular Sampler - Extract monomers and sample multi-molecule complexes
from Gaussian ONIOM or XYZ files.
"""

import re
import numpy as np
from collections import defaultdict, Counter
import os
import argparse
import sys

# Covalent radii (in Angstroms)
COVALENT_RADII = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 
    'S': 1.05, 'P': 1.07, 'F': 0.71, 'Cl': 0.99,
    'Br': 1.14, 'I': 1.33
}

def parse_gjf_file(filepath):
    """Parse Gaussian GJF file and extract atoms"""
    atoms = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    in_geometry = False
    for line in lines:
        line = line.strip()
        if not line:
            continue
        
        # Check for charge/multiplicity line (ONIOM: 0 1 0 1 0 1)
        if re.match(r'^\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+$', line):
            in_geometry = True
            continue
        
        # Standard format: 0 1 (single line)
        if re.match(r'^\d+\s+\d+$', line) and not in_geometry:
            in_geometry = True
            continue
        
        if in_geometry:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    element = parts[0]
                    
                    # Check if ONIOM format (has layer marker)
                    if len(parts) >= 5 and parts[1] in ['0', '-1', 'H', 'L']:
                        layer = 'H' if parts[1] in ['0', 'H'] else 'L'
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                    else:
                        # Standard XYZ format
                        layer = 'L'  # Default to L layer
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                    
                    atoms.append({
                        'element': element,
                        'layer': layer,
                        'x': x,
                        'y': y,
                        'z': z
                    })
                except (ValueError, IndexError):
                    continue
    
    return atoms

def parse_xyz_file(filepath):
    """Parse standard XYZ file"""
    atoms = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Skip first two lines (atom count and comment)
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >= 4:
            try:
                element = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                
                atoms.append({
                    'element': element,
                    'layer': 'L',  # Default to L layer
                    'x': x,
                    'y': y,
                    'z': z
                })
            except (ValueError, IndexError):
                continue
    
    return atoms

def bond_distance(e1, e2):
    """Calculate maximum bond distance between two elements"""
    r1 = COVALENT_RADII.get(e1, 0.77)
    r2 = COVALENT_RADII.get(e2, 0.77)
    return (r1 + r2) * 1.3

class UnionFind:
    """Union-Find data structure for molecular clustering"""
    
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n
    
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    
    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px != py:
            if self.rank[px] < self.rank[py]:
                px, py = py, px
            self.parent[py] = px
            if self.rank[px] == self.rank[py]:
                self.rank[px] += 1

def identify_molecules(atoms, layer_filter='L'):
    """Identify individual molecules using bond detection"""
    
    # Filter atoms by layer
    filtered_atoms = [a for a in atoms if layer_filter == 'all' or a['layer'] == layer_filter]
    
    if not filtered_atoms:
        return []
    
    # Build molecular graph using Union-Find
    uf = UnionFind(len(filtered_atoms))
    
    print(f"Building molecular graph ({len(filtered_atoms)} atoms)...")
    for i in range(len(filtered_atoms)):
        if i % 500 == 0 and i > 0:
            print(f"  Progress: {i}/{len(filtered_atoms)}")
        
        for j in range(i+1, len(filtered_atoms)):
            a1, a2 = filtered_atoms[i], filtered_atoms[j]
            
            # Calculate distance
            dist = np.sqrt((a1['x']-a2['x'])**2 + 
                          (a1['y']-a2['y'])**2 + 
                          (a1['z']-a2['z'])**2)
            
            # Check if within bond distance
            max_dist = bond_distance(a1['element'], a2['element'])
            if dist < max_dist:
                uf.union(i, j)
    
    # Group atoms by component
    components = defaultdict(list)
    for i in range(len(filtered_atoms)):
        components[uf.find(i)].append(filtered_atoms[i])
    
    # Create molecule objects
    molecules = []
    for comp_id, comp_atoms in components.items():
        elements = Counter(a['element'] for a in comp_atoms)
        
        # Filter out very small fragments (< 5 atoms)
        if len(comp_atoms) >= 5:
            center = np.mean([[a['x'], a['y'], a['z']] for a in comp_atoms], axis=0)
            molecules.append({
                'id': comp_id,
                'atoms': comp_atoms,
                'elements': dict(elements),
                'count': len(comp_atoms),
                'center': center
            })
    
    # Sort molecules by center position
    molecules.sort(key=lambda m: (m['center'][0], m['center'][1], m['center'][2]))
    
    return molecules

def write_xyz(mol_list, filename, comment=""):
    """Write molecules to XYZ file"""
    all_atoms = []
    for mol in mol_list:
        all_atoms.extend(mol['atoms'])
    
    with open(filename, 'w') as f:
        f.write(f"{len(all_atoms)}\n")
        f.write(f"{comment}\n")
        for atom in all_atoms:
            f.write(f"{atom['element']:2s} {atom['x']:12.6f} {atom['y']:12.6f} {atom['z']:12.6f}\n")

def sample_molecules(molecules, n_samples=20):
    """Sample monomers and multi-molecule complexes"""
    
    # Calculate pairwise distances
    print("Calculating pairwise distances...")
    n_mols = len(molecules)
    distances = {}
    
    for i in range(n_mols):
        for j in range(i+1, n_mols):
            dist = np.linalg.norm(molecules[i]['center'] - molecules[j]['center'])
            distances[(i, j)] = dist
    
    # Build sorted neighbor lists
    print("Building sorted neighbor lists...")
    neighbors = {i: [] for i in range(n_mols)}
    
    for (i, j), dist in distances.items():
        neighbors[i].append((j, dist))
        neighbors[j].append((i, dist))
    
    # Sort neighbors by distance
    for i in range(n_mols):
        neighbors[i].sort(key=lambda x: x[1])
    
    # Sample complexes
    samples = {
        'monomers': [],
        'dimers': [],
        'trimers': [],
        'tetramers': [],
        'pentamers': []
    }
    
    # All monomers
    for i, mol in enumerate(molecules):
        samples['monomers'].append([mol])
    
    # Multi-molecule complexes
    for n_mols_target, name in [(2, 'dimers'), (3, 'trimers'), 
                                 (4, 'tetramers'), (5, 'pentamers')]:
        sample_count = 0
        
        for start_idx in range(min(n_samples, n_mols)):
            selected_indices = [start_idx]
            
            # Take nearest neighbors
            for neighbor_idx, dist in neighbors[start_idx][:n_mols_target-1]:
                selected_indices.append(neighbor_idx)
            
            mol_list = [molecules[i] for i in selected_indices]
            samples[name].append(mol_list)
            sample_count += 1
            
            if sample_count >= n_samples:
                break
    
    return samples

def main():
    parser = argparse.ArgumentParser(description='Sample molecular structures from GJF/XYZ files')
    parser.add_argument('input_file', help='Input GJF or XYZ file')
    parser.add_argument('--output-dir', default='./molecular_samples', 
                       help='Output directory (default: ./molecular_samples)')
    parser.add_argument('--samples', type=int, default=20,
                       help='Number of samples per complex type (default: 20)')
    parser.add_argument('--layer', choices=['H', 'L', 'all'], default='L',
                       help='Layer to sample: H (high), L (low), or all (default: L)')
    
    args = parser.parse_args()
    
    # Parse input file
    print(f"Parsing {args.input_file}...")
    
    if args.input_file.endswith('.gjf'):
        atoms = parse_gjf_file(args.input_file)
    elif args.input_file.endswith('.xyz'):
        atoms = parse_xyz_file(args.input_file)
    else:
        print("Error: Unsupported file format. Use .gjf or .xyz")
        sys.exit(1)
    
    print(f"Total atoms: {len(atoms)}")
    
    # Identify molecules
    molecules = identify_molecules(atoms, layer_filter=args.layer)
    print(f"\nFound {len(molecules)} molecules")
    
    if not molecules:
        print("Error: No molecules found")
        sys.exit(1)
    
    # Sample molecules
    samples = sample_molecules(molecules, n_samples=args.samples)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save samples
    for sample_type, mol_lists in samples.items():
        sample_dir = os.path.join(args.output_dir, sample_type)
        os.makedirs(sample_dir, exist_ok=True)
        
        for i, mol_list in enumerate(mol_lists, 1):
            filename = os.path.join(sample_dir, f'{sample_type[:-1]}_{i:02d}.xyz')
            total_atoms = sum(m['count'] for m in mol_list)
            comment = f"{sample_type[:-1]} {i}: {len(mol_list)} molecules, {total_atoms} atoms"
            write_xyz(mol_list, filename, comment)
        
        print(f"  Saved {len(mol_lists)} {sample_type}")
    
    # Create summary
    summary_path = os.path.join(args.output_dir, 'sampling_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("=== Molecular Sampling Summary ===\n\n")
        f.write(f"Input file: {args.input_file}\n")
        f.write(f"Layer: {args.layer}\n")
        f.write(f"Total molecules: {len(molecules)}\n\n")
        
        f.write("Sampling Results:\n")
        f.write(f"  - Monomers: {len(samples['monomers'])} files\n")
        f.write(f"  - Dimers: {len(samples['dimers'])} samples\n")
        f.write(f"  - Trimers: {len(samples['trimers'])} samples\n")
        f.write(f"  - Tetramers: {len(samples['tetramers'])} samples\n")
        f.write(f"  - Pentamers: {len(samples['pentamers'])} samples\n")
    
    print(f"\n=== Summary ===")
    print(f"All files saved to: {args.output_dir}")
    print(f"Summary: {summary_path}")

if __name__ == '__main__':
    main()
