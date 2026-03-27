# Multiwfn Benzene Reference

Complete workflow for analyzing benzene using Multiwfn (v3.8), from Gaussian output to spectral plots.

## Prerequisites

Run a Gaussian DFT calculation first:
```bash
# benzene.com
%NProcShared=8
%Mem=8GB
#P B3LYP/cc-pVDZ Opt Freq TD(NStates=10)

Benzene DFT + TDDFT

0 1
C   0.000000   1.395000   0.000000
C   1.208543   0.697500   0.000000
C   1.208543  -0.697500   0.000000
C   0.000000  -1.395000   0.000000
C  -1.208543  -0.697500   0.000000
C  -1.208543   0.697500   0.000000
H   0.000000   2.479000   0.000000
H   2.150000   1.239500   0.000000
H   2.150000  -1.239500   0.000000
H   0.000000  -2.479000   0.000000
H  -2.150000  -1.239500   0.000000
H  -2.150000   1.239500   0.000000

```

Convert to `.fchk`:
```bash
formchk benzene.chk benzene.fchk
```

## Function 0: Orbital Visualization

```
Multiwfn benzene.fchk
0     # Orbital information
-1    # Show all orbital energies
5     # Generate isosurface → HOMO (orbital 15)
      # isovalue = 0.05 → benzene_homo.cub
5     # Generate isosurface → LUMO (orbital 17)
      # isovalue = 0.05 → benzene_lumo.cub
```

## Function 7: Population Analysis

### Mulliken
```
7 → 1 → 1    # Mulliken charges
```

Expected: C ~-0.089, H ~+0.089 (all equivalent by symmetry).

### Hirshfeld
```
7 → 2 → 1
```

### ADCH (recommended)
```
7 → 10 → 1
```

## Function 8: Orbital Composition

```
8 → 15 → 1 → 1    # HOMO composition by atom (Mulliken)
```

Benzene HOMO (π, e1g): each C contributes ~16.7%, H contributes ~0%.

## Function 9: Bond Order Analysis

### Wiberg/Mayer
```
9 → 1 → 1
```

Expected: C-C ~1.389, C-H ~0.962.

### Multi-center (aromaticity)
```
9 → 7 → 6    # 6-center bond order (select all 6 C atoms)
```

## Function 10: DOS/PDOS

```
10 → 1         # Total DOS, range -15 to 5 eV, σ=0.15 eV
10 → 2 → 1     # PDOS by atom (C vs H)
10 → 2 → 2     # PDOS by orbital type (s vs p)
```

## Function 11: UV-Vis Spectrum

```
11 → 1    # UV-Vis, σ=0.15 eV, range 2-8 eV, 500 points
```

Expected peaks:
| Peak | Energy (eV) | λ (nm) | Assignment |
|------|-------------|--------|------------|
| 1 | ~4.7 | ~264 | π→π* (1E1u) |
| 2 | ~6.2 | ~200 | π→π* (1B1u) |
| 3 | ~6.9 | ~180 | π→π* (1E2u) |

## Function 18: Excited State Analysis

```
18 → 1    # Transition contributions for state 1
```

## Function 20: RDG Weak Interaction

```
20 → 1    # RDG isosurface → benzene_rdg.cub
```

## Automation Script

```bash
#!/bin/bash
FCHK="benzene.fchk"

echo "7\n1\n1\n" | Multiwfn $FCHK > mulliken.log 2>&1
echo "9\n1\n1\n" | Multiwfn $FCHK > bondorder.log 2>&1
echo "11\n1\n0.15\n2\n8\n500\n" | Multiwfn $FCHK > uvvis.log 2>&1

echo "Done!"
```
