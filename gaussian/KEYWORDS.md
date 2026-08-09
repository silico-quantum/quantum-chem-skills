# Gaussian Keyword Reference

This page is a compact orientation guide, not a substitute for the
[official Gaussian keyword documentation](https://gaussian.com/keywords/).
Confirm version-specific syntax before production calculations.

## Job Types

| Keyword | Purpose | Common options |
|---------|---------|----------------|
| SP | Single-point energy | Default |
| Opt | Geometry optimization | `Opt=CalcFC`, `Opt=TS` |
| Freq | Frequency analysis | `Freq=Raman` |
| Scan | Potential-energy-surface scan | `Scan=Opt` |
| IRC | Intrinsic reaction coordinate | `IRC=CalcFC` |
| Polar | Polarizability | `Polar=EnOnly` |

## DFT Functionals

| Functional | Type | Characteristic |
|------------|------|----------------|
| B3LYP | Hybrid | 20% Hartree–Fock exchange |
| PBE0 | Hybrid | 25% Hartree–Fock exchange |
| M06-2X | Meta-hybrid GGA | 54% Hartree–Fock exchange |
| CAM-B3LYP | Range-separated hybrid | Long-range correction |
| ωB97X-D | Range-separated hybrid with dispersion | Long-range correction and dispersion |

Method suitability is system- and property-dependent; benchmark against an
appropriate reference rather than choosing a functional from this table alone.

## Basis Sets

### Pople family

- `6-31G`: split-valence double-zeta
- `6-31G*`: double-zeta with polarization on non-hydrogen atoms
- `6-31G**`: double-zeta with polarization on all atoms
- `6-31+G*`: double-zeta with diffuse functions and polarization
- `6-311G**`: split-valence triple-zeta with polarization

### Correlation-consistent family

- `cc-pVDZ`, `cc-pVTZ`, `cc-pVQZ`
- `aug-cc-pVTZ` with diffuse augmentation

### def2 family

- `def2-SVP`, `def2-TZVP`, `def2-QZVP`

## TDDFT Excited States

```
TD(NStates=10,Root=1,Singlets)
```

- `NStates`: number of excited states
- `Root`: state selected for state-specific operations such as optimization
- `Singlets` / `Triplets`: requested spin manifold

## Solvent Effects (SCRF)

```
SCRF=(SMD,Solvent=Water)
```

Common solvent names include `Water`, `Ethanol`, `Acetonitrile`, `DMF`,
`DMSO`, and `Toluene`.

## Optimization Options

- `Opt=CalcFC`: calculate force constants at the initial geometry
- `Opt=TS`: search for a transition state
- `Opt=Tight`: use tighter convergence criteria
- `Opt=MaxStep=N`: set the maximum optimization step

## SCF Convergence

- `SCF=QC`: use quadratic convergence
- `SCF=XQC`: start with the default procedure and fall back to quadratic convergence
- `SCF=Tight`: use tighter SCF convergence criteria

## Population Analysis (Pop)

- `Pop=Full`: request expanded orbital information
- `Pop=NBO`: request Natural Bond Orbital analysis when the required support is available
- `Pop=NPA`: request Natural Population Analysis when supported
- `Pop=MK`: fit Merz–Kollman electrostatic-potential charges

## Link 0 Commands

- `%chk=file.chk`: checkpoint file
- `%mem=60GB`: memory allocation
- `%nprocshared=64`: shared-memory worker count
