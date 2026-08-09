---
name: multiwfn
description: Use when analyzing molecular wavefunctions with Multiwfn, including orbital composition, population and bond-order analysis, spectra, density topology, excited states, and weak interactions.
license: MIT
compatibility: Requires a local Multiwfn installation and a supported wavefunction file; remote-copy examples assume an OpenSSH-compatible scp command.
---

# Multiwfn Wavefunction Analysis

## Overview

Multiwfn is a wavefunction-analysis program for molecular-orbital analysis,
spectrum generation, electron-density topology, and related tasks.

## Installation and executable discovery

Locate an existing installation before running a workflow:

```bash
command -v Multiwfn
```

If the executable is not already on `PATH`, point `MULTIWFN_HOME` at the
installation directory and add it explicitly:

```bash
export MULTIWFN_HOME="/path/to/Multiwfn_3.8"
export PATH="${MULTIWFN_HOME}:${PATH}"
command -v Multiwfn
```

This avoids relying on a machine-specific Homebrew or system path. The menu
sequences in this repository were recorded for Multiwfn 3.8; verify them
against the installed release. Multiwfn is developed by Tian Lu at the Beijing
Kein Research Center for Natural Sciences.

## Citation requirements

Publications that use Multiwfn should cite:

1. Tian Lu and Feiwu Chen, *J. Comput. Chem.* **33**, 580 (2012),
   DOI: [10.1002/jcc.22885](https://doi.org/10.1002/jcc.22885).
2. Tian Lu, *J. Chem. Phys.* **161**, 082503 (2024),
   DOI: [10.1063/5.0216272](https://doi.org/10.1063/5.0216272).

## Supported input files

| Format | Description |
|--------|-------------|
| `.fchk` / `.fch` | Gaussian formatted checkpoint file; recommended |
| `.wfn` | Wavefunction file |
| `.molden` | Molden-format file |
| `.xyz` | Molecular coordinates |
| `.cub` / `.cube` | Gaussian cube file |

Convert a binary Gaussian `.chk` file to `.fchk` with the matching `formchk`
utility before analysis. Coordinate-only and scalar-grid inputs do not contain
the full wavefunction information required by every menu function.

## Main-menu functions

### 0 - Molecular structure and orbitals

Open the interactive 3D molecular viewer and inspect molecular orbitals.

### 1 - Point-property calculation

Calculate properties at a point in space, including electron density,
gradients, and the Laplacian.

### 2 - Topology analysis (AIM)

Analyze electron-density critical points (CPs) for chemical-bond analysis.

### 3 - Line-property plot

Plot how a property changes along a line, such as the electron density along a
bond axis.

### 4 - Plane plot

Draw a contour plot on a selected plane, for example for electron density or
an orbital.

### 5 - Spatial grid data

Calculate three-dimensional grid data and generate a cube file.

### 6 - Inspect or modify a wavefunction

Inspect orbital information and energies, or modify the wavefunction.

### 7 - Population analysis

Calculate atomic charges with Mulliken, Hirshfeld, ADCH, and related schemes.

### 8 - Orbital-composition analysis

Analyze atomic or fragment contributions to molecular orbitals.

### 9 - Bond-order analysis

Calculate Mayer, Wiberg, and other bond-order measures.

### 10 - DOS and PDOS

Generate the density of states, projected density of states, or a
photoelectron spectrum.

### 11 - Spectrum generation

Generate IR, Raman, UV-Vis, ECD, and NMR spectra.

### 12 - Molecular-surface analysis

Calculate electrostatic potential, average local ionization energy, and other
properties on a molecular surface.

### 18 - Excited-state analysis

Analyze electronic excitations, hole-electron distributions, and charge
transfer.

### 20 - Weak-interaction visualization

Perform reduced-density-gradient (RDG) and noncovalent-interaction (NCI)
analyses to visualize intermolecular interactions.

## Common menu paths

| Task | Menu path |
|------|-----------|
| Inspect orbital energies | `0 -> q` |
| UV-Vis spectrum | `11 -> 1` |
| IR spectrum | `11 -> 2` |
| Atomic charges | `7 -> 1/2/3` |
| Orbital composition | `8 -> 1` |
| Bond order | `9 -> 1` |
| DOS | `10 -> 1` |
| Molecular surface | `12` |
| RDG analysis | `20 -> 1` |
| Excited-state analysis | `18` |

## Gaussian-to-Multiwfn workflow

1. Complete the DFT or TDDFT calculation with Gaussian.
2. Convert the `.chk` file to `.fchk` with `formchk`.
3. Analyze the `.fchk` file with Multiwfn.

## Copying a file from a remote server

Set the connection values in the shell rather than embedding private host or
account details in documentation:

```bash
: "${QC_USER:?Set QC_USER}"
: "${QC_HOST:?Set QC_HOST}"
: "${QC_PORT:?Set QC_PORT}"

scp -P "${QC_PORT}" \
  "${QC_USER}@${QC_HOST}:/path/to/molecule.fchk" ./
```

## Worked example

- [Benzene analysis workflow](references/benzene-analysis.md)

## Official resources

- [Multiwfn website](http://sobereva.com/multiwfn)
- [English-language Multiwfn forum](http://sobereva.com/wfnbbs)
- [Chinese-language wavefunction-analysis forum](http://bbs.keinsci.com/wfn)
