# Gaussian Input Templates

These templates preserve the common task patterns summarized in
[`SKILL.md`](SKILL.md). Replace every placeholder, choose methods and basis
sets for the target system and property, and confirm version-specific syntax
against the licensed Gaussian documentation before submitting production jobs.

## Ground-State Optimization and Frequency

```
%chk=opt_freq.chk
%mem=60GB
%nproc=64
# opt freq b3lyp/6-31G**

Molecule optimization
0 1
C  0.0  0.0  0.0
...
```

## TDDFT Excited States

```
%chk=tddft.chk
%mem=60GB
%nproc=64
# td=(nstates=10,singlets) b3lyp/6-31G*

S1 vertical excitations
0 1
C  0.0  0.0  0.0
...
```

## S1 Optimization and Vertical Emission

```
%chk=s1_opt.chk
%mem=60GB
%nproc=64
# td=(nstates=5,root=1) b3lyp/6-31G* opt

S1 optimization
0 1
... (ground-state optimized geometry)
```

```
%chk=emission.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* td=(root=1) geom=check

S1 vertical emission
0 1
... (S1 optimized geometry)
```

## Solvent Effect (SMD)

```
%chk=smd.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* scrf=(smd,solvent=acetonitrile)

SMD solvent model
0 1
...
```

## CBS-QB3 Composite Method

```
%chk=cbs-qb3.chk
%mem=60GB
%nproc=64
# cbs-qb3

CBS-QB3 calculation
0 1
...
```

## CASSCF

```
%chk=cas.chk
%mem=60GB
%nproc=64
# casSCF(10,8)/6-31G*

CAS(10,8) — 10 electrons, 8 orbitals in active space
0 1
...
```

## ModRedundant Scan

```
%chk=scan.chk
%mem=60GB
%nproc=64
# b3lyp/6-31G* opt=modredundant

Dihedral scan
0 1
C 0 0 0
C 1 2 1.4
...
1 2 3 4 90.0 S 20 5.0   # Scan dihedral 1-2-3-4, 20 steps, 5° each
```
