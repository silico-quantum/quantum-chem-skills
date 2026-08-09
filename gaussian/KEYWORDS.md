# Conservative Gaussian keyword orientation

Use this page only to locate a keyword family. The installed Gaussian manual is
the authority for syntax, defaults, availability, and interactions. Do not copy
a keyword into production input solely because it appears here.

| Need | Common family to verify | Acceptance consequence |
|---|---|---|
| Single-point energy | `SP` or the route default | Require SCF and normal termination |
| Geometry optimization | `Opt` | Require optimization convergence |
| Harmonic frequencies | `Freq` | Inspect negative modes and thermochemical conditions |
| Transition-state search | `Opt=TS` | Require one intended negative mode; use IRC if needed |
| Reaction path | `IRC` | Verify both directions and endpoint identities |
| Excited states | `TD` | Verify state manifold, root, and state character |
| Continuum solvent | `SCRF` | Record model and solvent printed by Gaussian |
| Population analysis | `Pop` | Report scheme and basis dependence |
| Wavefunction stability | `Stable` | Distinguish testing from `Stable=Opt` modification |
| User basis or ECP | `Gen`, `GenECP`, `Pseudo` | Validate element coverage and section layout |
| Geometry from checkpoint | `Geom=Check` or `Geom=AllCheck` | Verify which molecule fields are read |
| Guess from checkpoint | `Guess=Read` | Verify checkpoint compatibility |
| SCF recovery | `SCF` options | Diagnose first; change one justified option at a time |

Frequently used Link 0 commands include `%Chk`, `%OldChk`, `%Mem`, and
`%NProcShared`. Their values must agree with the actual scheduler allocation.

Use the official pages for the relevant installed revision:

- [Keyword index](https://gaussian.com/keywords/)
- [Link 0 commands](https://gaussian.com/link0/)
- [Optimization](https://gaussian.com/opt/)
- [Geometry input](https://gaussian.com/geom/)
- [Frequency calculations](https://gaussian.com/freq/)
- [TD excited states](https://gaussian.com/td/)
- [SCRF solvent models](https://gaussian.com/scrf/)
