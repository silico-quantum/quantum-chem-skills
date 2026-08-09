# Gaussian input templates

These are review templates, not recommended scientific protocols. Replace every
angle-bracket placeholder, choose the method and basis/ECP from a documented
protocol, and check syntax against the installed Gaussian revision. Each block
must end with a blank line in the actual input file.

## Five-stage S0-to-emission chain

Use a new log and checkpoint at every stage. `%OldChk` is the immutable accepted
source; `%Chk` is the new destination. Do not start a later stage until the
source stage passes its termination, convergence, and state-identity gates.

### 1. S0 optimization and frequency

```text
%Chk=<s0-optfreq-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> Opt Freq

<title>

<charge> <multiplicity>
<element> <x> <y> <z>
<remaining Cartesian coordinates>

```

Accept an S0 minimum only after normal termination, optimization convergence,
and frequency-mode inspection under a declared numerical-noise threshold.

### 2. Vertical excitation at the accepted S0 geometry

```text
%OldChk=<s0-optfreq-run-id>.chk
%Chk=<vertical-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> TD=(NStates=<count>,<spin-manifold>) Geom=Check Guess=Read

<title>

<charge> <multiplicity>

```

Characterize candidate states using excitation energy, oscillator strength,
dominant configurations or NTOs, symmetry when meaningful, and other relevant
state descriptors. Select a target state scientifically before optimization.

### 3. Target-state optimization

```text
%OldChk=<vertical-run-id>.chk
%Chk=<s1-opt-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> TD=(NStates=<count>,Root=<root>,<spin-manifold>) Opt Geom=Check Guess=Read

<title>

<charge> <multiplicity>

```

Track state character at every optimization step. The same root number does not
prove the same state. Record geometry step, excitation energy, oscillator
strength, leading configuration/NTO character, symmetry or localization
descriptor, and the pass/fail state-continuity decision.

### 4. Target-state frequency, when supported

Confirm that the installed revision supports the requested excited-state
frequency for the exact method and state. If it does:

```text
%OldChk=<s1-opt-run-id>.chk
%Chk=<s1-freq-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> TD=(NStates=<count>,Root=<root>,<spin-manifold>) Freq Geom=Check Guess=Read

<title>

<charge> <multiplicity>

```

Only a compatible accepted frequency analysis can support the label “verified
S1 minimum.” If it is unsupported or not run, report the result as an
`optimized S1 stationary geometry` with `minimum_not_verified`; do not infer a
minimum from optimization convergence alone.

### 5. Vertical emission at the accepted S1 geometry

```text
%OldChk=<s1-opt-run-id>.chk
%Chk=<emission-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> TD=(NStates=<count>,Root=<root>,<spin-manifold>) Geom=Check Guess=Read

<title>

<charge> <multiplicity>

```

For S1-to-S0 fluorescence, the target state's printed TD excitation energy at
the accepted S1 geometry is the vertical emission energy under this response
model. Report that value, its printed units and one controlled conversion,
oscillator strength, and state-character evidence. Do not subtract raw SCF
reference energies from different jobs to construct the emission energy.

## State-continuity ledger

For stages 2 through 5, retain a table like:

```text
stage,geometry_step,root,excitation_eV,oscillator_strength,leading_character,state_continuity,status
```

The exact character descriptor depends on the system. A discontinuous change
requires review or restart; normal termination does not override a root flip.

## Solvent calculation

```text
%Chk=<run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> <job-type> SCRF=(<model>,Solvent=<solvent>)

<title>

<charge> <multiplicity>
<element> <x> <y> <z>
<remaining Cartesian coordinates>

```

Record the exact solvent model and solvent identifier printed in the output.

## Relaxed ModRedundant dihedral scan

```text
%Chk=<scan-run-id>.chk
%Mem=<allocated-memory>
%NProcShared=<allocated-shared-memory-cores>
#P <method>/<basis> Opt=ModRedundant

<title>

<charge> <multiplicity>
<element> <x> <y> <z>
<remaining Cartesian coordinates>

D <atom-1> <atom-2> <atom-3> <atom-4> S <steps> <step-size-degrees>

```

For example, `D 1 2 3 4 S 20 5.0` requests 20 relaxed scan steps of
5 degrees from the input dihedral. Atom indices are one-based. Do not insert an
unverified starting angle before `S`; alter or freeze coordinates only with the
documented action syntax for the installed revision.

## Checkpoint conversion

Run conversion only after the source job is accepted:

```bash
test ! -e molecule.fchk && test ! -L molecule.fchk || exit 1
formchk molecule.chk molecule.fchk
```

This converts a binary checkpoint to a formatted checkpoint. The reverse is:

```bash
test ! -e molecule-restored.chk && \
  test ! -L molecule-restored.chk || exit 1
unfchk molecule.fchk molecule-restored.chk
```

This converts a formatted checkpoint to a machine-local binary checkpoint.
