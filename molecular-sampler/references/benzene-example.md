# Benzene Dimer Sampling Example

This example is a parser and sampling smoke test, not a claim that the supplied
geometry is an optimized benzene dimer.

## Input

Save the following single-frame, angstrom XYZ as `benzene_dimer.xyz`:

```text
24
Two separated benzene fragments; angstrom
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
C   0.000000   7.395000   0.000000
C   1.208543   6.697500   0.000000
C   1.208543   5.302500   0.000000
C   0.000000   4.605000   0.000000
C  -1.208543   5.302500   0.000000
C  -1.208543   6.697500   0.000000
H   0.000000   8.479000   0.000000
H   2.150000   7.239500   0.000000
H   2.150000   4.760500   0.000000
H   0.000000   3.521000   0.000000
H  -2.150000   4.760500   0.000000
H  -2.150000   7.239500   0.000000
```

The expected distance-inferred partition is two 12-atom `C6H6` fragments.

## Run

Use a new output path:

```bash
python3 molecular_sampler.py benzene_dimer.xyz \
  --xyz-units angstrom \
  --layer all \
  --expected-fragments 2 \
  --bond-scale 1.3 \
  --samples 20 \
  --output-dir benzene_dimer_samples_run01
```

Despite `--samples 20`, two source fragments can produce only two monomers and
one unique nearest-neighbor dimer. Trimer, tetramer, and pentamer directories
must be empty. Claiming 20 dimers from this input would be incorrect.

## Acceptance checks

```bash
test -s benzene_dimer_samples_run01/sampling_manifest.json
test -s benzene_dimer_samples_run01/sampling_summary.txt
test "$(find benzene_dimer_samples_run01/monomers -name '*.xyz' | wc -l)" -eq 2
test "$(find benzene_dimer_samples_run01/dimers -name '*.xyz' | wc -l)" -eq 1
shasum -a 256 benzene_dimer.xyz \
  benzene_dimer_samples_run01/sampling_manifest.json \
  benzene_dimer_samples_run01/monomers/*.xyz \
  benzene_dimer_samples_run01/dimers/*.xyz
```

Open both monomers and the dimer, confirm the two six-membered rings, and check
that no intermolecular contact was inferred as a covalent bond. The manifest
must list source indices `0..11` and `12..23` as two disjoint components.

The output XYZ files contain geometry only. Assign and validate each downstream
job's charge and multiplicity independently.
