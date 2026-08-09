# PySCF Version Notes

## Release Snapshot

As checked on 2026-08-09, the latest stable release published on PyPI is
**PySCF 2.14.0**, uploaded on 2026-07-18. Treat this as a dated snapshot rather
than a permanently current version claim.

- [PySCF on PyPI](https://pypi.org/project/pyscf/)
- [Official repository](https://github.com/pyscf/pyscf)
- [Official releases](https://github.com/pyscf/pyscf/releases)
- [Official documentation](https://pyscf.org/)
- [Upstream changelog](https://github.com/pyscf/pyscf/blob/master/CHANGELOG)

## Install or Update

```bash
# Latest stable release available to the selected package index
python -m pip install --upgrade pyscf

# Optional package containing additional and newer features
python -m pip install pyscf-forge

# Install the extension set exposed by the current PyPI release
python -m pip install 'pyscf[all]'
```

Verify the environment actually selected by the shell:

```bash
python -c "import platform, pyscf; print(platform.python_version(), pyscf.__version__)"
```

## Reproducible Environments

For a reproducible workflow, pin the tested release rather than relying on the
moving latest version:

```bash
python -m pip install 'pyscf==2.14.0'
```

Record at least the PySCF version, Python version, operating system and
architecture, numerical backend, thread settings, optional extensions, and
the exact calculation input. Repository examples may describe an older
environment and are not automatically evidence of compatibility with a newer
release.

## Compatibility

The PyPI metadata for 2.14.0 requires Python 3.7 or newer and lists Python
classifiers through 3.14. Wheel availability is platform- and
architecture-dependent. Check the selected release's metadata and the
upstream installation guide before choosing Python, NumPy, SciPy, JAX, GPU, or
extension versions for a production calculation.
