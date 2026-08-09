#!/usr/bin/env python3
"""
momap-run — One-command MOMAP wrapper with auto MPI fix + formchk + Slurm.

Handles:
  - OpenMPI 3.x --hostfile vs -machinefile patching (auto)
  - formchk generation for missing .fchk files (auto)
  - Local execution vs Slurm submission
  - Sub-commands: evc, spec, isc, full
"""
import sys
import os
import re
import shlex
import subprocess
import argparse
import shutil
import tempfile
from pathlib import Path

MOMAP_ROOT = os.environ.get('MOMAP_ROOT', '/opt/MOMAP-2024A')
MOMAP_BIN = os.path.join(MOMAP_ROOT, 'bin')
MOMAP_MPI_BIN = os.path.join(MOMAP_ROOT, 'bin', 'openmpi', 'bin')

PATCHED_SCRIPT = None
_SLURM_NAME = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$')

def find_momap_bin():
    """Find momap executable, prefer system PATH."""
    momap = shutil.which('momap')
    if momap:
        return momap
    momap = os.path.join(MOMAP_BIN, 'momap')
    if os.path.exists(momap):
        return momap
    raise RuntimeError("MOMAP not found. source /opt/MOMAP-2024A/env.sh first")

def patch_momap_for_mpi3(temp_parent=None):
    """Create a private patched launcher using ``--hostfile``.

    Slurm callers pass the shared work directory as ``temp_parent`` so compute
    nodes can read the launcher. Local callers use the system temporary
    directory and reuse the process-local cached path.
    """
    global PATCHED_SCRIPT
    if temp_parent is None and PATCHED_SCRIPT and os.path.exists(PATCHED_SCRIPT):
        return PATCHED_SCRIPT
    
    momap_orig = find_momap_bin()
    parent = None
    if temp_parent is not None:
        parent = Path(temp_parent).expanduser().resolve()
        if not parent.is_dir():
            raise NotADirectoryError(f"Temporary parent does not exist: {parent}")
    tmpdir = Path(tempfile.mkdtemp(prefix='.momap-patched-', dir=parent))
    os.chmod(tmpdir, 0o700)
    patched = tmpdir / 'momap_patched'
    
    with open(momap_orig) as f:
        content = f.read()
    
    content = content.replace('-machinefile', '--hostfile')
    
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, 'O_NOFOLLOW'):
        flags |= os.O_NOFOLLOW
    fd = os.open(patched, flags, 0o700)
    with os.fdopen(fd, 'w') as f:
        f.write(content)
    os.chmod(patched, 0o700)
    
    if temp_parent is None:
        PATCHED_SCRIPT = str(patched)
    return str(patched)

def ensure_fchk(logpath, workdir=None, target_stem=None):
    """Copy or generate the formatted checkpoint in ``workdir``.

    When ``workdir`` is omitted, this preserves the historical behavior of
    placing the file next to the Gaussian log. All paths passed to ``formchk``
    are absolute so changing its working directory cannot duplicate a relative
    path prefix.
    """
    logpath = Path(logpath).expanduser().resolve(strict=True)
    source_fchk = logpath.with_suffix('.fchk')
    source_chk = logpath.with_suffix('.chk')
    target_dir = Path(workdir).expanduser().resolve() if workdir else logpath.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    target_fchk = target_dir / f"{target_stem or logpath.stem}.fchk"

    if source_fchk.exists():
        if source_fchk.resolve() != target_fchk.resolve():
            shutil.copy2(source_fchk, target_fchk)
        return str(target_fchk)

    if not source_chk.exists():
        print(f"⚠️  No .fchk or .chk found for {logpath.name}")
        print(f"   MOMAP requires a formatted checkpoint for {logpath}")
        return None

    print(f"🔄 Generating {target_fchk.name} from {source_chk.name} ...")
    result = subprocess.run(
        ['formchk', str(source_chk), str(target_fchk)],
        cwd=str(target_dir), capture_output=True, text=True
    )
    if result.returncode == 0 and target_fchk.is_file():
        print(f"   ✅ {target_fchk.name}")
        return str(target_fchk)
    else:
        detail = result.stderr or f"formchk did not create {target_fchk}"
        print(f"   ❌ formchk failed: {detail}")
        return None


def stage_gaussian_inputs(logpath, workdir, target_stem=None):
    """Stage a Gaussian log and its matching fchk into a MOMAP workdir."""
    source_log = Path(logpath).expanduser().resolve(strict=True)
    target_dir = Path(workdir).expanduser().resolve()
    target_dir.mkdir(parents=True, exist_ok=True)
    target_log = target_dir / f"{target_stem or source_log.stem}{source_log.suffix}"

    if source_log.resolve() != target_log.resolve():
        shutil.copy2(source_log, target_log)

    staged_fchk = ensure_fchk(source_log, target_dir, target_stem=target_stem)
    if staged_fchk is None:
        raise FileNotFoundError(
            f"Cannot stage formatted checkpoint for Gaussian log {source_log}"
        )
    return target_log, Path(staged_fchk)

def create_nodefile(workdir, slots=4):
    """Create OpenMPI 3.x compatible hostfile."""
    nodefile = Path(workdir) / 'nodefile'
    with open(nodefile, 'w') as f:
        f.write(f"localhost slots={slots}\n")
    return str(nodefile)

def run_momap(inputfile, workdir=None, use_slurm=False, slurm_partition='X32Cv4', slurm_nprocs=4):
    """Run MOMAP calculation."""
    inputfile = Path(inputfile).expanduser().resolve(strict=True)
    cwd = (
        Path(workdir).expanduser().resolve()
        if workdir is not None
        else inputfile.parent
    )
    if not cwd.is_dir():
        raise NotADirectoryError(f"MOMAP workdir does not exist: {cwd}")
    
    if use_slurm:
        return submit_slurm(inputfile, cwd, slurm_partition, slurm_nprocs)
    else:
        return run_local(inputfile, cwd)

def run_local(inputfile, cwd):
    """Run MOMAP on the local system."""
    cwd = Path(cwd).expanduser().resolve()
    inputfile = Path(inputfile).expanduser()
    if not inputfile.is_absolute():
        inputfile = cwd / inputfile
    inputfile = inputfile.resolve(strict=True)
    patched = patch_momap_for_mpi3()
    create_nodefile(cwd)
    
    env = os.environ.copy()
    env['MOMAP_ROOT'] = MOMAP_ROOT
    env['MOMAP_LICENSE'] = os.path.join(MOMAP_ROOT, 'license', 'hzwtech.lic')
    
    cmd = [sys.executable, patched, '-i', str(inputfile)]
    print(f"🚀 Running: {' '.join(str(part) for part in cmd)}")
    
    result = subprocess.run(cmd, cwd=str(cwd), env=env, capture_output=False)
    return result.returncode

def submit_slurm(inputfile, cwd, partition, nprocs):
    """Submit MOMAP via Slurm to compute nodes."""
    if not isinstance(partition, str) or not partition:
        raise ValueError("Slurm partition must be a non-empty string")
    if any(not _SLURM_NAME.fullmatch(name) for name in partition.split(',')):
        raise ValueError(f"Invalid Slurm partition: {partition!r}")
    if isinstance(nprocs, bool) or not isinstance(nprocs, int) or not 1 <= nprocs <= 4096:
        raise ValueError("Slurm process count must be an integer from 1 to 4096")

    cwd = Path(cwd).expanduser().resolve()
    if not cwd.is_dir():
        raise NotADirectoryError(f"MOMAP workdir does not exist: {cwd}")
    inputfile = Path(inputfile).expanduser()
    if not inputfile.is_absolute():
        inputfile = cwd / inputfile
    inputfile = inputfile.resolve(strict=True)
    patched = Path(patch_momap_for_mpi3(cwd)).expanduser().resolve(strict=True)
    momap_env = Path(
        os.environ.get('MOMAP_ENV', os.path.join(MOMAP_ROOT, 'env.sh'))
    ).expanduser().resolve()
    slurm_script = Path(cwd) / 'momap_job.slurm'
    
    script = f"""#!/bin/bash
#SBATCH --job-name=momap
#SBATCH --partition={partition}
#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --time=12:00:00
#SBATCH --output=momap_%j.log

source {shlex.quote(str(momap_env))}
cd -- {shlex.quote(str(cwd))}
echo "localhost slots=$SLURM_NTASKS" > nodefile

# Patched momap
{shlex.quote(sys.executable)} {shlex.quote(str(patched))} -i {shlex.quote(str(inputfile))}
"""
    with open(slurm_script, 'w') as f:
        f.write(script)
    
    print(f"📤 Submitting Slurm job...")
    result = subprocess.run(['sbatch', str(slurm_script)], capture_output=True, text=True)
    print(result.stdout.strip())
    return result.returncode

def generate_evc_input(s0_log, s1_log, output='momap_evc.inp'):
    """Generate EVC input file."""
    content = f"""do_evc = 1

&evc
  ffreq(1) = "{s0_log}"
  ffreq(2) = "{s1_log}"
  sort_mode = 1
/
"""
    with open(output, 'w') as f:
        f.write(content)
    return output

def cmd_extract(args):
    """Re-export extract functionality."""
    from extract import main as extract_main
    sys.argv = ['extract'] + args.passthrough if hasattr(args, 'passthrough') else []
    # Delegate to extract module - just use subprocess
    cmd = [sys.executable, os.path.join(os.path.dirname(__file__), 'extract.py')] + args.passthrough
    os.execv(cmd[0], cmd)

def main():
    parser = argparse.ArgumentParser(description='MOMAP runner — one-command wrapper')
    parser.add_argument('input', nargs='?', help='MOMAP input file (momap.inp)')
    parser.add_argument('--slurm', '-s', action='store_true', help='Submit via Slurm')
    parser.add_argument('--partition', '-p', default='X32Cv4', help='Slurm partition')
    parser.add_argument('--nprocs', '-n', type=int, default=4, help='MPI processes')
    parser.add_argument('--workdir', '-C', default='.', help='Working directory')
    parser.add_argument('--patch-only', action='store_true', help='Only create patched momap script and exit')
    args = parser.parse_args()

    if args.patch_only:
        patched = patch_momap_for_mpi3()
        print(f"Patched MOMAP: {patched}")
        print(f"Usage: python {patched} -i momap.inp")
        return 0

    if not args.input:
        parser.print_help()
        return 1

    sys.exit(run_momap(args.input, args.workdir, args.slurm, args.partition, args.nprocs))

if __name__ == '__main__':
    main()
