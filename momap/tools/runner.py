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
import hashlib
import shlex
import stat
import subprocess
import argparse
import shutil
import tempfile
from pathlib import Path

MOMAP_ROOT = os.environ.get('MOMAP_ROOT')

PATCHED_SCRIPT = None
PATCHED_LAUNCHER_IDENTITY = None
_SLURM_NAME = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$')
_SHA256 = re.compile(r'^[0-9a-fA-F]{64}$')
SUPPORTED_MOMAP_BUILD = '2024A'


def _copy_file_exclusive(source, target):
    """Copy one source to a fresh regular target without following a target link."""
    source = Path(source).expanduser().resolve(strict=True)
    target = Path(target)
    if target.exists() or target.is_symlink():
        raise FileExistsError(f"Refusing to overwrite staged file: {target}")
    created = False
    try:
        with source.open('rb') as source_handle:
            target_handle = target.open('xb')
            created = True
            with target_handle:
                shutil.copyfileobj(source_handle, target_handle)
                target_handle.flush()
                os.fsync(target_handle.fileno())
        shutil.copystat(source, target, follow_symlinks=False)
    except Exception:
        if created:
            target.unlink(missing_ok=True)
        raise
    return target

def find_momap_bin():
    """Find momap executable, prefer system PATH."""
    momap = shutil.which('momap')
    if momap:
        return momap
    if MOMAP_ROOT:
        candidate = Path(MOMAP_ROOT).expanduser() / 'bin' / 'momap'
        if candidate.is_file():
            return str(candidate)
    raise RuntimeError(
        "MOMAP launcher not found. Source the licensed installation environment "
        "or set MOMAP_ROOT explicitly."
    )

def _validate_expected_identity(
    expected_build,
    expected_launcher_sha256,
    expected_version_banner=None,
):
    """Validate caller-supplied identity expectations without defaults."""
    if expected_build != SUPPORTED_MOMAP_BUILD:
        raise ValueError(
            f'Only the verified MOMAP {SUPPORTED_MOMAP_BUILD} contract is '
            f'supported; received expected_build={expected_build!r}'
        )
    if not isinstance(expected_launcher_sha256, str) or not _SHA256.fullmatch(
        expected_launcher_sha256
    ):
        raise ValueError('expected_launcher_sha256 must be exactly 64 hex digits')
    digest = expected_launcher_sha256.lower()
    if expected_version_banner is not None:
        if (
            not isinstance(expected_version_banner, str)
            or not expected_version_banner
            or '\n' in expected_version_banner
            or '\r' in expected_version_banner
        ):
            raise ValueError(
                'expected_version_banner must be one non-empty exact full line'
            )
        if expected_build not in expected_version_banner:
            raise ValueError(
                'expected_version_banner must identify the expected MOMAP build'
            )
    return digest


def _stable_read_nofollow(path):
    """Read one regular file without following links or accepting mutation."""
    lexical = Path(os.path.abspath(os.fspath(Path(path).expanduser())))
    before_lstat = os.stat(lexical, follow_symlinks=False)
    if stat.S_ISLNK(before_lstat.st_mode):
        raise OSError(f'Refusing symbolic-link launcher or evidence file: {lexical}')
    if not hasattr(os, 'O_NOFOLLOW'):
        raise RuntimeError('O_NOFOLLOW is required for launcher identity checks')
    fd = os.open(lexical, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode):
            raise OSError(f'Expected a regular file: {lexical}')
        chunks = []
        while True:
            chunk = os.read(fd, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(fd)
    finally:
        os.close(fd)
    after_lstat = os.stat(lexical, follow_symlinks=False)

    def identity(item):
        return (
            item.st_dev,
            item.st_ino,
            item.st_mode,
            item.st_size,
            item.st_mtime_ns,
            item.st_ctime_ns,
        )

    if identity(before_lstat) != identity(before) or identity(before) != identity(after):
        raise OSError(f'File identity changed while reading: {lexical}')
    if identity(after) != identity(after_lstat):
        raise OSError(f'File identity changed after reading: {lexical}')
    data = b''.join(chunks)
    if len(data) != before.st_size:
        raise OSError(f'File size changed while reading: {lexical}')
    return lexical, data, {
        'device': before.st_dev,
        'inode': before.st_ino,
        'size_bytes': before.st_size,
        'mtime_ns': before.st_mtime_ns,
        'ctime_ns': before.st_ctime_ns,
    }


def _file_description(path):
    lexical, data, metadata = _stable_read_nofollow(path)
    return {
        'path': str(lexical),
        'sha256': hashlib.sha256(data).hexdigest(),
        **metadata,
    }


def patch_momap_for_mpi3(
    *, expected_build, expected_launcher_sha256, temp_parent=None
):
    """Create a private hash-bound launcher using ``--hostfile``.

    Slurm callers pass the shared work directory as ``temp_parent`` so compute
    nodes can read the launcher. Local callers use a fresh private system
    temporary directory for each execution.
    """
    global PATCHED_SCRIPT, PATCHED_LAUNCHER_IDENTITY
    expected_launcher_sha256 = _validate_expected_identity(
        expected_build, expected_launcher_sha256
    )
    momap_orig = find_momap_bin()
    original_path, original_bytes, original_stat = _stable_read_nofollow(momap_orig)
    original_sha256 = hashlib.sha256(original_bytes).hexdigest()
    if original_sha256 != expected_launcher_sha256:
        raise ValueError(
            'MOMAP launcher SHA-256 mismatch: '
            f'expected {expected_launcher_sha256}, observed {original_sha256}'
        )
    try:
        content = original_bytes.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise ValueError('MOMAP launcher must be UTF-8 text for MPI patching') from exc
    replacement_count = content.count('-machinefile')
    patched_content = content.replace('-machinefile', '--hostfile')

    parent = None
    if temp_parent is not None:
        parent = Path(temp_parent).expanduser().resolve()
        if not parent.is_dir():
            raise NotADirectoryError(f"Temporary parent does not exist: {parent}")
    tmpdir = Path(tempfile.mkdtemp(prefix='.momap-patched-', dir=parent))
    os.chmod(tmpdir, 0o700)
    patched = tmpdir / 'momap_patched'
    
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, 'O_NOFOLLOW'):
        flags |= os.O_NOFOLLOW
    fd = os.open(patched, flags, 0o700)
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(patched_content)
        f.flush()
        os.fsync(f.fileno())
    os.chmod(patched, 0o700)

    patched_description = _file_description(patched)
    identity = {
        'build': expected_build,
        'original_launcher': {
            'path': str(original_path),
            'sha256': original_sha256,
            **original_stat,
        },
        'patched_launcher': patched_description,
        'mpi_patch': {
            'contract': 'replace_exact_-machinefile_with_--hostfile',
            'replacement_count': replacement_count,
            'original_preserved': True,
        },
    }
    PATCHED_SCRIPT = str(patched)
    PATCHED_LAUNCHER_IDENTITY = identity
    return identity


def _require_exact_banner(log_path, expected_version_banner):
    """Require one exact full-line version banner in a stable run log."""
    _, data, _ = _stable_read_nofollow(log_path)
    text = data.decode('utf-8', errors='replace')
    count = sum(
        line == expected_version_banner for line in text.splitlines()
    )
    if count != 1:
        raise ValueError(
            'MOMAP run log must contain exactly one exact version banner '
            f'{expected_version_banner!r}; found {count}'
        )
    return hashlib.sha256(data).hexdigest()


def validate_expected_identity(
    expected_build, expected_launcher_sha256, expected_version_banner
):
    """Public fail-closed validator for callers that parse 2024A artifacts."""
    return _validate_expected_identity(
        expected_build,
        expected_launcher_sha256,
        expected_version_banner,
    )


def validate_execution_evidence(
    execution,
    workdir,
    started_ns,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Revalidate build, launchers, MPI patch, and exact run-log banner."""
    expected_launcher_sha256 = _validate_expected_identity(
        expected_build,
        expected_launcher_sha256,
        expected_version_banner,
    )
    if not isinstance(execution, dict):
        raise ValueError('MOMAP execution evidence is required for stage acceptance')
    exit_code = execution.get('process_exit_code')
    if isinstance(exit_code, bool) or not isinstance(exit_code, int):
        raise ValueError('MOMAP execution evidence must record an integer exit code')
    identity = execution.get('executable_identity')
    if not isinstance(identity, dict):
        raise ValueError('MOMAP executable identity is required for stage acceptance')
    if identity.get('build') != expected_build:
        raise ValueError('MOMAP executable identity build does not match expectation')
    if identity.get('version_banner') != expected_version_banner:
        raise ValueError('MOMAP executable identity version banner does not match')

    original = identity.get('original_launcher')
    patched = identity.get('patched_launcher')
    patch_contract = identity.get('mpi_patch')
    if not all(isinstance(item, dict) for item in (original, patched, patch_contract)):
        raise ValueError('MOMAP launcher identity and MPI patch contract are required')
    if original.get('sha256') != expected_launcher_sha256:
        raise ValueError('Recorded original MOMAP launcher SHA-256 does not match')
    observed_original = _file_description(original.get('path', ''))
    observed_patched = _file_description(patched.get('path', ''))
    if observed_original['sha256'] != expected_launcher_sha256:
        raise ValueError('Current original MOMAP launcher SHA-256 does not match')
    if observed_patched['sha256'] != patched.get('sha256'):
        raise ValueError('Current patched MOMAP launcher SHA-256 does not match')
    if patch_contract.get('contract') != 'replace_exact_-machinefile_with_--hostfile':
        raise ValueError('Unsupported MOMAP MPI launcher patch contract')
    replacement_count = patch_contract.get('replacement_count')
    if (
        isinstance(replacement_count, bool)
        or not isinstance(replacement_count, int)
        or replacement_count < 0
    ):
        raise ValueError('Invalid MOMAP MPI launcher replacement count')
    _, original_bytes, _ = _stable_read_nofollow(original['path'])
    _, patched_bytes, _ = _stable_read_nofollow(patched['path'])
    if original_bytes.count(b'-machinefile') != replacement_count:
        raise ValueError('MOMAP MPI patch replacement count does not match source')
    if original_bytes.replace(b'-machinefile', b'--hostfile') != patched_bytes:
        raise ValueError('Patched MOMAP launcher does not satisfy recorded contract')
    patched_mode = stat.S_IMODE(os.stat(patched['path'], follow_symlinks=False).st_mode)
    patched_parent_mode = stat.S_IMODE(
        os.stat(Path(patched['path']).parent, follow_symlinks=False).st_mode
    )
    if patched_mode & 0o077 or patched_parent_mode & 0o077:
        raise ValueError('Patched MOMAP launcher and its directory must be private')

    workdir = Path(workdir).expanduser().resolve(strict=True)
    run_log = Path(execution.get('run_log', '')).expanduser()
    run_log_resolved = run_log.resolve(strict=True)
    if run_log_resolved.parent != workdir:
        raise ValueError('MOMAP run log must be inside its stage directory')
    run_log_description = _file_description(run_log)
    if run_log_description['size_bytes'] <= 0:
        raise ValueError('MOMAP run log must be non-empty')
    if run_log_description['mtime_ns'] < started_ns:
        raise ValueError('MOMAP run log is not fresh for this stage')
    if run_log_description['sha256'] != execution.get('run_log_sha256'):
        raise ValueError('MOMAP run log SHA-256 does not match execution evidence')
    _require_exact_banner(run_log, expected_version_banner)
    return {
        'process_exit_code': exit_code,
        'run_log': str(run_log_resolved),
        'run_log_sha256': run_log_description['sha256'],
        'executable_identity': identity,
    }

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
            _copy_file_exclusive(source_fchk, target_fchk)
        return str(target_fchk)

    if not source_chk.exists():
        print(f"⚠️  No .fchk or .chk found for {logpath.name}")
        print(f"   MOMAP requires a formatted checkpoint for {logpath}")
        return None

    if target_fchk.exists() or target_fchk.is_symlink():
        raise FileExistsError(
            f"Refusing to overwrite formatted checkpoint: {target_fchk}"
        )

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
        _copy_file_exclusive(source_log, target_log)

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

def run_momap(
    inputfile,
    workdir=None,
    use_slurm=False,
    slurm_partition=None,
    slurm_nprocs=4,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
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
        return submit_slurm(
            inputfile,
            cwd,
            slurm_partition,
            slurm_nprocs,
            expected_build=expected_build,
            expected_launcher_sha256=expected_launcher_sha256,
            expected_version_banner=expected_version_banner,
        )
    else:
        execution = run_local(
            inputfile,
            cwd,
            expected_build=expected_build,
            expected_launcher_sha256=expected_launcher_sha256,
            expected_version_banner=expected_version_banner,
        )
        return execution['process_exit_code']

def run_local(
    inputfile,
    cwd,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Run MOMAP on the local system."""
    expected_launcher_sha256 = _validate_expected_identity(
        expected_build,
        expected_launcher_sha256,
        expected_version_banner,
    )
    cwd = Path(cwd).expanduser().resolve()
    inputfile = Path(inputfile).expanduser()
    if not inputfile.is_absolute():
        inputfile = cwd / inputfile
    inputfile = inputfile.resolve(strict=True)
    identity = patch_momap_for_mpi3(
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
    )
    patched = identity['patched_launcher']['path']
    create_nodefile(cwd)
    
    env = os.environ.copy()
    if MOMAP_ROOT:
        env.setdefault('MOMAP_ROOT', MOMAP_ROOT)
    # Preserve MOMAP_LICENSE exactly as supplied by the licensed environment.
    
    cmd = [sys.executable, patched, '-i', str(inputfile)]
    print(f"🚀 Running: {' '.join(str(part) for part in cmd)}")

    run_log = cwd / 'momap_runner.log'
    if run_log.exists() or run_log.is_symlink():
        raise FileExistsError(f'Refusing to overwrite MOMAP run log: {run_log}')
    with run_log.open('xb') as log_handle:
        result = subprocess.run(
            cmd,
            cwd=str(cwd),
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
        log_handle.flush()
        os.fsync(log_handle.fileno())
    run_log_sha256 = _require_exact_banner(run_log, expected_version_banner)
    identity = {
        **identity,
        'version_banner': expected_version_banner,
    }
    return {
        'process_exit_code': result.returncode,
        'run_log': str(run_log),
        'run_log_sha256': run_log_sha256,
        'executable_identity': identity,
    }

def submit_slurm(
    inputfile,
    cwd,
    partition,
    nprocs,
    *,
    expected_build,
    expected_launcher_sha256,
    expected_version_banner,
):
    """Submit MOMAP via Slurm to compute nodes."""
    expected_launcher_sha256 = _validate_expected_identity(
        expected_build,
        expected_launcher_sha256,
        expected_version_banner,
    )
    if partition is not None and (
        not isinstance(partition, str) or not partition
    ):
        raise ValueError("Slurm partition must be a non-empty string when set")
    if partition and any(
        not _SLURM_NAME.fullmatch(name) for name in partition.split(',')
    ):
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
    identity = patch_momap_for_mpi3(
        expected_build=expected_build,
        expected_launcher_sha256=expected_launcher_sha256,
        temp_parent=cwd,
    )
    patched = Path(identity['patched_launcher']['path'])
    momap_env_value = os.environ.get('MOMAP_ENV')
    if momap_env_value is None and MOMAP_ROOT:
        candidate = Path(MOMAP_ROOT).expanduser() / 'env.sh'
        if candidate.is_file():
            momap_env_value = str(candidate)
    momap_env = (
        Path(momap_env_value).expanduser().resolve(strict=True)
        if momap_env_value
        else None
    )
    slurm_script = Path(cwd) / 'momap_job.slurm'
    partition_directive = (
        f"#SBATCH --partition={partition}\n" if partition else ""
    )
    environment_setup = (
        f"source {shlex.quote(str(momap_env))}\n" if momap_env else ""
    )

    script = f"""#!/bin/bash
#SBATCH --job-name=momap
{partition_directive}#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --time=12:00:00
#SBATCH --output=momap_%j.log

set -euo pipefail
{environment_setup}cd -- {shlex.quote(str(cwd))}
echo "localhost slots=$SLURM_NTASKS" > nodefile

# Hash-bound private MOMAP launcher. The Slurm job is unsuccessful unless the
# exact expected build banner appears once in its runner log.
{shlex.quote(sys.executable)} {shlex.quote(str(patched))} -i {shlex.quote(str(inputfile))} > momap_runner.log 2>&1
banner_count=$(grep -Fxc -- {shlex.quote(expected_version_banner)} momap_runner.log || true)
test "$banner_count" -eq 1
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
    parser.add_argument(
        '--partition', '-p', default=None,
        help='Optional Slurm partition; omit to use the site default',
    )
    parser.add_argument('--nprocs', '-n', type=int, default=4, help='MPI processes')
    parser.add_argument('--workdir', '-C', default='.', help='Working directory')
    parser.add_argument('--patch-only', action='store_true', help='Only create patched momap script and exit')
    parser.add_argument(
        '--expected-build',
        required=True,
        choices=[SUPPORTED_MOMAP_BUILD],
        help='Verified MOMAP build contract (currently 2024A only)',
    )
    parser.add_argument(
        '--expected-launcher-sha256',
        required=True,
        help='Expected SHA-256 of the original licensed momap launcher',
    )
    parser.add_argument(
        '--expected-version-banner',
        required=True,
        help='Exact full version-banner line expected once in the run log',
    )
    args = parser.parse_args()

    if args.patch_only:
        identity = patch_momap_for_mpi3(
            expected_build=args.expected_build,
            expected_launcher_sha256=args.expected_launcher_sha256,
        )
        patched = identity['patched_launcher']['path']
        print(f"Patched MOMAP: {patched}")
        print(f"Original SHA-256: {identity['original_launcher']['sha256']}")
        print(f"Patched SHA-256: {identity['patched_launcher']['sha256']}")
        print(f"Usage: {sys.executable} {patched} -i momap.inp")
        return 0

    if not args.input:
        parser.print_help()
        return 1

    sys.exit(
        run_momap(
            args.input,
            args.workdir,
            args.slurm,
            args.partition,
            args.nprocs,
            expected_build=args.expected_build,
            expected_launcher_sha256=args.expected_launcher_sha256,
            expected_version_banner=args.expected_version_banner,
        )
    )

if __name__ == '__main__':
    main()
