from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "xtb-cluster-md" / "scripts" / "build_cluster.py"
VALIDATOR = ROOT / "xtb-cluster-md" / "scripts" / "validate_md_run.py"
RUNNER = ROOT / "xtb-cluster-md" / "scripts" / "run_xtb_md.py"
TRAJECTORY_IO = ROOT / "xtb-cluster-md" / "scripts" / "trajectory_io.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("xtb_build_cluster", BUILDER)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_trajectory_io():
    spec = importlib.util.spec_from_file_location("xtb_trajectory_io", TRAJECTORY_IO)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sdf_record(title: str = "H2") -> str:
    return (
        f"{title}\n"
        "  contract-test       3D\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.7400    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n"
        "$$$$\n"
    )


def xyz_frame(comment: str, shift: float) -> str:
    return (
        "2\n"
        f"{comment}\n"
        f"He {shift:.3f} 0.0 0.0\n"
        f"He {2.0 + shift:.3f} 0.0 0.0\n"
    )


class XtbContractTests(unittest.TestCase):
    def test_gif_publication_never_overwrites_a_racing_output(self) -> None:
        trajectory_io = load_trajectory_io()

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "animation.gif"
            manifest = root / "animation.json"
            sentinel = b"concurrent user output\n"

            class OpenedGif:
                format = "GIF"
                n_frames = 1
                size = (2, 2)

                def __enter__(self):
                    return self

                def __exit__(self, *_args):
                    return False

            fake_pil = types.ModuleType("PIL")
            fake_pil.Image = types.SimpleNamespace(open=lambda _path: OpenedGif())

            class RacingFrame:
                def save(self, path, **_kwargs):
                    output.write_bytes(sentinel)
                    Path(path).write_bytes(b"GIF89a fixture")

            with mock.patch.dict(sys.modules, {"PIL": fake_pil}):
                with self.assertRaises(FileExistsError):
                    trajectory_io.save_verified_gif(
                        [RacingFrame()], output, manifest, 100, {}, {}
                    )

            self.assertEqual(output.read_bytes(), sentinel)
            self.assertFalse(manifest.exists())
            self.assertFalse((root / ".animation.gif.partial").exists())
            self.assertFalse((root / ".animation.json.partial").exists())

    def test_gif_publication_rolls_back_if_manifest_races(self) -> None:
        trajectory_io = load_trajectory_io()

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / "animation.gif"
            manifest = root / "animation.json"
            sentinel = b"concurrent user manifest\n"

            class OpenedGif:
                format = "GIF"
                n_frames = 1
                size = (2, 2)

                def __enter__(self):
                    return self

                def __exit__(self, *_args):
                    return False

            fake_pil = types.ModuleType("PIL")
            fake_pil.Image = types.SimpleNamespace(open=lambda _path: OpenedGif())

            class RacingFrame:
                def save(self, path, **_kwargs):
                    manifest.write_bytes(sentinel)
                    Path(path).write_bytes(b"GIF89a fixture")

            with mock.patch.dict(sys.modules, {"PIL": fake_pil}):
                with self.assertRaises(FileExistsError):
                    trajectory_io.save_verified_gif(
                        [RacingFrame()], output, manifest, 100, {}, {}
                    )

            self.assertFalse(output.exists())
            self.assertEqual(manifest.read_bytes(), sentinel)
            self.assertFalse((root / ".animation.gif.partial").exists())
            self.assertFalse((root / ".animation.json.partial").exists())

    def test_builder_rejects_multiple_sdf_records(self) -> None:
        builder = load_builder()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "two-records.sdf"
            path.write_text(sdf_record("first") + sdf_record("second"))

            with self.assertRaisesRegex(ValueError, "exactly one"):
                builder.read_sdf_coords(path)

    def test_builder_emits_connected_manifest_with_explicit_topology(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sdf = root / "h2.sdf"
            sdf.write_text(sdf_record())
            output = root / "h2_n3.xyz"
            manifest_path = root / "build.json"

            result = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
                    "--sdf",
                    str(sdf),
                    "--monomer-id",
                    "H2 test fixture",
                    "--coordinate-provenance",
                    "synthetic 3D unit-test fixture",
                    "--molecules",
                    "3",
                    "--box",
                    "5.0",
                    "--min-distance",
                    "1.0",
                    "--max-neighbor-distance",
                    "4.0",
                    "--seed",
                    "42",
                    "--out",
                    str(output),
                    "--manifest",
                    str(manifest_path),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(manifest["status"], "accepted")
            self.assertEqual(manifest["source"]["atom_labels"], ["H", "H"])
            self.assertEqual(manifest["source"]["bonds_zero_based"], [[0, 1, 1]])
            self.assertLessEqual(
                manifest["sampling"][
                    "largest_nearest_neighbor_distance_angstrom"
                ],
                4.0,
            )
            self.assertEqual(manifest["output"]["sha256"], sha256(output))
            self.assertFalse(
                Path(manifest["publication"]["incomplete_marker"]).exists()
            )
            self.assertFalse((root / ".h2_n3.xyz.partial").exists())

    def test_builder_hashes_the_exact_captured_sdf_bytes_it_parses(self) -> None:
        builder = load_builder()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sdf = root / "source.sdf"
            initial = sdf_record("captured")
            sdf.write_text(initial)
            output = root / "cluster.xyz"
            manifest_path = root / "build.json"
            original_parser = builder.parse_sdf_text

            def mutate_after_capture(text: str):
                record = original_parser(text)
                sdf.write_text(sdf_record("changed-after-capture"))
                return record

            argv = [
                str(BUILDER),
                "--sdf",
                str(sdf),
                "--monomer-id",
                "captured fixture",
                "--coordinate-provenance",
                "synthetic unit-test fixture",
                "--molecules",
                "1",
                "--box",
                "5.0",
                "--min-distance",
                "1.0",
                "--max-neighbor-distance",
                "4.0",
                "--seed",
                "42",
                "--out",
                str(output),
                "--manifest",
                str(manifest_path),
            ]
            with mock.patch.object(builder, "parse_sdf_text", mutate_after_capture):
                with mock.patch.object(sys, "argv", argv):
                    self.assertEqual(builder.main(), 0)

            manifest = json.loads(manifest_path.read_text())
            self.assertEqual(
                manifest["source"]["sha256"],
                hashlib.sha256(initial.encode("utf-8")).hexdigest(),
            )
            self.assertEqual(manifest["source"]["title"], "captured")
            self.assertNotEqual(manifest["source"]["sha256"], sha256(sdf))

    def test_runner_preserves_nonzero_exit_and_binds_log(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "md"
            run_dir.mkdir()
            fake_xtb = run_dir / "fake_xtb.py"
            fake_xtb.write_text(
                "#!/usr/bin/env python3\n"
                "import sys\n"
                "print('synthetic xTB failure')\n"
                "raise SystemExit(7)\n"
            )
            fake_xtb.chmod(0o755)
            input_file = run_dir / "input.xyz"
            input_file.write_text("1\nfixture\nHe 0 0 0\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    str(fake_xtb),
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 7, result.stdout + result.stderr)
            receipt = json.loads((run_dir / "run_receipt.json").read_text())
            self.assertEqual(receipt["status"], "completed")
            self.assertEqual(receipt["returncode"], 7)
            self.assertEqual(receipt["log"]["sha256"], sha256(run_dir / "md.log"))

    def test_runner_hashes_the_executable_before_launch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            fake_xtb = run_dir / "mutable_xtb.py"
            original = (
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "Path(__file__).write_text('# changed after launch\\n')\n"
                "raise SystemExit(7)\n"
            )
            fake_xtb.write_text(original)
            fake_xtb.chmod(0o755)
            input_file = run_dir / "input.xyz"
            input_file.write_text("1\nfixture\nHe 0 0 0\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    str(fake_xtb),
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 7, result.stdout + result.stderr)
            receipt = json.loads((run_dir / "run_receipt.json").read_text())
            self.assertEqual(
                receipt["executable"]["sha256"],
                hashlib.sha256(original.encode("utf-8")).hexdigest(),
            )
            self.assertNotEqual(receipt["executable"]["sha256"], sha256(fake_xtb))
            self.assertEqual(
                receipt["executable"]["hash_timing"],
                "captured before process launch",
            )

    def test_runner_executes_private_snapshots_and_detects_restored_source(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            input_file = run_dir / "input.xyz"
            original = b"1\nfixture\nHe 0 0 0\n"
            input_file.write_bytes(original)
            fake_xtb = run_dir / "mutating_xtb.py"
            fake_xtb.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "import sys\n"
                f"source = Path({str(input_file)!r})\n"
                "captured = source.read_bytes()\n"
                "source.write_bytes(b'changed during execution\\n')\n"
                "source.write_bytes(captured)\n"
                "print('executed-input=' + sys.argv[1])\n"
                "raise SystemExit(7)\n"
            )
            fake_xtb.chmod(0o755)

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    str(fake_xtb),
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 7, result.stdout + result.stderr)
            self.assertEqual(input_file.read_bytes(), original)
            receipt = json.loads((run_dir / "run_receipt.json").read_text())
            self.assertEqual(receipt["schema_version"], 3)
            record = receipt["input_files"][0]
            snapshot = Path(record["snapshot"]["path"])
            self.assertEqual(snapshot.read_bytes(), original)
            self.assertEqual(snapshot.stat().st_mode & 0o777, 0o400)
            self.assertEqual(receipt["argv"].count(record["snapshot"]["argv_token"]), 1)
            self.assertNotIn("input.xyz", receipt["argv"])
            self.assertFalse(record["source_unchanged_after_run"])

    @unittest.skipUnless(os.name == "posix" and hasattr(os, "fork"), "POSIX only")
    def test_runner_rejects_and_kills_background_descendant_before_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            input_file = run_dir / "input.xyz"
            input_file.write_text("1\nfixture\nHe 0 0 0\n")
            fake_xtb = run_dir / "background_xtb.py"
            fake_xtb.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "import os\n"
                "import time\n"
                "Path('xtb.trj').write_text('1\\ninitial\\nHe 0 0 0\\n')\n"
                "Path('xtbmdok').write_text('ok\\n')\n"
                "print('synthetic xTB log', flush=True)\n"
                "if os.fork() == 0:\n"
                "    time.sleep(0.5)\n"
                "    with Path('xtb.trj').open('a') as handle:\n"
                "        handle.write('late descendant write\\n')\n"
                "    os._exit(0)\n"
                "raise SystemExit(0)\n"
            )
            fake_xtb.chmod(0o755)

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    str(fake_xtb),
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )
            receipt_path = run_dir / "run_receipt.json"
            trajectory_before = (run_dir / "xtb.trj").read_bytes()
            receipt_before = receipt_path.read_bytes()
            time.sleep(0.7)

            self.assertNotEqual(result.returncode, 0)
            receipt = json.loads(receipt_before)
            self.assertEqual(receipt["status"], "rejected_descendants")
            self.assertFalse(receipt["process_group"]["had_no_survivors"])
            self.assertTrue(receipt["process_group"]["cleanup_confirmed"])
            self.assertEqual((run_dir / "xtb.trj").read_bytes(), trajectory_before)
            self.assertEqual(receipt_path.read_bytes(), receipt_before)

    @unittest.skipUnless(os.name == "posix" and hasattr(os, "fork"), "POSIX only")
    def test_runner_rejects_descendant_that_closes_inherited_sentinel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            input_file = run_dir / "input.xyz"
            input_file.write_text("1\nfixture\nHe 0 0 0\n")
            fake_xtb = run_dir / "background_xtb.py"
            fake_xtb.write_text(
                "#!/usr/bin/env python3\n"
                "from pathlib import Path\n"
                "import os\n"
                "import time\n"
                "Path('xtb.trj').write_text('1\\ninitial\\nHe 0 0 0\\n')\n"
                "Path('xtbmdok').write_text('ok\\n')\n"
                "print('synthetic xTB log', flush=True)\n"
                "os.closerange(3, 1024)\n"
                "if os.fork() == 0:\n"
                "    time.sleep(0.5)\n"
                "    with Path('xtb.trj').open('a') as handle:\n"
                "        handle.write('late descendant write\\n')\n"
                "    os._exit(0)\n"
                "raise SystemExit(0)\n"
            )
            fake_xtb.chmod(0o755)

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    str(fake_xtb),
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )
            receipt_path = run_dir / "run_receipt.json"
            trajectory_before = (run_dir / "xtb.trj").read_bytes()
            receipt_before = receipt_path.read_bytes()
            time.sleep(0.7)

            self.assertNotEqual(result.returncode, 0)
            receipt = json.loads(receipt_before)
            self.assertEqual(receipt["status"], "rejected_descendants")
            self.assertFalse(receipt["process_group"]["had_no_survivors"])
            self.assertTrue(receipt["process_group"]["cleanup_confirmed"])
            self.assertEqual((run_dir / "xtb.trj").read_bytes(), trajectory_before)
            self.assertEqual(receipt_path.read_bytes(), receipt_before)

    def test_animation_loader_scans_tail_after_render_limit(self) -> None:
        trajectory_io = load_trajectory_io()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            trajectory = root / "xtb.trj"
            trajectory.write_text(
                xyz_frame("frame 0", 0.0)
                + xyz_frame("frame 1", 0.1)
                + "2\nmalformed tail\nHe 0.2 0 0\nNe 2.2 0 0\n"
            )
            validation = root / "validation.json"
            validation.write_text(
                json.dumps(
                    {
                        "status": "accepted",
                        "inputs": {
                            "expected_molecules": 2,
                            "atoms_per_molecule": 1,
                        },
                        "trajectory": {
                            "path": str(trajectory.resolve()),
                            "sha256": sha256(trajectory),
                            "frame_count": 3,
                        },
                        "checks": {"synthetic": True},
                    }
                )
            )

            with self.assertRaisesRegex(ValueError, "ordering changed"):
                trajectory_io.load_accepted_sampled_trajectory(
                    trajectory, validation, 2, 1, stride=1, max_frames=1
                )

    def make_md_run(
        self, root: Path, *, exit_code: int = 0, include_md_exit: bool = True
    ) -> list[str]:
        run_dir = root / "md"
        run_dir.mkdir()
        input_xyz = run_dir / "input.xyz"
        input_xyz.write_text(xyz_frame("input", 0.0))
        md_input = run_dir / "md.inp"
        md_input.write_text("$md\n time=0.1\n$end\n")
        snapshot_dir = run_dir / ".xtb-input-snapshots"
        snapshot_dir.mkdir(mode=0o700)
        input_snapshot = snapshot_dir / "000-input.xyz"
        input_snapshot.write_bytes(input_xyz.read_bytes())
        input_snapshot.chmod(0o400)
        md_snapshot = snapshot_dir / "001-md.inp"
        md_snapshot.write_bytes(md_input.read_bytes())
        md_snapshot.chmod(0o400)
        trajectory = run_dir / "xtb.trj"
        trajectory.write_text(
            xyz_frame("frame 0", 0.0)
            + xyz_frame("frame 1", 0.1)
            + xyz_frame("frame 2", 0.2)
        )
        (run_dir / "xtbmdok").write_text("ok\n")
        marker = run_dir / "xtbmdok"
        md_exit = "normal exit of md()\n" if include_md_exit else ""
        fake_executable = run_dir / "xtb"
        fake_executable.write_text("#!/bin/sh\nexit 0\n")
        fake_executable.chmod(0o755)
        log = run_dir / "md.log"
        log.write_text(
            "xtb version 6.7.1\n"
            "Calculation Setup\n"
            "program call : xtb .xtb-input-snapshots/000-input.xyz --gfnff "
            "--chrg 0 --uhf 0 --md -I .xtb-input-snapshots/001-md.inp\n"
            "Molecular Dynamics\n"
            "MD time /ps : 0.10\n"
            "dt /fs : 1.00\n"
            "SCC accuracy : 2.00\n"
            "temperature /K : 300.00\n"
            "max steps : 100\n"
            "dumpstep(trj) /fs : 50.00 50\n"
            "H atoms mass (amu) : 4\n"
            "SHAKE on. # bonds : 0 all: F\n"
            "Berendsen THERMOSTAT on\n"
            + md_exit
            + "normal termination of xtb\n"
        )
        started_time_ns = min(
            path.stat().st_mtime_ns for path in (trajectory, marker, log)
        ) - 1
        receipt = run_dir / "run_receipt.json"
        receipt.write_text(
            json.dumps(
                {
                    "schema_version": 3,
                    "status": "completed",
                    "working_directory": str(run_dir.resolve()),
                    "argv": [
                        str(fake_executable.resolve()),
                        ".xtb-input-snapshots/000-input.xyz",
                        "--gfnff",
                        "--chrg",
                        "0",
                        "--uhf",
                        "0",
                        "--md",
                        "-I",
                        ".xtb-input-snapshots/001-md.inp",
                    ],
                    "returncode": exit_code,
                    "signal": None,
                    "started_time_ns": started_time_ns,
                    "executable": {
                        "requested": "xtb",
                        "resolved_path": str(fake_executable.resolve()),
                        "sha256": sha256(fake_executable),
                        "size_bytes": fake_executable.stat().st_size,
                        "mtime_ns": fake_executable.stat().st_mtime_ns,
                        "hash_timing": "captured before process launch",
                    },
                    "input_files": [
                        {
                            "source": {
                                "path": str(input_xyz.resolve()),
                                "sha256": sha256(input_xyz),
                                "size_bytes": input_xyz.stat().st_size,
                                "device": input_xyz.stat().st_dev,
                                "inode": input_xyz.stat().st_ino,
                                "mtime_ns": input_xyz.stat().st_mtime_ns,
                                "ctime_ns": input_xyz.stat().st_ctime_ns,
                            },
                            "source_after_run": {
                                "path": str(input_xyz.resolve()),
                                "sha256": sha256(input_xyz),
                                "size_bytes": input_xyz.stat().st_size,
                                "device": input_xyz.stat().st_dev,
                                "inode": input_xyz.stat().st_ino,
                                "mtime_ns": input_xyz.stat().st_mtime_ns,
                                "ctime_ns": input_xyz.stat().st_ctime_ns,
                            },
                            "snapshot": {
                                "path": str(input_snapshot.resolve()),
                                "argv_token": ".xtb-input-snapshots/000-input.xyz",
                                "sha256": sha256(input_snapshot),
                                "size_bytes": input_snapshot.stat().st_size,
                                "mode": "0o400",
                            },
                            "source_unchanged_after_run": True,
                        },
                        {
                            "source": {
                                "path": str(md_input.resolve()),
                                "sha256": sha256(md_input),
                                "size_bytes": md_input.stat().st_size,
                                "device": md_input.stat().st_dev,
                                "inode": md_input.stat().st_ino,
                                "mtime_ns": md_input.stat().st_mtime_ns,
                                "ctime_ns": md_input.stat().st_ctime_ns,
                            },
                            "source_after_run": {
                                "path": str(md_input.resolve()),
                                "sha256": sha256(md_input),
                                "size_bytes": md_input.stat().st_size,
                                "device": md_input.stat().st_dev,
                                "inode": md_input.stat().st_ino,
                                "mtime_ns": md_input.stat().st_mtime_ns,
                                "ctime_ns": md_input.stat().st_ctime_ns,
                            },
                            "snapshot": {
                                "path": str(md_snapshot.resolve()),
                                "argv_token": ".xtb-input-snapshots/001-md.inp",
                                "sha256": sha256(md_snapshot),
                                "size_bytes": md_snapshot.stat().st_size,
                                "mode": "0o400",
                            },
                            "source_unchanged_after_run": True,
                        },
                    ],
                    "log": {
                        "path": str(log.resolve()),
                        "sha256": sha256(log),
                        "size_bytes": log.stat().st_size,
                        "mtime_ns": log.stat().st_mtime_ns,
                    },
                    "generated_artifacts": {
                        "trajectory": {
                            "path": str(trajectory.resolve()),
                            "regular_file": True,
                            "sha256": sha256(trajectory),
                            "size_bytes": trajectory.stat().st_size,
                            "mtime_ns": trajectory.stat().st_mtime_ns,
                        },
                        "success_marker": {
                            "path": str(marker.resolve()),
                            "regular_file": True,
                            "sha256": sha256(marker),
                            "size_bytes": marker.stat().st_size,
                            "mtime_ns": marker.stat().st_mtime_ns,
                        },
                    },
                }
            )
            + "\n"
        )
        report = run_dir / "validation.json"
        return [
            sys.executable,
            str(VALIDATOR),
            "--run-dir",
            str(run_dir),
            "--receipt",
            "run_receipt.json",
            "--log",
            "md.log",
            "--trajectory",
            "xtb.trj",
            "--input-xyz",
            "input.xyz",
            "--md-input",
            "md.inp",
            "--expected-method",
            "gfnff",
            "--expected-charge",
            "0",
            "--expected-uhf",
            "0",
            "--expected-xtb-version",
            "6.7.1",
            "--expected-molecules",
            "2",
            "--atoms-per-molecule",
            "1",
            "--requested-time-ps",
            "0.1",
            "--requested-step-fs",
            "1.0",
            "--requested-dump-fs",
            "50.0",
            "--requested-temperature-k",
            "300.0",
            "--requested-scc-accuracy",
            "2.0",
            "--requested-hydrogen-mass",
            "4.0",
            "--requested-shake-bonds",
            "0",
            "--expect-thermostat",
            "on",
            "--minimum-intermolecular-distance",
            "1.0",
            "--maximum-system-extent",
            "100.0",
            "--report",
            str(report),
        ]

    def test_validator_rejects_nonzero_process_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            command = self.make_md_run(Path(tmp), exit_code=7)
            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((Path(tmp) / "md" / "validation.json").read_text())
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["process_exit_zero"])

    def test_validator_accepts_fully_bound_md_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            command = self.make_md_run(Path(tmp))
            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((Path(tmp) / "md" / "validation.json").read_text())
            self.assertEqual(report["status"], "accepted")
            self.assertTrue(all(report["checks"].values()))
            self.assertEqual(report["trajectory"]["distinct_coordinate_frames"], 3)
            self.assertAlmostEqual(
                report["trajectory"][
                    "global_minimum_inter_distance_angstrom"
                ],
                2.0,
            )

    def test_validator_requires_md_specific_normal_exit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            command = self.make_md_run(Path(tmp), include_md_exit=False)
            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((Path(tmp) / "md" / "validation.json").read_text())
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["md_normal_exit"])

    def test_validator_rejects_an_input_changed_after_the_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            (root / "md" / "md.inp").write_text("$md\n time=999\n$end\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_input_binding"])

    def test_validator_rejects_source_mutated_and_restored_during_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["input_files"][0]["source_unchanged_after_run"] = False
            receipt_path.write_text(json.dumps(receipt) + "\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_input_binding"])

    def test_validator_rejects_source_identity_change_even_if_flag_says_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["input_files"][0]["source_after_run"]["ctime_ns"] += 1
            receipt_path.write_text(json.dumps(receipt) + "\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_input_binding"])

    def test_validator_rejects_symlink_substituted_for_private_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            snapshot_path = Path(receipt["input_files"][0]["snapshot"]["path"])
            replacement = root / "md" / ".xtb-input-snapshots" / "replacement.xyz"
            replacement.write_bytes(snapshot_path.read_bytes())
            replacement.chmod(0o400)
            snapshot_path.unlink()
            snapshot_path.symlink_to(replacement.name)

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_input_binding"])

    def test_validator_rejects_world_traversable_snapshot_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            (root / "md" / ".xtb-input-snapshots").chmod(0o755)

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_input_binding"])

    def test_validator_rejects_method_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            command[command.index("--expected-method") + 1] = "gfn2"

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["argv_method_matches_expected"])

    def test_validator_rejects_charge_or_uhf_mismatch(self) -> None:
        for option, value, check in (
            ("--expected-charge", "1", "argv_charge_matches_expected"),
            ("--expected-uhf", "2", "argv_uhf_matches_expected"),
        ):
            with self.subTest(option=option), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                command = self.make_md_run(root)
                command[command.index(option) + 1] = value

                result = subprocess.run(command, text=True, capture_output=True)

                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                report = json.loads((root / "md" / "validation.json").read_text())
                self.assertFalse(report["checks"][check])

    def test_validator_rejects_xtb_version_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            command[command.index("--expected-xtb-version") + 1] = "6.6.1"

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["log_xtb_version_matches_expected"])

    def test_validator_rejects_program_call_not_equal_to_receipt_argv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            log_path = root / "md" / "md.log"
            changed = log_path.read_text().replace("--chrg 0", "--chrg 1")
            log_path.write_text(changed)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["log"]["sha256"] = sha256(log_path)
            receipt["log"]["size_bytes"] = log_path.stat().st_size
            receipt["log"]["mtime_ns"] = log_path.stat().st_mtime_ns
            receipt_path.write_text(json.dumps(receipt) + "\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["program_call_matches_receipt_argv"])

    def test_validator_rejects_executable_identity_not_bound_to_argv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            receipt["executable"]["resolved_path"] = "/different/xtb"
            receipt_path.write_text(json.dumps(receipt) + "\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(
                report["checks"]["receipt_executable_argv_binding"]
            )

    def test_validator_rejects_a_trajectory_not_bound_as_fresh_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            command = self.make_md_run(root)
            receipt_path = root / "md" / "run_receipt.json"
            receipt = json.loads(receipt_path.read_text())
            trajectory_mtime = receipt["generated_artifacts"]["trajectory"][
                "mtime_ns"
            ]
            receipt["started_time_ns"] = trajectory_mtime + 10_000_000_000
            receipt_path.write_text(json.dumps(receipt) + "\n")

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
            report = json.loads((root / "md" / "validation.json").read_text())
            self.assertFalse(report["checks"]["receipt_trajectory_binding"])

    def test_validator_rejects_nan_physical_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            command = self.make_md_run(Path(tmp))
            index = command.index("--minimum-intermolecular-distance") + 1
            command[index] = "nan"

            result = subprocess.run(command, text=True, capture_output=True)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("physical thresholds", result.stdout + result.stderr)

    def test_runner_rejects_stale_conventional_md_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            input_file = run_dir / "input.xyz"
            input_file.write_text("1\nfixture\nHe 0 0 0\n")
            (run_dir / "xtb.trj").write_text("stale\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER),
                    "--run-dir",
                    str(run_dir),
                    "--input-file",
                    "input.xyz",
                    "--",
                    sys.executable,
                    "input.xyz",
                    "--md",
                ],
                text=True,
                capture_output=True,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("refusing to overwrite", result.stdout + result.stderr)
            self.assertFalse((run_dir / "run_receipt.json").exists())


if __name__ == "__main__":
    unittest.main()
