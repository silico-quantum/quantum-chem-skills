from __future__ import annotations

import json
import hashlib
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "multiwfn" / "scripts" / "run_batch.py"


class MultiwfnBatchTests(unittest.TestCase):
    def make_executable(self, directory: Path, body: str) -> Path:
        executable = directory / "fake_multiwfn.py"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "print('MULTIWFN VERSION 3.8', flush=True)\n"
            + body,
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return executable

    def run_batch(
        self,
        workdir: Path,
        executable: Path,
        timeout: float = 2.0,
        *,
        markers: tuple[str, str, str] = ("LOAD OK", "TASK OK", "CLEAN EXIT"),
        outputs: tuple[str, ...] = ("orbital.cub", "summary.txt"),
        report: str | Path = "batch.json",
        expected_executable_sha256: str | None = None,
        version_marker: str = "MULTIWFN VERSION 3.8",
    ) -> subprocess.CompletedProcess[str]:
        source = workdir / "source.fchk"
        source.write_text("formatted checkpoint fixture\n", encoding="utf-8")
        responses = workdir / "commands.in"
        responses.write_text("1\n2\n0\n", encoding="utf-8")
        command = [
            sys.executable,
            str(RUNNER),
            str(source),
            "--workdir",
            str(workdir),
            "--responses",
            str(responses),
            "--executable",
            str(executable),
            "--expected-executable-sha256",
            expected_executable_sha256
            or hashlib.sha256(executable.read_bytes()).hexdigest(),
            "--timeout-seconds",
            str(timeout),
            "--load-marker",
            markers[0],
            "--version-marker",
            version_marker,
            "--task-marker",
            markers[1],
            "--exit-marker",
            markers[2],
            "--final-response",
            "0",
            "--log",
            "multiwfn.log",
            "--report",
            str(report),
        ]
        for output in outputs:
            command.extend(("--output", output))
        return subprocess.run(
            command,
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )

    def test_accepts_the_exact_piloted_executable_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdin.read()\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )
            expected_hash = hashlib.sha256(executable.read_bytes()).hexdigest()

            result = self.run_batch(workdir, executable)

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))
            self.assertEqual(report["executable"]["expected_sha256"], expected_hash)
            self.assertEqual(report["executable"]["sha256_before"], expected_hash)
            self.assertEqual(report["executable"]["sha256_after"], expected_hash)
            self.assertTrue(report["checks"]["executable_matches_pilot"])
            self.assertEqual(
                report["markers"]["version"], "MULTIWFN VERSION 3.8"
            )
            self.assertTrue(report["checks"]["version_marker"])

    def test_rejects_an_executable_hash_not_saved_by_the_pilot(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "Path('launched').write_text('unexpected launch\\n')\n",
            )

            result = self.run_batch(
                workdir, executable, expected_executable_sha256="0" * 64
            )

            self.assertNotEqual(result.returncode, 0)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))
            self.assertFalse(report["checks"]["executable_matches_pilot"])
            self.assertIsNone(report["executable"]["sha256_after"])
            self.assertFalse((workdir / "launched").exists())
            self.assertFalse((workdir / "multiwfn.log").exists())

    def test_rejects_a_nonunique_exact_version_banner(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdin.read()\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('MULTIWFN VERSION 3.8')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(report["checks"]["version_marker"])
            self.assertEqual(report["status"], "rejected")

    def test_rejects_an_executable_changed_during_execution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdin.read()\n"
                "Path(__file__).write_text('# changed during execution\\n')\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(report["checks"]["executable_unchanged_after_process"])
            self.assertEqual(report["status"], "rejected")

    @unittest.skipUnless(os.name == "posix" and hasattr(os, "fork"), "POSIX only")
    def test_rejects_and_kills_background_descendant_before_artifact_capture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import os\n"
                "import sys\n"
                "import time\n"
                "sys.stdin.read()\n"
                "Path('orbital.cub').write_text('initial cube\\n')\n"
                "Path('summary.txt').write_text('initial summary\\n')\n"
                "print('LOAD OK', flush=True)\n"
                "print('TASK OK', flush=True)\n"
                "print('CLEAN EXIT', flush=True)\n"
                "if os.fork() == 0:\n"
                "    time.sleep(0.5)\n"
                "    with Path('orbital.cub').open('a') as handle:\n"
                "        handle.write('late descendant write\\n')\n"
                "    os._exit(0)\n"
                "raise SystemExit(0)\n",
            )

            result = self.run_batch(workdir, executable)
            report_path = workdir / "batch.json"
            artifact_before = (workdir / "orbital.cub").read_bytes()
            report_before = report_path.read_bytes()
            time.sleep(0.7)

            self.assertNotEqual(result.returncode, 0)
            report = json.loads(report_before)
            self.assertIn("checks", report, report)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["process_group_had_no_survivors"])
            self.assertTrue(report["checks"]["process_group_cleanup_confirmed"])
            self.assertEqual((workdir / "orbital.cub").read_bytes(), artifact_before)
            self.assertEqual(report_path.read_bytes(), report_before)

    @unittest.skipUnless(os.name == "posix" and hasattr(os, "fork"), "POSIX only")
    def test_rejects_background_descendant_that_closes_inherited_sentinel(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import os\n"
                "import sys\n"
                "import time\n"
                "sys.stdin.read()\n"
                "Path('orbital.cub').write_text('initial cube\\n')\n"
                "Path('summary.txt').write_text('initial summary\\n')\n"
                "print('LOAD OK', flush=True)\n"
                "print('TASK OK', flush=True)\n"
                "print('CLEAN EXIT', flush=True)\n"
                "os.closerange(3, 1024)\n"
                "if os.fork() == 0:\n"
                "    time.sleep(0.5)\n"
                "    with Path('orbital.cub').open('a') as handle:\n"
                "        handle.write('late descendant write\\n')\n"
                "    os._exit(0)\n"
                "raise SystemExit(0)\n",
            )

            result = self.run_batch(workdir, executable)
            report_path = workdir / "batch.json"
            artifact_before = (workdir / "orbital.cub").read_bytes()
            report_before = report_path.read_bytes()
            time.sleep(0.7)

            self.assertNotEqual(result.returncode, 0)
            report = json.loads(report_before)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["process_group_had_no_survivors"])
            self.assertTrue(report["checks"]["process_group_cleanup_confirmed"])
            self.assertEqual((workdir / "orbital.cub").read_bytes(), artifact_before)
            self.assertEqual(report_path.read_bytes(), report_before)

    def test_report_path_must_remain_inside_workdir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "work"
            workdir.mkdir()
            external_report = root / "escaped.json"
            executable = self.make_executable(
                workdir,
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(
                workdir,
                executable,
                outputs=("orbital.cub",),
                report=external_report,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(external_report.exists())

    def test_child_reads_captured_input_and_response_snapshots(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            original_input = b"formatted checkpoint fixture\n"
            original_responses = b"1\n2\n0\n"
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "Path('source.fchk').write_bytes(b'replaced input\\n')\n"
                "Path('commands.in').write_bytes(b'replaced responses\\n')\n"
                "actual_input = Path(sys.argv[1]).read_bytes()\n"
                "actual_responses = sys.stdin.buffer.read()\n"
                "Path('orbital.cub').write_bytes(actual_input)\n"
                "Path('summary.txt').write_bytes(actual_responses)\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual((workdir / "orbital.cub").read_bytes(), original_input)
            self.assertEqual((workdir / "summary.txt").read_bytes(), original_responses)
            self.assertEqual(report["input_sha256"], report["input_snapshot_sha256"])
            self.assertEqual(
                report["responses_sha256"], report["responses_snapshot_sha256"]
            )
            self.assertNotEqual(report["command"][1], str(workdir / "source.fchk"))

    def test_child_cannot_mutate_snapshots_and_still_be_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "input_snapshot = Path(sys.argv[1])\n"
                "responses_snapshot = input_snapshot.parent / 'responses.in'\n"
                "input_snapshot.chmod(0o600)\n"
                "responses_snapshot.chmod(0o600)\n"
                "input_snapshot.write_bytes(b'MUTATED INPUT\\n')\n"
                "responses_snapshot.write_bytes(b'MUTATED RESPONSES\\n')\n"
                "sys.stdin.buffer.read()\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["input_snapshot_unchanged"])
            self.assertFalse(report["checks"]["responses_snapshot_unchanged"])

    def test_accepts_all_markers_and_all_fresh_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdin.read()\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(report["status"], "process_validated")
            self.assertEqual(report["process_status"], "accepted")
            self.assertEqual(
                report["scientific_status"], "pending_independent_validation"
            )
            self.assertEqual(len(report["outputs"]), 2)

    def test_symlink_analysis_output_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "work"
            workdir.mkdir()
            external = root / "external.cub"
            external.write_text("unrelated external artifact\n", encoding="utf-8")
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "import sys\n"
                "sys.stdin.read()\n"
                f"Path('orbital.cub').symlink_to({str(external)!r})\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["all_expected_outputs_non_empty"])
            self.assertEqual(
                external.read_text(encoding="utf-8"), "unrelated external artifact\n"
            )

    def test_missing_clean_exit_marker_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD OK')\nprint('TASK OK')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["clean_exit_marker"])

    def test_timeout_terminates_and_rejects_the_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "import time\ntime.sleep(5)\n",
            )

            result = self.run_batch(workdir, executable, timeout=0.1)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(report["status"], "rejected")
            self.assertFalse(report["checks"]["completed_before_timeout"])

    def test_log_cannot_be_declared_as_an_analysis_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "print('LOAD OK')\nprint('TASK OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(
                workdir, executable, outputs=("multiwfn.log",)
            )
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(report["status"], "rejected")
            self.assertIn("distinct", " ".join(report["errors"]))

    def test_markers_must_be_nonoverlapping_and_ordered(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('TASK OK')\nprint('LOAD OK')\nprint('CLEAN EXIT')\n",
            )

            result = self.run_batch(workdir, executable)
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(report["checks"]["markers_unique_and_ordered"])

        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            executable = self.make_executable(
                workdir,
                "from pathlib import Path\n"
                "Path('orbital.cub').write_text('cube\\n')\n"
                "Path('summary.txt').write_text('summary\\n')\n"
                "print('LOAD TASK CLEAN')\n",
            )

            result = self.run_batch(
                workdir,
                executable,
                markers=("LOAD", "LOAD TASK", "LOAD TASK CLEAN"),
            )
            report = json.loads((workdir / "batch.json").read_text(encoding="utf-8"))

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("overlap", " ".join(report["errors"]))


if __name__ == "__main__":
    unittest.main()
