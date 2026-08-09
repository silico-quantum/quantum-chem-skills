from __future__ import annotations

import importlib.util
import inspect
import json
import math
import os
import shlex
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "momap" / "tools"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MomapExtractionContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.extract = load_module("momap_contract_extract", TOOLS / "extract.py")

    def test_s1_total_energy_uses_last_scf_and_last_state_one_excitation(self) -> None:
        extractor = getattr(self.extract, "extract_s1_total_energy", None)
        self.assertIsNotNone(
            extractor,
            "extract.py must expose a fail-closed S1 total-energy extractor",
        )

        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "s1.log"
            log.write_text(
                textwrap.dedent(
                    """\
                     SCF Done:  E(RB3LYP) =  -100.000000000     A.U.
                     Excited State   1:      Singlet-A      1.5000 eV
                     SCF Done:  E(RB3LYP) =   -99.000000000     A.U.
                     Excited State   1:      Singlet-A      2.0000 eV
                    """
                ),
                encoding="utf-8",
            )

            total = extractor(log)

        self.assertAlmostEqual(total, -99.0 + 2.0 / self.extract.HA2EV, places=10)

    def test_s1_total_energy_fails_when_state_one_excitation_is_missing(self) -> None:
        extractor = getattr(self.extract, "extract_s1_total_energy", None)
        self.assertIsNotNone(extractor)

        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "s1.log"
            log.write_text(
                " SCF Done:  E(RB3LYP) =  -100.000000000     A.U.\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "state 1|S1|excitation"):
                extractor(log)

    def test_transition_endpoints_use_state_one_from_first_and_last_blocks(self) -> None:
        extractor = getattr(
            self.extract, "extract_state1_transition_endpoints", None
        )
        self.assertIsNotNone(
            extractor,
            "extract.py must select state 1 independently in each TD block",
        )

        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "s1.log"
            log.write_text(
                textwrap.dedent(
                    """\
                     Ground to excited state transition electric dipole moments (Au):
                          state          X           Y           Z        Dip. S.      Osc.
                         1        1.0000      2.0000      2.0000      9.0000      0.1000
                         2        9.0000      0.0000      0.0000     81.0000      0.9000

                     Ground to excited state transition electric dipole moments (Au):
                          state          X           Y           Z        Dip. S.      Osc.
                         1        0.0000      3.0000      4.0000     25.0000      0.2000
                         2        8.0000      0.0000      0.0000     64.0000      0.8000

                    """
                ),
                encoding="utf-8",
            )

            absorption, emission = extractor(log)

        self.assertEqual(absorption["state"], 1)
        self.assertEqual(emission["state"], 1)
        self.assertAlmostEqual(absorption["magnitude_au"], 3.0)
        self.assertAlmostEqual(emission["magnitude_au"], 5.0)
        self.assertAlmostEqual(
            absorption["magnitude_debye"], 3.0 * self.extract.AU2DEBYE
        )
        self.assertAlmostEqual(
            emission["magnitude_debye"], 5.0 * self.extract.AU2DEBYE
        )


class MomapRunnerContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.runner = load_module("momap_contract_runner", TOOLS / "runner.py")

    def tearDown(self) -> None:
        self.runner.PATCHED_SCRIPT = None

    def test_gaussian_log_and_fchk_are_staged_into_workdir(self) -> None:
        stage_inputs = getattr(self.runner, "stage_gaussian_inputs", None)
        self.assertIsNotNone(
            stage_inputs,
            "runner.py must stage both Gaussian log and fchk files",
        )

        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            workdir = Path(tmp) / "work"
            source.mkdir()
            log = source / "s1.log"
            fchk = source / "s1.fchk"
            log.write_text("Gaussian log\n", encoding="utf-8")
            fchk.write_text("Formatted checkpoint\n", encoding="utf-8")

            staged_log, staged_fchk = stage_inputs(log, workdir)

            self.assertEqual(staged_log.resolve(), (workdir / "s1.log").resolve())
            self.assertEqual(staged_fchk.resolve(), (workdir / "s1.fchk").resolve())
            self.assertEqual(staged_log.read_text(), "Gaussian log\n")
            self.assertEqual(staged_fchk.read_text(), "Formatted checkpoint\n")

    def test_formchk_uses_absolute_paths_and_generates_in_workdir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source"
            workdir = Path(tmp) / "work"
            source.mkdir()
            log = source / "calculation.log"
            chk = source / "calculation.chk"
            log.write_text("Gaussian log\n", encoding="utf-8")
            chk.write_text("Binary checkpoint stand-in\n", encoding="utf-8")
            seen = {}

            def fake_formchk(command, **kwargs):
                seen["command"] = command
                seen["cwd"] = kwargs["cwd"]
                Path(command[2]).write_text(
                    "Generated formatted checkpoint\n", encoding="utf-8"
                )
                return SimpleNamespace(returncode=0, stderr="")

            with mock.patch.object(
                self.runner.subprocess, "run", side_effect=fake_formchk
            ):
                staged_log, staged_fchk = self.runner.stage_gaussian_inputs(
                    log, workdir, target_stem="s1"
                )

            self.assertEqual(seen["command"][0], "formchk")
            self.assertTrue(Path(seen["command"][1]).is_absolute())
            self.assertTrue(Path(seen["command"][2]).is_absolute())
            self.assertEqual(Path(seen["cwd"]), workdir.resolve())
            self.assertEqual(staged_log.name, "s1.log")
            self.assertEqual(staged_fchk.name, "s1.fchk")
            self.assertTrue(staged_fchk.is_file())

    def test_run_momap_resolves_input_before_changing_workdir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "work"
            workdir.mkdir()
            input_file = workdir / "momap.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            previous = Path.cwd()
            os.chdir(root)
            try:
                with mock.patch.object(
                    self.runner, "run_local", return_value=0
                ) as run_local:
                    self.runner.run_momap(Path("work") / "momap.inp")
            finally:
                os.chdir(previous)

        passed_input, passed_cwd = run_local.call_args.args
        self.assertEqual(Path(passed_input), input_file.resolve())
        self.assertEqual(Path(passed_cwd), workdir.resolve())

    def test_patched_launcher_uses_a_private_secure_temp_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            original = Path(tmp) / "momap"
            original.write_text(
                "#!/usr/bin/env python3\nprint('-machinefile')\n", encoding="utf-8"
            )

            with mock.patch.object(
                self.runner, "find_momap_bin", return_value=str(original)
            ):
                patched = Path(self.runner.patch_momap_for_mpi3())

            mode = stat.S_IMODE(patched.stat().st_mode)
            parent_mode = stat.S_IMODE(patched.parent.stat().st_mode)
            self.assertNotEqual(patched.parent, Path(tempfile.gettempdir()) / "momap_patched")
            self.assertEqual(mode & 0o077, 0)
            self.assertEqual(parent_mode & 0o077, 0)
            self.assertIn("--hostfile", patched.read_text(encoding="utf-8"))

    def test_slurm_rejects_injected_partition_and_invalid_process_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            input_file = workdir / "momap.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            with mock.patch.object(
                self.runner,
                "patch_momap_for_mpi3",
                return_value=str(workdir / "patched-momap"),
            ), mock.patch.object(
                self.runner.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=0, stdout="Submitted 1\n"),
            ):
                with self.assertRaises(ValueError):
                    self.runner.submit_slurm(
                        input_file,
                        workdir,
                        "gpu\n#SBATCH --time=99:00:00",
                        4,
                    )
                with self.assertRaises(ValueError):
                    self.runner.submit_slurm(input_file, workdir, "gpu", 0)

    def test_slurm_shell_quotes_paths_and_uses_current_python(self) -> None:
        with tempfile.TemporaryDirectory(prefix="momap work ") as tmp:
            workdir = Path(tmp)
            input_file = workdir / "momap input.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            patched = workdir / "patched momap"
            patched.write_text("#!/usr/bin/env python3\n", encoding="utf-8")

            with mock.patch.object(
                self.runner, "patch_momap_for_mpi3", return_value=str(patched)
            ) as patch_launcher, mock.patch.object(
                self.runner.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=0, stdout="Submitted 1\n"),
            ):
                rc = self.runner.submit_slurm(input_file, workdir, "gpu", 8)

            script = (workdir / "momap_job.slurm").read_text(encoding="utf-8")

        self.assertEqual(rc, 0)
        self.assertIn(f"cd -- {shlex.quote(str(workdir.resolve()))}", script)
        self.assertIn(shlex.quote(sys.executable), script)
        self.assertIn(shlex.quote(str(input_file.resolve())), script)
        self.assertIn(shlex.quote(str(patched.resolve())), script)
        patch_launcher.assert_called_once_with(workdir.resolve())


class TadfPipelineContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if str(TOOLS) not in sys.path:
            sys.path.insert(0, str(TOOLS))
        cls.tadf = load_module("momap_contract_tadf", TOOLS / "tadf.py")
        cls.stage4_path = ROOT / "tadf-screening" / "stage4_momap.py"
        cls.stage4 = load_module("momap_contract_stage4", cls.stage4_path)

    def test_molecule_id_cannot_escape_output_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "output"
            escaped_dir = root / "escaped"

            result = self.tadf.process_molecule(
                "../escaped",
                root / "missing-s0.log",
                root / "missing-s1.log",
                None,
                output_dir,
            )

        self.assertFalse(result["success"])
        self.assertIn("mol_id", result["error"])
        self.assertFalse(escaped_dir.exists())

    def test_stage4_reads_json_result_file_not_mixed_stdout(self) -> None:
        source = self.stage4_path.read_text(encoding="utf-8")
        self.assertIn("--json-output", source)
        self.assertNotIn("json.loads(result.stdout)", source)

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "output"
            output_dir.mkdir()
            seen_command = []

            def fake_run(command, **kwargs):
                seen_command.extend(command)
                output_arg = Path(command[command.index("--json-output") + 1])
                output_arg.write_text(
                    json.dumps({"mol_id": "m1", "success": True}),
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(
                    command,
                    returncode=0,
                    stdout="progress on stdout\nnot JSON\n",
                    stderr="",
                )

            with mock.patch.object(
                self.stage4.subprocess, "run", side_effect=fake_run
            ):
                result = self.stage4.run_single(
                    "m1", "s0.log", "s1.log", None, str(output_dir)
                )

        self.assertTrue(result["success"])
        self.assertIn("--json-output", seen_command)
        self.assertNotIn("--json", seen_command)

    def test_stage4_reads_failure_payload_from_json_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "output"

            def fake_run(command, **kwargs):
                output_arg = Path(command[command.index("--json-output") + 1])
                output_arg.write_text(
                    json.dumps({
                        "mol_id": "m1",
                        "success": False,
                        "error": "missing state 1 excitation",
                    }),
                    encoding="utf-8",
                )
                return subprocess.CompletedProcess(
                    command,
                    returncode=1,
                    stdout="mixed progress that is not JSON\n",
                    stderr="",
                )

            with mock.patch.object(
                self.stage4.subprocess, "run", side_effect=fake_run
            ):
                result = self.stage4.run_single(
                    "m1", "s0.log", "s1.log", None, str(output_dir)
                )

        self.assertFalse(result["success"])
        self.assertEqual(result["error"], "missing state 1 excitation")

    def test_tadf_cli_writes_json_output_file(self) -> None:
        source = (TOOLS / "tadf.py").read_text(encoding="utf-8")
        self.assertIn("--json-output", source)
        self.assertNotIn("parser.add_argument('--json'", source)

        with tempfile.TemporaryDirectory() as tmp:
            json_output = Path(tmp) / "result.json"
            argv = [
                "tadf.py",
                "m1",
                "--s0",
                "s0.log",
                "--s1",
                "s1.log",
                "--json-output",
                str(json_output),
            ]
            with mock.patch.object(
                self.tadf,
                "process_molecule",
                return_value={"mol_id": "m1", "success": True},
            ), mock.patch.object(sys, "argv", argv):
                rc = self.tadf.main()

            payload = json.loads(json_output.read_text(encoding="utf-8"))

        self.assertEqual(rc, 0)
        self.assertEqual(payload, {"mol_id": "m1", "success": True})

    def test_tadf_cli_writes_failure_json_when_pipeline_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            json_output = Path(tmp) / "failure.json"
            argv = [
                "tadf.py",
                "m1",
                "--s0",
                "s0.log",
                "--s1",
                "s1.log",
                "--json-output",
                str(json_output),
            ]
            with mock.patch.object(
                self.tadf,
                "process_molecule",
                side_effect=RuntimeError("missing evc.out"),
            ), mock.patch.object(sys, "argv", argv):
                try:
                    rc = self.tadf.main()
                except RuntimeError as exc:
                    self.fail(f"CLI leaked pipeline exception instead of writing JSON: {exc}")

            payload = json.loads(json_output.read_text(encoding="utf-8"))

        self.assertEqual(rc, 1)
        self.assertFalse(payload["success"])
        self.assertIn("missing evc.out", payload["error"])

    def test_no_soc_input_is_reported_as_not_computed(self) -> None:
        signature = inspect.signature(self.tadf.process_molecule)
        self.assertIn("hso_cm1", signature.parameters)
        source = (TOOLS / "tadf.py").read_text(encoding="utf-8")
        self.assertNotIn("Hso     = 1.0 cm-1", source)

        gaussian = textwrap.dedent(
            """\
             SCF Done:  E(RB3LYP) =  -100.000000000     A.U.
             Normal termination of Gaussian 16
             Normal termination of Gaussian 16
            """
        )
        s1_gaussian = gaussian + textwrap.dedent(
            """\
             Excited State   1:      Singlet-A      2.0000 eV
             Ground to excited state transition electric dipole moments (Au):
                  state          X           Y           Z        Dip. S.      Osc.
                 1        1.0000      2.0000      2.0000      9.0000      0.1000

             Ground to excited state transition electric dipole moments (Au):
                  state          X           Y           Z        Dip. S.      Osc.
                 1        0.0000      3.0000      4.0000     25.0000      0.2000

            """
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            logs = root / "logs"
            logs.mkdir()
            for name, content in (
                ("s0", gaussian),
                ("s1", s1_gaussian),
                ("t1", gaussian.replace("-100.000000000", "-99.950000000")),
            ):
                (logs / f"{name}.log").write_text(content, encoding="utf-8")
                (logs / f"{name}.fchk").write_text(
                    f"{name} formatted checkpoint\n", encoding="utf-8"
                )

            def fake_momap(input_file, workdir):
                input_name = Path(input_file).name
                if input_name == "momap_evc_s1.inp":
                    (Path(workdir) / "evc.out").write_text(
                        "1 # num of atoms\n3 # num of modes\n", encoding="utf-8"
                    )
                elif input_name == "momap_spec.inp":
                    (Path(workdir) / "spec.tvcf.spec.dat").write_text(
                        "0 0 0 440 0 0.5 0\n"
                        "0 0 0 460 0 1.0 0\n"
                        "0 0 0 480 0 0.5 0\n",
                        encoding="utf-8",
                    )
                return True

            with mock.patch.object(
                self.tadf, "run_momap_in_dir", side_effect=fake_momap
            ):
                result = self.tadf.process_molecule(
                    "m1",
                    logs / "s0.log",
                    logs / "s1.log",
                    logs / "t1.log",
                    root / "output",
                    hso_cm1=None,
                )

                result_with_soc = self.tadf.process_molecule(
                    "m2",
                    logs / "s0.log",
                    logs / "s1.log",
                    logs / "t1.log",
                    root / "output",
                    hso_cm1=12.5,
                )

            mol_dir = root / "output" / "m1"
            self.assertEqual(result["isc"]["status"], "not_computed")
            self.assertFalse((mol_dir / "momap_isc.inp").exists())
            self.assertTrue((mol_dir / "s0.log").exists())
            self.assertTrue((mol_dir / "s0.fchk").exists())
            self.assertAlmostEqual(result["EDMA"], 3.0 * self.tadf.AU2DEBYE)
            self.assertAlmostEqual(result["EDME"], 5.0 * self.tadf.AU2DEBYE)
            expected_s1 = -100.0 + 2.0 / self.tadf.HA2EV
            self.assertAlmostEqual(result["E_S1"], expected_s1)
            self.assertAlmostEqual(
                result["delta_EST_eV"],
                (expected_s1 - -99.95) * self.tadf.HA2EV,
            )
            self.assertEqual(result_with_soc["isc"]["status"], "computed")
            isc_input = (root / "output" / "m2" / "momap_isc.inp").read_text(
                encoding="utf-8"
            )
            self.assertIn("Hso     = 12.5 cm-1", isc_input)

    def test_core_spectrum_failure_or_missing_peak_fails_closed(self) -> None:
        s0_text = textwrap.dedent(
            """\
             SCF Done:  E(RB3LYP) =  -100.000000000     A.U.
             Normal termination of Gaussian 16
             Normal termination of Gaussian 16
            """
        )
        s1_text = s0_text + textwrap.dedent(
            """\
             Excited State   1:      Singlet-A      2.0000 eV
             Ground to excited state transition electric dipole moments (Au):
                  state          X           Y           Z        Dip. S.      Osc.
                 1        1.0000      0.0000      0.0000      1.0000      0.1000

             Ground to excited state transition electric dipole moments (Au):
                  state          X           Y           Z        Dip. S.      Osc.
                 1        0.0000      1.0000      0.0000      1.0000      0.2000

            """
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for stem, content in (("s0", s0_text), ("s1", s1_text)):
                (root / f"{stem}.log").write_text(content, encoding="utf-8")
                (root / f"{stem}.fchk").write_text(
                    "formatted checkpoint\n", encoding="utf-8"
                )

            def run_case(mol_id, spectrum_mode):
                def fake_momap(input_file, workdir):
                    name = Path(input_file).name
                    if name == "momap_evc_s1.inp":
                        (Path(workdir) / "evc.out").write_text(
                            "1 # num of atoms\n3 # num of modes\n",
                            encoding="utf-8",
                        )
                        return True
                    if name == "momap_spec.inp":
                        if spectrum_mode == "no_peak":
                            (Path(workdir) / "spec.tvcf.spec.dat").write_text(
                                "# no numeric spectrum rows\n", encoding="utf-8"
                            )
                            return True
                        return False
                    return True

                with mock.patch.object(
                    self.tadf, "run_momap_in_dir", side_effect=fake_momap
                ):
                    return self.tadf.process_molecule(
                        mol_id,
                        root / "s0.log",
                        root / "s1.log",
                        None,
                        root / "output",
                    )

            failed_run = run_case("failed-spectrum", "failed")
            missing_peak = run_case("missing-peak", "no_peak")

        for result in (failed_run, missing_peak):
            with self.subTest(mol_id=result["mol_id"]):
                self.assertFalse(result["success"])
                self.assertRegex(result["error"].lower(), r"spectrum|peak")


if __name__ == "__main__":
    unittest.main()
