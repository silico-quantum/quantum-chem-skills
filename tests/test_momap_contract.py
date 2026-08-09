from __future__ import annotations

import importlib.util
import hashlib
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
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "momap" / "tools"
MOMAP_BUILD = "2024A"
MOMAP_BANNER = "MOMAP version 2024A contract-test build"


def momap_expected_kwargs() -> dict[str, str]:
    launcher_bytes = (
        "#!/usr/bin/env python3\n"
        f"print({MOMAP_BANNER!r})\n"
        "# -machinefile\n"
    ).encode("utf-8")
    return {
        "expected_build": MOMAP_BUILD,
        "expected_launcher_sha256": hashlib.sha256(launcher_bytes).hexdigest(),
        "expected_version_banner": MOMAP_BANNER,
    }


def momap_cli_flags() -> list[str]:
    expected = momap_expected_kwargs()
    return [
        "--expected-build",
        expected["expected_build"],
        "--expected-launcher-sha256",
        expected["expected_launcher_sha256"],
        "--expected-version-banner",
        expected["expected_version_banner"],
    ]


def write_fake_momap_execution(
    stage: Path,
    *,
    banner: str = MOMAP_BANNER,
    build: str = MOMAP_BUILD,
    exit_code: int = 0,
) -> dict:
    """Write deterministic fake launcher evidence for pipeline contract tests."""
    original = stage / "contract-test-momap"
    original_bytes = (
        "#!/usr/bin/env python3\n"
        f"print({MOMAP_BANNER!r})\n"
        "# -machinefile\n"
    ).encode("utf-8")
    original.write_bytes(original_bytes)
    original.chmod(0o700)
    private_dir = stage / ".contract-test-momap-private"
    private_dir.mkdir(exist_ok=True)
    private_dir.chmod(0o700)
    patched = private_dir / "momap-patched"
    patched_bytes = original_bytes.replace(b"-machinefile", b"--hostfile")
    patched.write_bytes(patched_bytes)
    patched.chmod(0o700)
    run_log = stage / "momap_runner.log"
    run_log.write_text(banner + "\n", encoding="utf-8")
    return {
        "process_exit_code": exit_code,
        "run_log": str(run_log),
        "run_log_sha256": hashlib.sha256(run_log.read_bytes()).hexdigest(),
        "executable_identity": {
            "build": build,
            "version_banner": banner,
            "original_launcher": {
                "path": str(original),
                "sha256": hashlib.sha256(original_bytes).hexdigest(),
            },
            "patched_launcher": {
                "path": str(patched),
                "sha256": hashlib.sha256(patched_bytes).hexdigest(),
            },
            "mpi_patch": {
                "contract": "replace_exact_-machinefile_with_--hostfile",
                "replacement_count": 1,
                "original_preserved": True,
            },
        },
    }


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gaussian_chain(
    state: str,
    *,
    energy: float = -100.0,
    excitation_ev: float = 2.0,
    charge: int = 0,
    multiplicity: int | None = None,
    optimization_coordinates: tuple[tuple[float, float, float], ...] = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
    frequency_coordinates: tuple[tuple[float, float, float], ...] | None = None,
) -> str:
    """Return a minimal two-job Gaussian opt/freq chain for contract tests."""
    if multiplicity is None:
        multiplicity = 3 if state == "T1" else 1
    if frequency_coordinates is None:
        frequency_coordinates = optimization_coordinates

    def orientation(coordinates: tuple[tuple[float, float, float], ...]) -> str:
        rows = "\n".join(
            f"{index:19d}{atomic_number:11d}{0:12d}"
            f"{x:16.6f}{y:12.6f}{z:12.6f}"
            for index, (atomic_number, (x, y, z)) in enumerate(
                zip((6, 1), coordinates), start=1
            )
        )
        return textwrap.dedent(
            f"""\
             Standard orientation:
             ---------------------------------------------------------------------
             Center     Atomic      Atomic             Coordinates (Angstroms)
             Number     Number       Type             X           Y           Z
             ---------------------------------------------------------------------
            {rows}
             ---------------------------------------------------------------------
            """
        )

    def td_block(x: float, y: float, z: float, osc: float) -> str:
        dip_sq = x * x + y * y + z * z
        return textwrap.dedent(
            f"""\
             Excited State   1:      Singlet-A      {excitation_ev:.4f} eV
             This state for optimization and/or second-order correction.
             Ground to excited state transition electric dipole moments (Au):
                  state          X           Y           Z        Dip. S.      Osc.
                 1        {x:.4f}      {y:.4f}      {z:.4f}      {dip_sq:.4f}      {osc:.4f}

            """
        )

    first_td = td_block(1.0, 2.0, 2.0, 0.1) if state == "S1" else ""
    final_td = td_block(0.0, 3.0, 4.0, 0.2) if state == "S1" else ""
    return (
        f" Charge = {charge} Multiplicity = {multiplicity}\n"
        + orientation(optimization_coordinates)
        + f" SCF Done:  E(RB3LYP) =  {energy:.9f}     A.U.\n"
        + first_td
        + " Optimization completed.\n"
        + " Normal termination of Gaussian 16 at contract-test\n"
        + " --Link1--\n"
        + f" Charge = {charge} Multiplicity = {multiplicity}\n"
        + orientation(frequency_coordinates)
        + f" SCF Done:  E(RB3LYP) =  {energy:.9f}     A.U.\n"
        + final_td
        + " Harmonic frequencies (cm**-1)\n"
        + " Frequencies --   100.0  200.0  300.0\n"
        + " Normal termination of Gaussian 16 at contract-test\n"
    )


def gaussian_fchk(
    *,
    charge: int = 0,
    multiplicity: int = 1,
    coordinates_bohr: tuple[float, ...] | None = (
        0.0,
        0.0,
        0.0,
        1.8897261246257702,
        0.0,
        0.0,
    ),
) -> str:
    coordinate_section = ""
    if coordinates_bohr is not None:
        payload = " ".join(f"{value:.12E}" for value in coordinates_bohr)
        coordinate_section = (
            "Current cartesian coordinates              R   N=           "
            f"{len(coordinates_bohr)}\n {payload}\n"
        )
    return textwrap.dedent(
        f"""\
        Contract test formatted checkpoint
        Number of atoms                            I                2
        Charge                                     I                {charge}
        Multiplicity                               I                {multiplicity}
        Atomic numbers                             I   N=           2
                   6           1
        {coordinate_section}
        """
    )


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

    def test_gaussian_chain_requires_final_normal_opt_freq_and_state_continuity(self) -> None:
        validator = getattr(self.extract, "validate_gaussian_log_contract", None)
        self.assertIsNotNone(
            validator,
            "extract.py must validate the complete Gaussian opt/freq state chain",
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "s1.log"
            fchk = root / "s1.fchk"
            valid = gaussian_chain("S1")
            log.write_text(valid, encoding="utf-8")
            fchk.write_text(gaussian_fchk(), encoding="utf-8")

            metadata = validator(log, fchk, state_label="S1")
            self.assertEqual(metadata["normal_job_count"], 2)
            self.assertEqual(metadata["atomic_numbers"], [6, 1])
            self.assertEqual(metadata["target_root"], 1)

            invalid_logs = {
                "fatal tail": valid + " Error termination via Lnk1e\n",
                "missing optimization": valid.replace(" Optimization completed.\n", ""),
                "missing frequency": valid.replace(" Frequencies --   100.0  200.0  300.0\n", ""),
                "triplet mislabeled S1": valid.replace("Singlet-A", "Triplet-A"),
                "root switch": valid.replace(
                    "Excited State   1:      Singlet-A", "Excited State   2:      Singlet-A", 1
                ),
                "atom reorder": valid.replace(
                    "     2          1           0", "     2          8           0", 1
                ),
            }
            for case, content in invalid_logs.items():
                with self.subTest(case=case):
                    log.write_text(content, encoding="utf-8")
                    with self.assertRaises(ValueError):
                        validator(log, fchk, state_label="S1")

            log.write_text(valid, encoding="utf-8")
            fchk.write_text(gaussian_fchk().replace("           6           1", "           6           8"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "atom|Atomic"):
                validator(log, fchk, state_label="S1")

    def test_gaussian_contract_binds_final_frequency_geometry_to_fchk(self) -> None:
        validator = self.extract.validate_gaussian_log_contract
        final_coordinates = ((0.0, 0.0, 0.0), (1.25, -0.5, 0.125))
        matching_bohr = tuple(
            coordinate / 0.529177210903
            for atom in final_coordinates
            for coordinate in atom
        )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "s0.log"
            fchk = root / "s0.fchk"
            log.write_text(
                gaussian_chain(
                    "S0",
                    frequency_coordinates=final_coordinates,
                ),
                encoding="utf-8",
            )
            fchk.write_text(
                gaussian_fchk(coordinates_bohr=matching_bohr),
                encoding="utf-8",
            )

            metadata = validator(log, fchk, state_label="S0")
            self.assertEqual(metadata["geometry_angstrom"], [list(atom) for atom in final_coordinates])

            fchk.write_text(gaussian_fchk(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "geometry|coordinate"):
                validator(log, fchk, state_label="S0")

    def test_gaussian_contract_rejects_missing_fchk_cartesian_coordinates(self) -> None:
        validator = self.extract.validate_gaussian_log_contract
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "s0.log"
            fchk = root / "s0.fchk"
            log.write_text(gaussian_chain("S0"), encoding="utf-8")
            fchk.write_text(
                gaussian_fchk(coordinates_bohr=None),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "cartesian coordinates"):
                validator(log, fchk, state_label="S0")

    def test_gaussian_contract_rejects_nonfinite_fchk_cartesian_coordinates(self) -> None:
        validator = self.extract.validate_gaussian_log_contract
        nonfinite_coordinates = (0.0, 0.0, 0.0, math.nan, 0.0, 0.0)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            log = root / "s0.log"
            fchk = root / "s0.fchk"
            log.write_text(gaussian_chain("S0"), encoding="utf-8")
            fchk.write_text(
                gaussian_fchk(coordinates_bohr=nonfinite_coordinates),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "cartesian coordinates"):
                validator(log, fchk, state_label="S0")

    def test_spectrum_input_generator_refuses_existing_target(self) -> None:
        generator = self.extract.generate_spec_tvcf_input
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "momap_spec.inp"
            output.write_text("protected evidence\n", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                generator(
                    {"Ead": 0.1, "EDMA": 1.0, "EDME": 1.0}, output
                )
            self.assertEqual(
                output.read_text(encoding="utf-8"), "protected evidence\n"
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

    def test_gaussian_staging_refuses_existing_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            workdir = root / "work"
            source.mkdir()
            workdir.mkdir()
            (source / "s1.log").write_text("new log\n", encoding="utf-8")
            (source / "s1.fchk").write_text("new fchk\n", encoding="utf-8")
            (workdir / "s1.log").write_text("protected log\n", encoding="utf-8")
            (workdir / "s1.fchk").write_text("protected fchk\n", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                self.runner.stage_gaussian_inputs(
                    source / "s1.log", workdir, target_stem="s1"
                )

            self.assertEqual((workdir / "s1.log").read_text(), "protected log\n")
            self.assertEqual((workdir / "s1.fchk").read_text(), "protected fchk\n")

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
                    self.runner,
                    "run_local",
                    return_value={"process_exit_code": 0},
                ) as run_local:
                    self.runner.run_momap(
                        Path("work") / "momap.inp", **momap_expected_kwargs()
                    )
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
                identity = self.runner.patch_momap_for_mpi3(
                    expected_build=MOMAP_BUILD,
                    expected_launcher_sha256=hashlib.sha256(
                        original.read_bytes()
                    ).hexdigest(),
                )
                patched = Path(identity["patched_launcher"]["path"])

            mode = stat.S_IMODE(patched.stat().st_mode)
            parent_mode = stat.S_IMODE(patched.parent.stat().st_mode)
            self.assertNotEqual(patched.parent, Path(tempfile.gettempdir()) / "momap_patched")
            self.assertEqual(mode & 0o077, 0)
            self.assertEqual(parent_mode & 0o077, 0)
            self.assertIn("--hostfile", patched.read_text(encoding="utf-8"))

    def test_launcher_identity_binds_build_original_and_private_patch_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "momap"
            original.write_text(
                "#!/usr/bin/env python3\n"
                f"print({MOMAP_BANNER!r})\n"
                "# -machinefile\n",
                encoding="utf-8",
            )
            expected_sha256 = hashlib.sha256(original.read_bytes()).hexdigest()

            with mock.patch.object(
                self.runner, "find_momap_bin", return_value=str(original)
            ):
                identity = self.runner.patch_momap_for_mpi3(
                    expected_build=MOMAP_BUILD,
                    expected_launcher_sha256=expected_sha256,
                )

            patched = Path(identity["patched_launcher"]["path"])
            self.assertEqual(identity["build"], MOMAP_BUILD)
            self.assertEqual(
                identity["original_launcher"]["sha256"], expected_sha256
            )
            self.assertEqual(
                identity["patched_launcher"]["sha256"],
                hashlib.sha256(patched.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                identity["mpi_patch"]["contract"],
                "replace_exact_-machinefile_with_--hostfile",
            )
            self.assertEqual(identity["mpi_patch"]["replacement_count"], 1)
            self.assertEqual(stat.S_IMODE(patched.parent.stat().st_mode), 0o700)
            self.assertEqual(stat.S_IMODE(patched.stat().st_mode), 0o700)

    def test_launcher_identity_rejects_wrong_hash_symlink_and_wrong_build(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            original = root / "momap-real"
            original.write_text("#!/usr/bin/env python3\n# -machinefile\n")
            launcher_link = root / "momap"
            launcher_link.symlink_to(original)
            expected_sha256 = hashlib.sha256(original.read_bytes()).hexdigest()

            for case, found_path, build, digest in (
                ("wrong hash", original, MOMAP_BUILD, "0" * 64),
                ("symlink", launcher_link, MOMAP_BUILD, expected_sha256),
                ("wrong build", original, "2023B", expected_sha256),
            ):
                with self.subTest(case=case), mock.patch.object(
                    self.runner, "find_momap_bin", return_value=str(found_path)
                ), self.assertRaises((OSError, ValueError)):
                    self.runner.patch_momap_for_mpi3(
                        expected_build=build,
                        expected_launcher_sha256=digest,
                    )

    def test_local_run_requires_one_exact_version_banner_and_returns_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            input_file = workdir / "momap.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            launcher = workdir / "momap"
            launcher.write_text(
                "#!/usr/bin/env python3\n"
                f"print({MOMAP_BANNER!r})\n"
                "# -machinefile\n",
                encoding="utf-8",
            )
            launcher_sha256 = hashlib.sha256(launcher.read_bytes()).hexdigest()

            with mock.patch.object(
                self.runner, "find_momap_bin", return_value=str(launcher)
            ):
                execution = self.runner.run_local(
                    input_file,
                    workdir,
                    expected_build=MOMAP_BUILD,
                    expected_launcher_sha256=launcher_sha256,
                    expected_version_banner=MOMAP_BANNER,
                )

            self.assertEqual(execution["process_exit_code"], 0)
            self.assertEqual(
                execution["executable_identity"]["version_banner"],
                MOMAP_BANNER,
            )
            self.assertEqual(
                Path(execution["run_log"]).read_text(encoding="utf-8"),
                MOMAP_BANNER + "\n",
            )

    def test_local_run_rejects_missing_wrong_or_duplicate_version_banner(self) -> None:
        for case, printed_lines, expected_banner in (
            ("missing", [], MOMAP_BANNER),
            ("wrong", ["MOMAP version 2023B"], MOMAP_BANNER),
            ("duplicate", [MOMAP_BANNER, MOMAP_BANNER], MOMAP_BANNER),
        ):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                workdir = Path(tmp)
                input_file = workdir / "momap.inp"
                input_file.write_text("do_evc = 1\n", encoding="utf-8")
                launcher = workdir / "momap"
                body = "".join(f"print({line!r})\n" for line in printed_lines)
                launcher.write_text(
                    "#!/usr/bin/env python3\n" + body + "# -machinefile\n",
                    encoding="utf-8",
                )
                launcher_sha256 = hashlib.sha256(
                    launcher.read_bytes()
                ).hexdigest()
                with mock.patch.object(
                    self.runner, "find_momap_bin", return_value=str(launcher)
                ), self.assertRaisesRegex(ValueError, "banner|version"):
                    self.runner.run_local(
                        input_file,
                        workdir,
                        expected_build=MOMAP_BUILD,
                        expected_launcher_sha256=launcher_sha256,
                        expected_version_banner=expected_banner,
                    )

    def test_slurm_rejects_injected_partition_and_invalid_process_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            input_file = workdir / "momap.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            with mock.patch.object(
                self.runner,
                "patch_momap_for_mpi3",
                return_value={
                    "patched_launcher": {"path": str(workdir / "patched-momap")}
                },
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
                        **momap_expected_kwargs(),
                    )
                with self.assertRaises(ValueError):
                    self.runner.submit_slurm(
                        input_file,
                        workdir,
                        "gpu",
                        0,
                        **momap_expected_kwargs(),
                    )

    def test_slurm_shell_quotes_paths_and_uses_current_python(self) -> None:
        with tempfile.TemporaryDirectory(prefix="momap work ") as tmp:
            workdir = Path(tmp)
            input_file = workdir / "momap input.inp"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            patched = workdir / "patched momap"
            patched.write_text("#!/usr/bin/env python3\n", encoding="utf-8")

            with mock.patch.object(
                self.runner,
                "patch_momap_for_mpi3",
                return_value={"patched_launcher": {"path": str(patched)}},
            ) as patch_launcher, mock.patch.object(
                self.runner.subprocess,
                "run",
                return_value=SimpleNamespace(returncode=0, stdout="Submitted 1\n"),
            ):
                rc = self.runner.submit_slurm(
                    input_file,
                    workdir,
                    "gpu",
                    8,
                    **momap_expected_kwargs(),
                )

            script = (workdir / "momap_job.slurm").read_text(encoding="utf-8")

        self.assertEqual(rc, 0)
        self.assertIn(f"cd -- {shlex.quote(str(workdir.resolve()))}", script)
        self.assertIn(shlex.quote(sys.executable), script)
        self.assertIn(shlex.quote(str(input_file.resolve())), script)
        self.assertIn(shlex.quote(str(patched)), script)
        patch_launcher.assert_called_once_with(
            expected_build=MOMAP_BUILD,
            expected_launcher_sha256=momap_expected_kwargs()[
                "expected_launcher_sha256"
            ],
            temp_parent=workdir.resolve(),
        )


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
                **momap_expected_kwargs(),
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
                self.assertFalse(
                    output_arg.exists(),
                    "Stage 4 must give tadf.py a fresh JSON path for exclusive publication",
                )
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
                    "m1",
                    "s0.log",
                    "s1.log",
                    None,
                    str(output_dir),
                    **momap_expected_kwargs(),
                )

        self.assertTrue(result["success"])
        self.assertIn("--json-output", seen_command)
        self.assertIn("--expected-build", seen_command)
        self.assertIn("--expected-launcher-sha256", seen_command)
        self.assertIn("--expected-version-banner", seen_command)
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
                    "m1",
                    "s0.log",
                    "s1.log",
                    None,
                    str(output_dir),
                    **momap_expected_kwargs(),
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
                *momap_cli_flags(),
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

    def test_tadf_cli_refuses_to_overwrite_json_output_or_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            protected = root / "s0.log"
            protected.write_text("original Gaussian evidence\n", encoding="utf-8")
            argv = [
                "tadf.py",
                "m1",
                "--s0",
                str(protected),
                "--s1",
                "s1.log",
                "--json-output",
                str(protected),
                *momap_cli_flags(),
            ]
            with mock.patch.object(
                self.tadf,
                "process_molecule",
                return_value={"mol_id": "m1", "success": False},
            ), mock.patch.object(sys, "argv", argv):
                with self.assertRaises((FileExistsError, ValueError)):
                    self.tadf.main()

            self.assertEqual(
                protected.read_text(encoding="utf-8"),
                "original Gaussian evidence\n",
            )

    def test_tadf_cli_rejects_dangling_json_output_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            outside = root / "redirected.json"
            requested = root / "requested.json"
            requested.symlink_to(outside)
            argv = [
                "tadf.py",
                "m1",
                "--s0",
                "s0.log",
                "--s1",
                "s1.log",
                "--json-output",
                str(requested),
                *momap_cli_flags(),
            ]
            with mock.patch.object(sys, "argv", argv):
                with self.assertRaises((FileExistsError, ValueError)):
                    self.tadf.main()
            self.assertTrue(requested.is_symlink())
            self.assertFalse(outside.exists())

    def test_pipeline_rejects_cross_state_charge_or_atom_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "s0.log").write_text(gaussian_chain("S0"), encoding="utf-8")
            (root / "s0.fchk").write_text(gaussian_fchk(), encoding="utf-8")
            (root / "s1.log").write_text(
                gaussian_chain("S1", charge=1), encoding="utf-8"
            )
            (root / "s1.fchk").write_text(
                gaussian_fchk(charge=1), encoding="utf-8"
            )
            with mock.patch.object(
                self.tadf,
                "run_momap_in_dir",
                side_effect=AssertionError("MOMAP must not run for mixed states"),
            ):
                result = self.tadf.process_molecule(
                    "mixed",
                    root / "s0.log",
                    root / "s1.log",
                    None,
                    root / "output",
                    **momap_expected_kwargs(),
                )
            self.assertFalse(result["success"])
            self.assertRegex(result["error"].lower(), r"charge|atom|same system")

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
                *momap_cli_flags(),
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

        gaussian = gaussian_chain("S0")
        s1_gaussian = gaussian_chain("S1")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            logs = root / "logs"
            logs.mkdir()
            for name, content in (
                ("s0", gaussian),
                ("s1", s1_gaussian),
                ("t1", gaussian_chain("T1", energy=-99.95)),
            ):
                (logs / f"{name}.log").write_text(content, encoding="utf-8")
                (logs / f"{name}.fchk").write_text(
                    gaussian_fchk(multiplicity=3 if name == "t1" else 1),
                    encoding="utf-8",
                )

            def fake_momap(input_file, workdir, **_expected):
                stage = Path(workdir)
                if stage.name.startswith("evc_"):
                    (stage / "evc.out").write_text(
                        "1 # num of atoms\n3 # num of modes\n"
                        "Normal finish of evc calculation\n",
                        encoding="utf-8",
                    )
                    (stage / "evc.cart.dat").write_text(
                        f"{stage.name} Duschinsky data\n", encoding="utf-8"
                    )
                elif stage.name == "spectrum":
                    (stage / "spec.tvcf.log").write_text(
                        "Normal finish of spec_tvcf calculation\n",
                        encoding="utf-8",
                    )
                    (stage / "spec.tvcf.spec.dat").write_text(
                        "0.10354 2.81782 22727.27 440.0 0 0.5 0\n"
                        "0.09904 2.69531 21739.13 460.0 0 1.0 0\n"
                        "0.09494 2.58300 20833.33 480.0 0 0.5 0\n",
                        encoding="utf-8",
                    )
                elif stage.name == "isc":
                    (stage / "isc.tvcf.log").write_text(
                        "Normal finish of isc_tvcf calculation\n",
                        encoding="utf-8",
                    )
                    (stage / "isc.tvcf.fo.dat").write_text(
                        "ISC rate is 1.0E+06 s-1\n"
                        "RISC rate is 2.0E+04 s-1\n",
                        encoding="utf-8",
                    )
                return write_fake_momap_execution(stage)

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
                    **momap_expected_kwargs(),
                )

                result_with_soc = self.tadf.process_molecule(
                    "m2",
                    logs / "s0.log",
                    logs / "s1.log",
                    logs / "t1.log",
                    root / "output",
                    hso_cm1=12.5,
                    **momap_expected_kwargs(),
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
            isc_input = (
                root / "output" / "m2" / "isc" / "momap_isc.inp"
            ).read_text(encoding="utf-8")
            self.assertRegex(isc_input, r"Hso\s*=\s*12\.5(?:0+)?\s+cm-1")

    def test_core_spectrum_failure_or_missing_peak_fails_closed(self) -> None:
        s0_text = gaussian_chain("S0")
        s1_text = gaussian_chain("S1")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for stem, content in (("s0", s0_text), ("s1", s1_text)):
                (root / f"{stem}.log").write_text(content, encoding="utf-8")
                (root / f"{stem}.fchk").write_text(
                    gaussian_fchk(), encoding="utf-8"
                )

            def run_case(mol_id, spectrum_mode):
                def fake_momap(input_file, workdir, **_expected):
                    stage = Path(workdir)
                    if stage.name.startswith("evc_"):
                        (stage / "evc.out").write_text(
                            "1 # num of atoms\n3 # num of modes\n"
                            "Normal finish of evc calculation\n",
                            encoding="utf-8",
                        )
                        (stage / "evc.cart.dat").write_text(
                            f"{stage.name} Duschinsky data\n",
                            encoding="utf-8",
                        )
                        return write_fake_momap_execution(stage)
                    if stage.name == "spectrum":
                        if spectrum_mode == "no_peak":
                            (stage / "spec.tvcf.log").write_text(
                                "Normal finish of spec_tvcf calculation\n",
                                encoding="utf-8",
                            )
                            (stage / "spec.tvcf.spec.dat").write_text(
                                "# no numeric spectrum rows\n", encoding="utf-8"
                            )
                            return write_fake_momap_execution(stage)
                        return write_fake_momap_execution(stage, exit_code=1)
                    return write_fake_momap_execution(stage)

                with mock.patch.object(
                    self.tadf, "run_momap_in_dir", side_effect=fake_momap
                ):
                    return self.tadf.process_molecule(
                        mol_id,
                        root / "s0.log",
                        root / "s1.log",
                        None,
                        root / "output",
                        **momap_expected_kwargs(),
                    )

            failed_run = run_case("failed-spectrum", "failed")
            missing_peak = run_case("missing-peak", "no_peak")

        for result in (failed_run, missing_peak):
            with self.subTest(mol_id=result["mol_id"]):
                self.assertFalse(result["success"])
                self.assertRegex(result["error"].lower(), r"spectrum|peak")


class MomapStageArtifactContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        if str(TOOLS) not in sys.path:
            sys.path.insert(0, str(TOOLS))
        cls.tadf = load_module("momap_stage_contract_tadf", TOOLS / "tadf.py")
        cls.oled = load_module("momap_stage_contract_oled", TOOLS / "oled.py")
        cls.plot = load_module("momap_stage_contract_plot", TOOLS / "plot.py")

    def test_evc_acceptance_requires_fresh_nonempty_outputs_and_counts(self) -> None:
        validator = getattr(self.tadf, "validate_evc_stage", None)
        self.assertIsNotNone(
            validator,
            "tadf.py must expose a fail-closed EVC artifact validator",
        )

        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            started_ns = time.time_ns()
            (stage / "evc.out").write_text(
                "2 # num of atoms\n6 # num of modes\n"
                "Normal finish of evc calculation\n",
                encoding="utf-8",
            )
            (stage / "evc.cart.dat").write_text(
                "validated Duschinsky data\n", encoding="utf-8"
            )
            execution = write_fake_momap_execution(stage)

            def validate_current():
                return validator(
                    stage,
                    started_ns,
                    execution,
                    **momap_expected_kwargs(),
                )

            info = validate_current()
            self.assertEqual(info["natoms"], 2)
            self.assertEqual(info["nmodes"], 6)

            (stage / "evc.cart.dat").write_text("", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "empty|non-empty"):
                validate_current()

            (stage / "evc.cart.dat").write_text("data\n", encoding="utf-8")
            (stage / "evc.out").write_text(
                "0 # num of atoms\n0 # num of modes\n"
                "Normal finish of evc calculation\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "atom|mode"):
                validate_current()

            (stage / "evc.out").write_text(
                "2 # num of atoms\n6 # num of modes\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "normal completion|Normal finish"):
                validate_current()

            (stage / "evc.out").write_text(
                "2 # num of atoms\n6 # num of modes\n"
                "Normal finish of evc calculation\n",
                encoding="utf-8",
            )
            old_ns = started_ns - 10_000_000_000
            os.utime(stage / "evc.cart.dat", ns=(old_ns, old_ns))
            with self.assertRaisesRegex(ValueError, "fresh|older"):
                validate_current()

    def test_2024a_specific_parsers_require_explicit_build_binding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            spectrum = stage / "spec.tvcf.spec.dat"
            spectrum.write_text(
                "0.10354 2.81782 22727.27 440.0 0.1 0.5 0.0\n"
                "0.09904 2.69531 21739.13 460.0 0.2 1.0 0.0\n"
                "0.09494 2.58300 20833.33 480.0 0.1 0.5 0.0\n",
                encoding="utf-8",
            )
            rates = stage / "isc.tvcf.fo.dat"
            rates.write_text(
                "rate is 1.25E+06 s-1\nrate is 2.50E+04 s-1\n",
                encoding="utf-8",
            )

            with self.assertRaises((TypeError, ValueError)):
                self.tadf.parse_spec_output(spectrum)
            parsed = self.tadf.parse_spec_output(
                spectrum, expected_build=MOMAP_BUILD
            )
            self.assertEqual(
                parsed["spectrum_contract"], "MOMAP_2024A_exact_7_columns"
            )
            with self.assertRaises(ValueError):
                self.tadf.parse_spec_output(spectrum, expected_build="2023B")

            with self.assertRaises((TypeError, ValueError)):
                self.oled.parse_isc_output(rates)
            parsed_rates = self.oled.parse_isc_output(
                rates, expected_build=MOMAP_BUILD
            )
            self.assertEqual(
                parsed_rates["rate_parse_contract"],
                "MOMAP_2024A_ordered_first_ISC_second_RISC",
            )
            with self.assertRaises(ValueError):
                self.oled.parse_isc_output(rates, expected_build="2023B")

    def test_stage_receipt_requires_and_records_executable_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            input_file = stage / "momap.inp"
            output_file = stage / "evc.out"
            input_file.write_text("do_evc = 1\n", encoding="utf-8")
            output_file.write_text(
                "Normal finish of evc calculation\n", encoding="utf-8"
            )
            execution = write_fake_momap_execution(stage)

            receipt_path = self.tadf.write_stage_receipt(
                stage,
                "evc_s1_s0",
                time.time_ns(),
                [input_file],
                [output_file, Path(execution["run_log"])],
                {"process_exit_code": 0},
                executable_identity=execution["executable_identity"],
            )
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            self.assertEqual(receipt["executable_identity"]["build"], MOMAP_BUILD)
            self.assertEqual(
                receipt["executable_identity"]["version_banner"], MOMAP_BANNER
            )
            self.assertEqual(
                receipt["executable_identity"]["original_launcher"]["sha256"],
                momap_expected_kwargs()["expected_launcher_sha256"],
            )

    def test_spectrum_acceptance_requires_strict_complete_physical_table(self) -> None:
        validator = getattr(self.tadf, "validate_spectrum_stage", None)
        self.assertIsNotNone(validator)
        valid_rows = (
            "0.10354 2.81782 22727.27 440.0 0.1 0.5 0.0\n"
            "0.09904 2.69531 21739.13 460.0 0.2 1.0 0.0\n"
            "0.09494 2.58300 20833.33 480.0 0.1 0.5 0.0\n"
        )

        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)

            def validate(table: str, log: str = "Normal finish of spec_tvcf calculation\n"):
                started_ns = time.time_ns()
                (stage / "spec.tvcf.log").write_text(log, encoding="utf-8")
                (stage / "spec.tvcf.spec.dat").write_text(table, encoding="utf-8")
                execution = write_fake_momap_execution(stage)
                return validator(
                    stage,
                    started_ns,
                    execution,
                    **momap_expected_kwargs(),
                )

            accepted = validate(valid_rows)
            self.assertEqual(accepted["data_points"], 3)
            self.assertAlmostEqual(accepted["peak_wavelength"], 460.0)

            invalid_tables = {
                "single truncated row": valid_rows.splitlines(keepends=True)[0],
                "too few columns": (
                    "0.10354 2.81782 22727.27 440.0 0.1 0.5\n"
                    + "".join(valid_rows.splitlines(keepends=True)[1:])
                ),
                "extra column": valid_rows.replace(
                    "0.1 0.5 0.0\n", "0.1 0.5 0.0 99\n", 1
                ),
                "trailing garbage": valid_rows + "partial trailing output\n",
                "non-finite": valid_rows.replace("2.69531", "NaN", 1),
                "non-monotonic energy": (
                    "0.10354 2.80000 22582.0 442.80 0.1 0.5 0.0\n"
                    "0.10400 2.90000 23390.0 427.53 0.2 1.0 0.0\n"
                    "0.09900 2.70000 21776.0 459.20 0.1 0.5 0.0\n"
                ),
                "inconsistent eV nm": valid_rows.replace("480.0", "800.0", 1),
                "inconsistent Hartree": valid_rows.replace("0.10354", "999", 1),
                "inconsistent wavenumber": valid_rows.replace("22727.27", "1", 1),
                "negative absorption": valid_rows.replace("0.1 0.5 0.0", "-0.1 0.5 0.0", 1),
                "negative emission": valid_rows.replace("0.2 1.0 0.0", "0.2 -1.0 0.0", 1),
                "negative final intensity": valid_rows.replace("0.1 0.5 0.0", "0.1 0.5 -1.0", 1),
                "fatal diagnostic comment": "# FATAL: truncated output\n" + valid_rows,
            }
            for case, table in invalid_tables.items():
                with self.subTest(case=case), self.assertRaises(ValueError):
                    validate(table)

            with self.assertRaisesRegex(ValueError, "fatal|diagnostic"):
                validate(
                    valid_rows,
                    "FATAL: worker aborted\n"
                    "Normal finish of spec_tvcf calculation\n",
                )

            for invalid_marker in (
                "did not reach Normal finish of spec_tvcf calculation\n",
                "Normal finish of spec_tvcf calculation\n"
                "Normal finish of spec_tvcf calculation\n",
            ):
                with self.subTest(marker=invalid_marker), self.assertRaisesRegex(
                    ValueError, "marker|completion|exactly"
                ):
                    validate(valid_rows, invalid_marker)

    def test_isc_acceptance_requires_fresh_marker_finite_rates_and_units(self) -> None:
        validator = getattr(self.tadf, "validate_isc_stage", None)
        self.assertIsNotNone(
            validator,
            "tadf.py must expose a fail-closed ISC artifact validator",
        )

        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)

            def write_valid():
                started = time.time_ns()
                (stage / "isc.tvcf.log").write_text(
                    "Normal finish of isc_tvcf calculation\n", encoding="utf-8"
                )
                (stage / "isc.tvcf.fo.dat").write_text(
                    "ISC rate is 1.25E+06 s-1\n"
                    "RISC rate is 2.50E+04 s^-1\n",
                    encoding="utf-8",
                )
                return started, write_fake_momap_execution(stage)

            started_ns, execution = write_valid()
            rates = validator(
                stage, started_ns, execution, **momap_expected_kwargs()
            )
            self.assertEqual(rates["k_ISC_s-1"], 1.25e6)
            self.assertEqual(rates["k_RISC_s-1"], 2.50e4)

            started_ns = time.time_ns()
            (stage / "isc.tvcf.log").write_text(
                "Normal finish of isc_tvcf calculation\n", encoding="utf-8"
            )
            (stage / "isc.tvcf.fo.dat").write_text(
                "rate is 1.25E+06 s-1\nrate is 2.50E+04 s-1\n",
                encoding="utf-8",
            )
            execution = write_fake_momap_execution(stage)
            ordered_rates = validator(
                stage, started_ns, execution, **momap_expected_kwargs()
            )
            self.assertEqual(ordered_rates["k_ISC_s-1"], 1.25e6)
            self.assertEqual(ordered_rates["k_RISC_s-1"], 2.50e4)

            for case, log_text, rate_text in (
                ("empty", "Normal finish of isc_tvcf calculation\n", ""),
                (
                    "missing marker",
                    "calculation ended without marker\n",
                    "ISC rate is 1.25E+06 s-1\n",
                ),
                (
                    "missing unit",
                    "Normal finish of isc_tvcf calculation\n",
                    "ISC rate is 1.25E+06\n",
                ),
                (
                    "non-finite",
                    "Normal finish of isc_tvcf calculation\n",
                    "ISC rate is NaN s-1\n",
                ),
            ):
                with self.subTest(case=case):
                    started_ns = time.time_ns()
                    (stage / "isc.tvcf.log").write_text(
                        log_text, encoding="utf-8"
                    )
                    (stage / "isc.tvcf.fo.dat").write_text(
                        rate_text, encoding="utf-8"
                    )
                    execution = write_fake_momap_execution(stage)
                    with self.assertRaises(ValueError):
                        validator(
                            stage,
                            started_ns,
                            execution,
                            **momap_expected_kwargs(),
                        )

            started_ns, execution = write_valid()
            old_ns = started_ns - 10_000_000_000
            os.utime(stage / "isc.tvcf.fo.dat", ns=(old_ns, old_ns))
            with self.assertRaisesRegex(ValueError, "fresh|older"):
                validator(
                    stage,
                    started_ns,
                    execution,
                    **momap_expected_kwargs(),
                )

    def test_isc_rate_records_reject_duplicates_and_mixed_contracts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            stage = Path(tmp)
            log = stage / "isc.tvcf.log"
            rates = stage / "isc.tvcf.fo.dat"

            invalid_outputs = {
                "duplicate same ISC": (
                    "ISC rate is 1.0E+06 s-1\n"
                    "ISC rate is 1.0E+06 s-1\n"
                    "RISC rate is 2.0E+04 s-1\n"
                ),
                "duplicate conflicting ISC": (
                    "ISC rate is 1.0E+06 s-1\n"
                    "ISC rate is 3.0E+06 s-1\n"
                    "RISC rate is 2.0E+04 s-1\n"
                ),
                "duplicate RISC": (
                    "ISC rate is 1.0E+06 s-1\n"
                    "RISC rate is 2.0E+04 s-1\n"
                    "RISC rate is 2.0E+04 s-1\n"
                ),
                "explicit plus unlabeled": (
                    "ISC rate is 1.0E+06 s-1\n"
                    "RISC rate is 2.0E+04 s-1\n"
                    "rate is 9.0E+09 s-1\n"
                ),
                "partial explicit plus fallback": (
                    "ISC rate is 1.0E+06 s-1\n"
                    "rate is 2.0E+04 s-1\n"
                ),
                "three unlabeled": (
                    "rate is 1.0E+06 s-1\n"
                    "rate is 2.0E+04 s-1\n"
                    "rate is 9.0E+09 s-1\n"
                ),
            }

            for case, content in invalid_outputs.items():
                with self.subTest(case=case):
                    started_ns = time.time_ns()
                    log.write_text(
                        "Normal finish of isc_tvcf calculation\n",
                        encoding="utf-8",
                    )
                    rates.write_text(content, encoding="utf-8")
                    execution = write_fake_momap_execution(stage)
                    with self.assertRaises(ValueError):
                        self.tadf.validate_isc_stage(
                            stage,
                            started_ns,
                            execution,
                            **momap_expected_kwargs(),
                        )
                    with self.assertRaises(ValueError):
                        self.oled.parse_isc_output(rates)

    def test_standalone_isc_generator_validates_inputs_and_refuses_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evc = root / "evc.cart.dat"
            evc.write_text("Duschinsky data\n", encoding="utf-8")
            output = root / "momap_isc.inp"

            generated = self.oled.generate_isc_input(
                evc, 0.01, 5.0, output, temp=300, tmax=5000
            )
            self.assertEqual(Path(generated), output.resolve())

            with self.assertRaises(FileExistsError):
                self.oled.generate_isc_input(
                    evc, 0.01, 5.0, output, temp=300, tmax=5000
                )

            output.unlink()
            for kwargs in (
                {"ead_au": 0.0},
                {"ead_au": float("nan")},
                {"temp": 0.0},
                {"temp": float("inf")},
                {"tmax": 0.0},
            ):
                params = {
                    "evc_dat": evc,
                    "ead_au": 0.01,
                    "hso_cm1": 5.0,
                    "output": output,
                    "temp": 300,
                    "tmax": 5000,
                }
                params.update(kwargs)
                with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                    self.oled.generate_isc_input(**params)

            evc.write_text("", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "empty|non-empty"):
                self.oled.generate_isc_input(evc, 0.01, 5.0, output)

    def test_plot_output_path_must_not_exist(self) -> None:
        guard = getattr(self.plot, "ensure_new_output_path", None)
        self.assertIsNotNone(
            guard,
            "plot.py must expose the same no-overwrite guard used by both renderers",
        )
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "spectrum.png"
            self.assertEqual(guard(output), output.resolve())
            output.write_bytes(b"existing image")
            with self.assertRaises(FileExistsError):
                guard(output)

    def test_pipeline_uses_separate_stage_directories_receipts_and_signed_gap(self) -> None:
        gaussian = gaussian_chain("S0")
        s1_gaussian = gaussian_chain("S1")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            logs = root / "logs"
            logs.mkdir()
            for name, content in (
                ("s0", gaussian),
                ("s1", s1_gaussian),
                ("t1", gaussian_chain("T1", energy=-99.95)),
            ):
                (logs / f"{name}.log").write_text(content, encoding="utf-8")
                (logs / f"{name}.fchk").write_text(
                    gaussian_fchk(multiplicity=3 if name == "t1" else 1),
                    encoding="utf-8",
                )

            def fake_momap(input_file, workdir, **_expected):
                stage = Path(workdir)
                if stage.name.startswith("evc_"):
                    (stage / "evc.out").write_text(
                        "1 # num of atoms\n3 # num of modes\n"
                        "Normal finish of evc calculation\n",
                        encoding="utf-8",
                    )
                    (stage / "evc.cart.dat").write_text(
                        f"{stage.name} Duschinsky data\n", encoding="utf-8"
                    )
                elif stage.name == "spectrum":
                    (stage / "spec.tvcf.log").write_text(
                        "Normal finish of spec_tvcf calculation\n", encoding="utf-8"
                    )
                    (stage / "spec.tvcf.spec.dat").write_text(
                        "0.10354 2.81782 22727.27 440.0 0 0.5 0\n"
                        "0.09904 2.69531 21739.13 460.0 0 1.0 0\n"
                        "0.09494 2.58300 20833.33 480.0 0 0.5 0\n",
                        encoding="utf-8",
                    )
                elif stage.name == "isc":
                    (stage / "isc.tvcf.log").write_text(
                        "Normal finish of isc_tvcf calculation\n", encoding="utf-8"
                    )
                    (stage / "isc.tvcf.fo.dat").write_text(
                        "ISC rate is 1.0E+06 s-1\n"
                        "RISC rate is 2.0E+04 s-1\n",
                        encoding="utf-8",
                    )
                return write_fake_momap_execution(stage)

            with mock.patch.object(
                self.tadf, "run_momap_in_dir", side_effect=fake_momap
            ):
                result = self.tadf.process_molecule(
                    "m1",
                    logs / "s0.log",
                    logs / "s1.log",
                    logs / "t1.log",
                    root / "output",
                    hso_cm1=12.5,
                    **momap_expected_kwargs(),
                )

            mol_dir = root / "output" / "m1"
            self.assertTrue(result["success"])
            self.assertEqual(result["isc"]["status"], "computed")
            self.assertEqual(result["E_T1_energy_type"], "Gaussian SCF total energy")
            self.assertGreater(result["delta_EST_signed_eV"], 0)
            self.assertEqual(result["delta_EST_eV"], result["delta_EST_signed_eV"])
            self.assertNotIn("abs(E_S1 - E_T1)", (TOOLS / "tadf.py").read_text())
            self.assertTrue((mol_dir / "evc_s1_s0" / "stage_receipt.json").is_file())
            self.assertTrue((mol_dir / "spectrum" / "stage_receipt.json").is_file())
            self.assertTrue((mol_dir / "evc_s1_t1" / "stage_receipt.json").is_file())
            self.assertTrue((mol_dir / "isc" / "stage_receipt.json").is_file())
            for receipt_path in mol_dir.glob("*/stage_receipt.json"):
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                self.assertEqual(receipt["status"], "accepted")
                self.assertEqual(receipt["metadata"]["process_exit_code"], 0)
                identity = receipt["executable_identity"]
                self.assertEqual(identity["build"], MOMAP_BUILD)
                self.assertEqual(identity["version_banner"], MOMAP_BANNER)
                self.assertEqual(
                    identity["original_launcher"]["sha256"],
                    momap_expected_kwargs()["expected_launcher_sha256"],
                )
                self.assertRegex(
                    identity["patched_launcher"]["sha256"], r"^[0-9a-f]{64}$"
                )
                self.assertTrue(
                    any(
                        artifact["path"] == "momap_runner.log"
                        for artifact in receipt["outputs"]
                    )
                )
                for artifact in [*receipt["inputs"], *receipt["outputs"]]:
                    self.assertRegex(artifact["sha256"], r"^[0-9a-f]{64}$")
                    self.assertGreater(artifact["size_bytes"], 0)
            self.assertNotEqual(
                (mol_dir / "evc_s1_s0" / "evc.cart.dat").read_text(),
                (mol_dir / "evc_s1_t1" / "evc.cart.dat").read_text(),
            )

    def test_invalid_isc_artifacts_fail_the_requested_pipeline(self) -> None:
        gaussian = gaussian_chain("S0")
        s1_gaussian = gaussian_chain("S1")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, content in (
                ("s0", gaussian),
                ("s1", s1_gaussian),
                ("t1", gaussian_chain("T1", energy=-99.95)),
            ):
                (root / f"{name}.log").write_text(content, encoding="utf-8")
                (root / f"{name}.fchk").write_text(
                    gaussian_fchk(multiplicity=3 if name == "t1" else 1),
                    encoding="utf-8",
                )

            for mode in (
                "missing",
                "empty",
                "old",
                "no_unit",
                "parse_failure",
                "missing_marker",
            ):
                def fake_momap(input_file, workdir, *, _mode=mode, **_expected):
                    stage = Path(workdir)
                    if stage.name.startswith("evc_"):
                        (stage / "evc.out").write_text(
                            "1 # num of atoms\n3 # num of modes\n"
                            "Normal finish of evc calculation\n",
                            encoding="utf-8",
                        )
                        (stage / "evc.cart.dat").write_text(
                            "Duschinsky data\n", encoding="utf-8"
                        )
                    elif stage.name == "spectrum":
                        (stage / "spec.tvcf.log").write_text(
                            "Normal finish of spec_tvcf calculation\n",
                            encoding="utf-8",
                        )
                        (stage / "spec.tvcf.spec.dat").write_text(
                            "0.10354 2.81782 22727.27 440.0 0 0.5 0\n"
                            "0.09904 2.69531 21739.13 460.0 0 1.0 0\n"
                            "0.09494 2.58300 20833.33 480.0 0 0.5 0\n",
                            encoding="utf-8",
                        )
                    elif stage.name == "isc" and _mode != "missing":
                        log_text = (
                            "calculation stopped without completion\n"
                            if _mode == "missing_marker"
                            else "Normal finish of isc_tvcf calculation\n"
                        )
                        (stage / "isc.tvcf.log").write_text(
                            log_text, encoding="utf-8"
                        )
                        rate_text = {
                            "empty": "",
                            "no_unit": (
                                "ISC rate is 1.0E+06\nRISC rate is 2.0E+04\n"
                            ),
                            "parse_failure": (
                                "ISC rate is nonsense s-1\n"
                                "RISC rate is 2.0E+04 s-1\n"
                            ),
                        }.get(
                            _mode,
                            "ISC rate is 1.0E+06 s-1\n"
                            "RISC rate is 2.0E+04 s-1\n",
                        )
                        rate_file = stage / "isc.tvcf.fo.dat"
                        rate_file.write_text(rate_text, encoding="utf-8")
                        if _mode == "old":
                            old_ns = time.time_ns() - 10_000_000_000
                            os.utime(rate_file, ns=(old_ns, old_ns))
                    return write_fake_momap_execution(stage)

                with self.subTest(mode=mode), mock.patch.object(
                    self.tadf, "run_momap_in_dir", side_effect=fake_momap
                ):
                    result = self.tadf.process_molecule(
                        f"m-{mode}",
                        root / "s0.log",
                        root / "s1.log",
                        root / "t1.log",
                        root / "output",
                        hso_cm1=12.5,
                        **momap_expected_kwargs(),
                    )

                self.assertFalse(result["success"])
                self.assertEqual(result["isc"]["status"], "failed")
                self.assertIn("ISC", result["error"])

    def test_requested_isc_rejects_nonpositive_signed_state_gap(self) -> None:
        gaussian = gaussian_chain("S0")
        s1_gaussian = gaussian_chain("S1")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, content in (
                ("s0", gaussian),
                ("s1", s1_gaussian),
                ("t1", gaussian_chain("T1", energy=-99.8)),
            ):
                (root / f"{name}.log").write_text(content, encoding="utf-8")
                (root / f"{name}.fchk").write_text(
                    gaussian_fchk(multiplicity=3 if name == "t1" else 1),
                    encoding="utf-8",
                )

            with mock.patch.object(
                self.tadf,
                "run_momap_in_dir",
                side_effect=AssertionError("MOMAP must not run for invalid ordering"),
            ):
                result = self.tadf.process_molecule(
                    "bad-order",
                    root / "s0.log",
                    root / "s1.log",
                    root / "t1.log",
                    root / "output",
                    hso_cm1=12.5,
                    **momap_expected_kwargs(),
                )

        self.assertFalse(result["success"])
        self.assertEqual(result["isc"]["status"], "failed")
        self.assertLess(result["delta_EST_signed_eV"], 0)
        self.assertEqual(result["state_order_status"], "unexpected_S1_not_above_T1")
        self.assertIn("Signed S1-T1 gap is non-positive", result["error"])
        self.assertNotIn("abs(E_S1 - E_T1)", (TOOLS / "tadf.py").read_text())

    def test_requested_isc_failure_keeps_global_success_false(self) -> None:
        source = (TOOLS / "tadf.py").read_text(encoding="utf-8")
        self.assertIn("validate_isc_stage", source)
        self.assertNotRegex(
            source,
            r"'status':\s*'computed'\s*if\s+ok_isc",
            "A zero process status alone must never produce computed",
        )


if __name__ == "__main__":
    unittest.main()
