from __future__ import annotations

import importlib.util
import json
import py_compile
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PYSCF_SKILL = ROOT / "pyscf"
RUNNER = PYSCF_SKILL / "scripts" / "run_safe_dft_tda.py"


def load_runner() -> types.ModuleType:
    if not RUNNER.is_file():
        raise AssertionError(f"supported runner is missing: {RUNNER}")
    spec = importlib.util.spec_from_file_location("run_safe_dft_tda", RUNNER)
    if spec is None or spec.loader is None:
        raise AssertionError("could not create a module spec for the runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PySCFContractTests(unittest.TestCase):
    def test_supported_runner_imports_without_optional_packages(self) -> None:
        module = load_runner()
        self.assertTrue(callable(module.main))
        self.assertNotIn("pyscf", module.__dict__)

    def test_config_is_explicit_and_closed_shell_only(self) -> None:
        module = load_runner()
        valid = {
            "atom": "H 0 0 0; H 0 0 0.74",
            "unit": "Angstrom",
            "basis": "sto-3g",
            "charge": 0,
            "spin": 0,
            "xc": "pbe0",
            "grid_level": 3,
            "conv_tol": 1e-9,
            "max_cycle": 80,
            "nstates": 2,
        }
        normalized = module.validate_config(valid)
        self.assertEqual(normalized, valid)

        for key, value in (
            ("unit", "unknown"),
            ("spin", 1),
            ("nstates", 0),
            ("conv_tol", float("nan")),
        ):
            broken = dict(valid)
            broken[key] = value
            with self.subTest(key=key, value=value):
                with self.assertRaises(ValueError):
                    module.validate_config(broken)

    def test_scf_acceptance_requires_finite_orbitals_and_checkpoint(self) -> None:
        module = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "scf.chk"
            checkpoint.write_bytes(b"checkpoint")
            mf = types.SimpleNamespace(
                converged=True,
                mol=types.SimpleNamespace(nelectron=2),
                e_tot=-1.1,
                mo_energy=[-0.5, 0.2],
                mo_coeff=[[1.0, 0.0], [0.0, 1.0]],
                mo_occ=[2.0, 0.0],
            )
            module.validate_scf_result(mf, -1.1, checkpoint)

            mf.mo_coeff[0][0] = float("nan")
            with self.assertRaises(ValueError):
                module.validate_scf_result(mf, -1.1, checkpoint)

            mf.mo_coeff[0][0] = 1.0
            checkpoint.write_bytes(b"")
            with self.assertRaises(ValueError):
                module.validate_scf_result(mf, -1.1, checkpoint)

    def test_scf_kernel_energy_must_match_final_mean_field_energy(self) -> None:
        module = load_runner()
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "scf.chk"
            checkpoint.write_bytes(b"checkpoint")
            mf = types.SimpleNamespace(
                converged=True,
                mol=types.SimpleNamespace(nelectron=2),
                e_tot=-1.1,
                mo_energy=[-0.5, 0.2],
                mo_coeff=[[1.0, 0.0], [0.0, 1.0]],
                mo_occ=[2.0, 0.0],
            )

            with self.assertRaisesRegex(ValueError, "does not match"):
                module.validate_scf_result(mf, -1.100001, checkpoint)

    def test_scf_occupations_must_be_physical_and_sum_to_electron_count(self) -> None:
        module = load_runner()
        invalid_occupations = (
            ([2.0, float("nan")], "non-finite"),
            ([-0.1, 2.1], "range"),
            ([2.0, 1.0], "electron count"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "scf.chk"
            checkpoint.write_bytes(b"checkpoint")
            for occupations, message in invalid_occupations:
                with self.subTest(occupations=occupations):
                    mf = types.SimpleNamespace(
                        converged=True,
                        mol=types.SimpleNamespace(nelectron=2),
                        e_tot=-1.1,
                        mo_energy=[-0.5, 0.2],
                        mo_coeff=[[1.0, 0.0], [0.0, 1.0]],
                        mo_occ=occupations,
                    )
                    with self.assertRaisesRegex(ValueError, message):
                        module.validate_scf_result(mf, -1.1, checkpoint)

    def test_scf_orbital_arrays_require_consistent_physical_shapes(self) -> None:
        module = load_runner()
        invalid_arrays = (
            (
                "two-dimensional orbital energies",
                [[-0.5, 0.2]],
                [[1.0, 0.0], [0.0, 1.0]],
                [2.0, 0.0],
            ),
            (
                "two-dimensional occupations",
                [-0.5, 0.2],
                [[1.0, 0.0], [0.0, 1.0]],
                [[2.0, 0.0]],
            ),
            (
                "occupation length mismatch",
                [-0.5, 0.2],
                [[1.0, 0.0], [0.0, 1.0]],
                [1.0, 1.0, 0.0],
            ),
            (
                "one-dimensional coefficients",
                [-0.5, 0.2],
                [1.0, 0.0],
                [2.0, 0.0],
            ),
            (
                "coefficient orbital-count mismatch",
                [-0.5, 0.2],
                [[1.0], [0.0]],
                [2.0, 0.0],
            ),
            (
                "ragged coefficient matrix",
                [-0.5, 0.2],
                [[1.0, 0.0], [0.0]],
                [2.0, 0.0],
            ),
            (
                "zero AO dimension",
                [-0.5, 0.2],
                [],
                [2.0, 0.0],
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            checkpoint = Path(temporary) / "scf.chk"
            checkpoint.write_bytes(b"checkpoint")
            for label, energies, coefficients, occupations in invalid_arrays:
                with self.subTest(label=label):
                    mf = types.SimpleNamespace(
                        converged=True,
                        mol=types.SimpleNamespace(nelectron=2),
                        e_tot=-1.1,
                        mo_energy=energies,
                        mo_coeff=coefficients,
                        mo_occ=occupations,
                    )
                    with self.assertRaisesRegex(ValueError, "shape"):
                        module.validate_scf_result(mf, -1.1, checkpoint)

    def test_scf_stability_gate_requires_internal_and_external_stability(self) -> None:
        module = load_runner()

        class MeanField:
            def __init__(self, internal: bool, external: bool):
                self.status = (internal, external)
                self.calls = []

            def stability(self, **kwargs):
                self.calls.append(kwargs)
                return object(), object(), *self.status

        stable = MeanField(True, True)
        self.assertEqual(
            module.validate_scf_stability(stable),
            {"internal": True, "external": True},
        )
        self.assertEqual(
            stable.calls,
            [{"internal": True, "external": True, "return_status": True}],
        )

        for internal, external in ((False, True), (True, False)):
            with self.subTest(internal=internal, external=external):
                with self.assertRaisesRegex(ValueError, "stability"):
                    module.validate_scf_stability(MeanField(internal, external))

    def test_tda_acceptance_requires_exact_finite_root_vectors(self) -> None:
        module = load_runner()
        module.validate_tda_result(
            converged=[True, True],
            energies=[0.2, 0.3],
            oscillator_strengths=[0.0, 0.1],
            nstates=2,
        )

        invalid_cases = (
            ([True], [0.2], [0.1], "missing requested root"),
            ([True, False], [0.2, 0.3], [0.1, 0.2], "unconverged root"),
            ([True, True], [0.2], [0.1, 0.2], "energy shape mismatch"),
            ([True, True], [0.2, 0.3], [0.1], "strength shape mismatch"),
            ([True, True], [0.2, float("nan")], [0.1, 0.2], "non-finite energy"),
            ([True, True], [0.0, 0.3], [0.1, 0.2], "zero excitation energy"),
            ([True, True], [-0.2, 0.3], [0.1, 0.2], "negative excitation energy"),
            ([True, True], [0.2, 0.3], [-0.1, 0.2], "negative strength"),
        )
        for converged, energies, strengths, label in invalid_cases:
            with self.subTest(label=label):
                with self.assertRaises(ValueError):
                    module.validate_tda_result(
                        converged, energies, strengths, nstates=2
                    )

    def test_run_directory_is_fresh_and_manifest_transitions_are_atomic(self) -> None:
        module = load_runner()
        config_bytes = b'{"unit":"Angstrom"}\n'
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "run"
            manifest = module.initialize_run(
                output_dir=output_dir,
                config_bytes=config_bytes,
                config_sha256=module.sha256_bytes(config_bytes),
            )
            partial = output_dir / "run_manifest.partial.json"
            self.assertEqual(json.loads(partial.read_text())["state"], "running")
            with self.assertRaises(FileExistsError):
                module.initialize_run(
                    output_dir=output_dir,
                    config_bytes=config_bytes,
                    config_sha256=module.sha256_bytes(config_bytes),
                )

            checkpoint = output_dir / "scf.chk"
            log_file = output_dir / "pyscf.log"
            checkpoint.write_bytes(b"checkpoint")
            log_file.write_text("calculation log\n", encoding="utf-8")
            stream = types.SimpleNamespace(flush=lambda: None)
            module.finalize_accepted_run(
                output_dir=output_dir,
                manifest=manifest,
                results={"energy_hartree": -1.1},
                artifacts=(checkpoint, log_file),
                log_path=log_file,
                mol=types.SimpleNamespace(stdout=stream),
                mf=types.SimpleNamespace(stdout=stream),
                pyscf_version="test-version",
                elapsed_seconds=1.25,
            )

            accepted_path = output_dir / "run_manifest.json"
            accepted = json.loads(accepted_path.read_text())
            self.assertEqual(accepted["state"], "accepted")
            self.assertEqual(accepted["pyscf_version"], "test-version")
            self.assertEqual(accepted["elapsed_seconds"], 1.25)
            self.assertIn("results.json", accepted["artifacts"])
            self.assertRegex(
                accepted["artifacts"]["results.json"]["sha256"], r"^[0-9a-f]{64}$"
            )
            self.assertFalse(partial.exists())

            failed_output = Path(temporary) / "failed-run"
            failed_manifest = module.initialize_run(
                output_dir=failed_output,
                config_bytes=config_bytes,
                config_sha256=module.sha256_bytes(config_bytes),
            )
            module.finalize_failed_run(
                output_dir=failed_output,
                manifest=failed_manifest,
                error=RuntimeError("SCF did not converge"),
                pyscf_version="test-version",
                elapsed_seconds=2.5,
            )
            failed = json.loads(
                (failed_output / "run_manifest.partial.json").read_text()
            )
            self.assertEqual(failed["state"], "failed")
            self.assertEqual(failed["failure"]["type"], "RuntimeError")
            self.assertEqual(failed["elapsed_seconds"], 2.5)
            self.assertFalse((failed_output / "run_manifest.json").exists())
            self.assertFalse((failed_output / "results.json").exists())

    def test_buffered_log_is_flushed_before_hash_and_rechecked_after_publish(self) -> None:
        module = load_runner()

        class BufferedLog:
            def __init__(self, path: Path):
                self.path = path
                self.pending = bytearray()
                self.flush_count = 0

            def write(self, text: str) -> None:
                self.pending.extend(text.encode("utf-8"))

            def flush(self) -> None:
                self.flush_count += 1
                if self.pending:
                    with self.path.open("ab") as handle:
                        handle.write(self.pending)
                    self.pending.clear()

        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "run"
            config_bytes = b'{"unit":"Angstrom"}\n'
            manifest = module.initialize_run(
                output_dir,
                config_bytes,
                module.sha256_bytes(config_bytes),
            )
            checkpoint = output_dir / "scf.chk"
            checkpoint.write_bytes(b"checkpoint")
            log_file = output_dir / "pyscf.log"
            stream = BufferedLog(log_file)
            stream.write("buffered PySCF log\n")
            mol = types.SimpleNamespace(stdout=stream)
            mf = types.SimpleNamespace(stdout=stream)

            accepted = module.finalize_accepted_run(
                output_dir=output_dir,
                manifest=manifest,
                results={"energy_hartree": -1.1},
                artifacts=(checkpoint, log_file),
                log_path=log_file,
                mol=mol,
                mf=mf,
                pyscf_version="test-version",
                elapsed_seconds=1.25,
            )

            self.assertEqual(stream.flush_count, 2)
            self.assertEqual(log_file.read_bytes(), b"buffered PySCF log\n")
            self.assertEqual(
                accepted["artifacts"]["pyscf.log"]["sha256"],
                module.sha256_bytes(b"buffered PySCF log\n"),
            )

    def test_log_change_during_manifest_publish_revokes_acceptance(self) -> None:
        module = load_runner()
        config_bytes = b'{"unit":"Angstrom"}\n'
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "run"
            manifest = module.initialize_run(
                output_dir,
                config_bytes,
                module.sha256_bytes(config_bytes),
            )
            checkpoint = output_dir / "scf.chk"
            log_file = output_dir / "pyscf.log"
            checkpoint.write_bytes(b"checkpoint")
            log_file.write_bytes(b"initial log\n")
            stream = types.SimpleNamespace(flush=lambda: None)
            real_write = module._write_json_atomic

            def write_then_mutate(path, payload):
                real_write(path, payload)
                if Path(path).name == "run_manifest.json":
                    with log_file.open("ab") as handle:
                        handle.write(b"late log line\n")

            with mock.patch.object(module, "_write_json_atomic", write_then_mutate):
                with self.assertRaisesRegex(RuntimeError, "changed after"):
                    module.finalize_accepted_run(
                        output_dir=output_dir,
                        manifest=manifest,
                        results={"energy_hartree": -1.1},
                        artifacts=(checkpoint, log_file),
                        log_path=log_file,
                        mol=types.SimpleNamespace(stdout=stream),
                        mf=types.SimpleNamespace(stdout=stream),
                        pyscf_version="test-version",
                        elapsed_seconds=1.25,
                    )

            self.assertFalse((output_dir / "run_manifest.json").exists())
            self.assertFalse((output_dir / "results.json").exists())

    def test_stdout_flush_failure_cannot_publish_acceptance(self) -> None:
        module = load_runner()

        class BrokenStream:
            def flush(self) -> None:
                raise OSError("fixture flush failure")

        config_bytes = b'{"unit":"Angstrom"}\n'
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary) / "run"
            manifest = module.initialize_run(
                output_dir,
                config_bytes,
                module.sha256_bytes(config_bytes),
            )
            checkpoint = output_dir / "scf.chk"
            log_file = output_dir / "pyscf.log"
            checkpoint.write_bytes(b"checkpoint")
            log_file.write_bytes(b"initial log\n")

            with self.assertRaisesRegex(RuntimeError, "flush"):
                module.finalize_accepted_run(
                    output_dir=output_dir,
                    manifest=manifest,
                    results={"energy_hartree": -1.1},
                    artifacts=(checkpoint, log_file),
                    log_path=log_file,
                    mol=types.SimpleNamespace(stdout=BrokenStream()),
                    mf=types.SimpleNamespace(stdout=BrokenStream()),
                    pyscf_version="test-version",
                    elapsed_seconds=1.25,
                )

            self.assertFalse((output_dir / "run_manifest.json").exists())
            self.assertFalse((output_dir / "results.json").exists())

    def test_legacy_clis_fail_closed_before_optional_imports(self) -> None:
        legacy = sorted((PYSCF_SKILL / "tools").glob("*.py"))
        legacy = [path for path in legacy if path.name != "_legacy_guard.py"]
        legacy.append(PYSCF_SKILL / "scripts" / "dft_calculation.py")
        legacy.append(PYSCF_SKILL / "references" / "benzene-dft-tddft.py")
        self.assertTrue(legacy)

        for path in legacy:
            with self.subTest(path=path.name):
                completed = subprocess.run(
                    [sys.executable, str(path), "--help"],
                    capture_output=True,
                    check=False,
                    text=True,
                )
                combined = (completed.stdout + completed.stderr).lower()
                self.assertNotEqual(completed.returncode, 0)
                self.assertIn("quarantined", combined)
                self.assertIn("run_safe_dft_tda.py", combined)

    def test_dependency_import_failure_is_recorded_not_accepted(self) -> None:
        module = load_runner()
        config = {
            "atom": "H 0 0 0; H 0 0 0.74",
            "unit": "Angstrom",
            "basis": "sto-3g",
            "charge": 0,
            "spin": 0,
            "xc": "pbe0",
            "grid_level": 3,
            "conv_tol": 1e-9,
            "max_cycle": 80,
            "nstates": 2,
        }
        real_import = __import__

        def block_pyscf(name, *args, **kwargs):
            if name == "pyscf":
                raise ModuleNotFoundError("test blocked PySCF")
            return real_import(name, *args, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            config_path = temporary_path / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            output_dir = temporary_path / "run"
            with mock.patch("builtins.__import__", side_effect=block_pyscf):
                with self.assertRaises(ModuleNotFoundError):
                    module.run(config_path, output_dir)

            failed = json.loads(
                (output_dir / "run_manifest.partial.json").read_text()
            )
            self.assertEqual(failed["state"], "failed")
            self.assertEqual(failed["pyscf_version"], "not_loaded")
            self.assertFalse((output_dir / "run_manifest.json").exists())
            self.assertFalse((output_dir / "results.json").exists())

    def test_nonphysical_tda_values_cannot_publish_results(self) -> None:
        module = load_runner()
        config = {
            "atom": "H 0 0 0; H 0 0 0.74",
            "unit": "Angstrom",
            "basis": "sto-3g",
            "charge": 0,
            "spin": 0,
            "xc": "pbe0",
            "grid_level": 3,
            "conv_tol": 1e-9,
            "max_cycle": 80,
            "nstates": 2,
        }

        def fake_pyscf(energies, strengths):
            fake_module = types.ModuleType("pyscf")
            fake_module.__version__ = "test-version"

            class FakeMolecule:
                nelectron = 2
                spin = 0
                natm = 2

                def __init__(self, **kwargs):
                    Path(kwargs["output"]).write_text(
                        "fake PySCF log\n", encoding="utf-8"
                    )

            class FakeTDA:
                def __init__(self):
                    self.nstates = 0
                    self.converged = [True, True]
                    self.e = energies

                def kernel(self):
                    return self.e

                def oscillator_strength(self):
                    return strengths

            class FakeRKS:
                def __init__(self, molecule):
                    self.mol = molecule
                    self.grids = types.SimpleNamespace(level=None)
                    self.converged = True
                    self.e_tot = -1.1
                    self.mo_energy = [-0.5, 0.2]
                    self.mo_coeff = [[1.0, 0.0], [0.0, 1.0]]
                    self.mo_occ = [2.0, 0.0]
                    self.chkfile = ""

                def kernel(self):
                    Path(self.chkfile).write_bytes(b"checkpoint")
                    return self.e_tot

                def stability(self, **kwargs):
                    return self.mo_coeff, self.mo_coeff, True, True

                def TDA(self):
                    return FakeTDA()

            fake_module.gto = types.SimpleNamespace(M=FakeMolecule)
            fake_module.dft = types.SimpleNamespace(RKS=FakeRKS)
            fake_module.lib = types.SimpleNamespace(
                param=types.SimpleNamespace(HARTREE2EV=27.211386245988)
            )
            return fake_module

        cases = (
            ([0.0, 0.3], [0.1, 0.2], "zero-energy"),
            ([-0.2, 0.3], [0.1, 0.2], "negative-energy"),
            ([0.2, 0.3], [-0.1, 0.2], "negative-strength"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            config_path = temporary_path / "config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            for energies, strengths, label in cases:
                output_dir = temporary_path / label
                with self.subTest(label=label):
                    with mock.patch.dict(
                        sys.modules,
                        {"pyscf": fake_pyscf(energies, strengths)},
                    ):
                        with self.assertRaises(ValueError):
                            module.run(config_path, output_dir)
                    failed = json.loads(
                        (output_dir / "run_manifest.partial.json").read_text()
                    )
                    self.assertEqual(failed["state"], "failed")
                    self.assertFalse((output_dir / "run_manifest.json").exists())
                    self.assertFalse((output_dir / "results.json").exists())

    def test_all_pyscf_python_sources_compile(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            cache = Path(temporary)
            for source in sorted(PYSCF_SKILL.rglob("*.py")):
                with self.subTest(source=source.relative_to(PYSCF_SKILL)):
                    py_compile.compile(
                        str(source),
                        cfile=str(cache / f"{source.stem}-{len(str(source))}.pyc"),
                        doraise=True,
                    )

    def test_skill_names_supported_runner_and_exact_acceptance_checks(self) -> None:
        text = (PYSCF_SKILL / "SKILL.md").read_text(encoding="utf-8")
        for required in (
            "scripts/run_safe_dft_tda.py",
            "len(td.converged) == td.nstates",
            "len(td.e) == td.nstates",
            "len(oscillator_strength) == td.nstates",
            "all(value > 0.0 for value in td.e)",
            "all(value >= 0.0 for value in oscillator_strength)",
            "run_manifest.partial.json",
            "run_manifest.json",
            "results.json",
        ):
            with self.subTest(required=required):
                self.assertIn(required, text)


if __name__ == "__main__":
    unittest.main()
