from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "molecular-orbital-analysis"
    / "scripts"
    / "generate_orbital_cubes.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("generate_orbital_cubes", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MO = load_module()


def cube_text(
    values: list[float | str],
    *,
    atom_records: list[str] | None = None,
    grid: tuple[int, int, int] = (2, 2, 2),
) -> str:
    records = atom_records if atom_records is not None else ["8 0.0 0.0 0.0 0.0"]
    lines = [
        "orbital cube",
        "generated fixture",
        f"{len(records)} 0.0 0.0 0.0",
        f"{grid[0]} 0.2 0.0 0.0",
        f"{grid[1]} 0.0 0.2 0.0",
        f"{grid[2]} 0.0 0.0 0.2",
        *records,
        " ".join(str(value) for value in values),
    ]
    return "\n".join(lines) + "\n"


class FakeVector:
    def __init__(self, values: list[float]):
        self.values = tuple(values)
        self.shape = (len(values),)

    @property
    def flat(self):
        return iter(self.values)

    def __iter__(self):
        return iter(self.values)

    def __getitem__(self, index: int) -> float:
        return self.values[index]


class FakeMatrix:
    def __init__(self, rows: list[list[float]]):
        if not rows or not rows[0]:
            raise ValueError("fixture matrix must be non-empty")
        width = len(rows[0])
        if any(len(row) != width for row in rows):
            raise ValueError("fixture matrix must be rectangular")
        self.rows = tuple(tuple(row) for row in rows)
        self.shape = (len(rows), width)
        self.ndim = 2

    @property
    def flat(self):
        return (value for row in self.rows for value in row)

    def __getitem__(self, key):
        row_selector, column = key
        if row_selector != slice(None):
            raise TypeError("fixture supports only full-column slicing")
        return [row[column] for row in self.rows]


def spatial_channels(
    *,
    coefficients: list[list[float]] | None = None,
    energies: list[float] | None = None,
    occupations: list[float] | None = None,
):
    return {
        "spatial": (
            FakeMatrix(coefficients or [[1.0, 0.0], [0.0, 1.0]]),
            FakeVector(energies or [-0.5, 0.2]),
            FakeVector(occupations or [2.0, 0.0]),
        )
    }


def unrestricted_channels(
    *,
    alpha_occupations: list[float] | None = None,
    beta_occupations: list[float] | None = None,
):
    return {
        "alpha": (
            FakeMatrix([[1.0, 0.0], [0.0, 1.0]]),
            FakeVector([-0.5, 0.2]),
            FakeVector(
                [1.0, 0.0] if alpha_occupations is None else alpha_occupations
            ),
        ),
        "beta": (
            FakeMatrix([[1.0, 0.0], [0.0, 1.0]]),
            FakeVector([-0.4, 0.3]),
            FakeVector(
                [1.0, 0.0] if beta_occupations is None else beta_occupations
            ),
        ),
    }


class FakeMolecule:
    natm = 1
    nelectron = 2


class FakeCubeGenerator:
    def __init__(self, payloads: list[list[float | str]]):
        self.payloads = list(payloads)
        self.calls: list[Path] = []

    def orbital(self, _mol, path, _coefficients, **_kwargs):
        output = Path(path)
        self.calls.append(output)
        payload = self.payloads[len(self.calls) - 1]
        output.write_text(cube_text(payload), encoding="utf-8")


class FakeMeanField:
    def __init__(
        self,
        *,
        on_kernel=None,
        checkpoint_content: bytes | None = b"chk",
        stale_checkpoint: bool = False,
    ):
        self.converged = True
        self.mo_coeff = FakeMatrix([[1.0, 0.0], [0.0, 1.0]])
        self.mo_energy = FakeVector([-0.5, 0.2])
        self.mo_occ = FakeVector([2.0, 0.0])
        self.max_cycle = None
        self.chkfile = None
        self.on_kernel = on_kernel
        self.checkpoint_content = checkpoint_content
        self.stale_checkpoint = stale_checkpoint

    def kernel(self) -> float:
        if self.on_kernel is not None:
            self.on_kernel()
        if self.checkpoint_content is not None and self.chkfile is not None:
            checkpoint = Path(self.chkfile)
            checkpoint.write_bytes(self.checkpoint_content)
            if self.stale_checkpoint:
                os.utime(checkpoint, ns=(1, 1))
        return -75.0


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run_args(xyz: Path, output_dir: Path):
    return SimpleNamespace(
        xyz=xyz,
        output_dir=output_dir,
        unit="Angstrom",
        charge=0,
        spin=0,
        basis="fixture-basis",
        method="rks",
        xc="fixture-xc",
        max_cycle=20,
        grid_points=2,
        margin_bohr=3.0,
        orbital=[],
    )


def fake_pyscf_modules(
    generator,
    *,
    log_content: bytes | None = b"PySCF fixture log\n",
    stale_log: bool = False,
    captured: dict | None = None,
):
    fake_pyscf = types.ModuleType("pyscf")
    fake_pyscf.__version__ = "fixture-version"
    fake_gto = types.ModuleType("pyscf.gto")

    def make_molecule(**kwargs):
        if captured is not None:
            captured["atom"] = kwargs["atom"]
        if log_content is not None:
            log_path = Path(kwargs["output"])
            log_path.write_bytes(log_content)
            if stale_log:
                os.utime(log_path, ns=(1, 1))
        return FakeMolecule()

    fake_gto.M = make_molecule
    fake_tools = types.ModuleType("pyscf.tools")
    fake_tools.cubegen = generator
    fake_pyscf.gto = fake_gto
    fake_pyscf.tools = fake_tools
    return fake_pyscf, fake_gto, fake_tools


class CubeValidationTests(unittest.TestCase):
    def test_complete_cube_records_exact_payload_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "complete.cube"
            cube.write_text(cube_text([0.1] * 8), encoding="utf-8")

            record = MO.validate_cube_file(cube, expected_atoms=1)

            self.assertEqual(record["grid"], [2, 2, 2])
            self.assertEqual(record["value_count"], 8)
            self.assertEqual(record["value_min"], 0.1)
            self.assertEqual(record["value_max"], 0.1)

    def test_truncated_cube_payload_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "truncated.cube"
            cube.write_text(cube_text([0.1] * 7), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "payload"):
                MO.validate_cube_file(cube, expected_atoms=1)

    def test_non_finite_cube_payload_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "nan.cube"
            cube.write_text(
                cube_text([0.1, 0.2, 0.3, "nan", 0.5, 0.6, 0.7, 0.8]),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "non-finite"):
                MO.validate_cube_file(cube, expected_atoms=1)

    def test_missing_cube_atom_record_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / "missing-atom.cube"
            cube.write_text(
                "\n".join(
                    [
                        "orbital cube",
                        "generated fixture",
                        "1 0.0 0.0 0.0",
                        "1 0.2 0.0 0.0",
                        "1 0.0 0.2 0.0",
                        "1 0.0 0.0 0.2",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "atom record"):
                MO.validate_cube_file(cube, expected_atoms=1)


class ElectronicResultValidationTests(unittest.TestCase):
    def test_finite_consistent_arrays_are_accepted(self) -> None:
        summary = MO.validate_electronic_result(
            -75.0, spatial_channels(), electron_count=2
        )
        self.assertEqual(summary, {"spatial": {"ao_count": 2, "orbital_count": 2}})

    def test_non_finite_total_energy_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "total energy"):
            MO.validate_electronic_result(
                math.nan, spatial_channels(), electron_count=2
            )

    def test_non_finite_orbital_array_is_rejected(self) -> None:
        channels = spatial_channels(coefficients=[[1.0, math.inf], [0.0, 1.0]])
        with self.assertRaisesRegex(ValueError, "coefficients"):
            MO.validate_electronic_result(-75.0, channels, electron_count=2)

    def test_mismatched_orbital_shapes_are_rejected(self) -> None:
        channels = spatial_channels(energies=[-0.5])
        with self.assertRaisesRegex(ValueError, "shape"):
            MO.validate_electronic_result(-75.0, channels, electron_count=2)

    def test_negative_occupation_is_rejected_for_each_orbital_layout(self) -> None:
        cases = (
            (spatial_channels(occupations=[2.0, -0.1]), "spatial"),
            (
                unrestricted_channels(alpha_occupations=[1.0, -0.1]),
                "alpha",
            ),
            (
                unrestricted_channels(beta_occupations=[1.0, -0.1]),
                "beta",
            ),
        )
        for channels, expected_channel in cases:
            with self.subTest(channel=expected_channel):
                with self.assertRaisesRegex(ValueError, expected_channel):
                    MO.validate_electronic_result(
                        -75.0, channels, electron_count=2
                    )

    def test_occupation_sum_must_match_molecular_electron_count(self) -> None:
        cases = (
            spatial_channels(occupations=[1.0, 0.0]),
            unrestricted_channels(
                alpha_occupations=[1.0, 0.0],
                beta_occupations=[0.0, 0.0],
            ),
        )
        for channels in cases:
            with self.subTest(channels=tuple(channels)):
                with self.assertRaisesRegex(ValueError, "electron count"):
                    MO.validate_electronic_result(
                        -75.0, channels, electron_count=2
                    )

    def test_occupation_sum_allows_only_small_numerical_roundoff(self) -> None:
        MO.validate_electronic_result(
            -75.0,
            spatial_channels(occupations=[1.999999995, 0.0]),
            electron_count=2,
        )
        with self.assertRaisesRegex(ValueError, "electron count"):
            MO.validate_electronic_result(
                -75.0,
                spatial_channels(occupations=[1.99999, 0.0]),
                electron_count=2,
            )


class CubePublicationTests(unittest.TestCase):
    def test_all_requests_are_preflighted_before_cube_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            generator = FakeCubeGenerator([[0.1] * 8])

            with self.assertRaisesRegex(IndexError, "exceeds"):
                MO.generate_cube_artifacts(
                    FakeMolecule(),
                    generator,
                    spatial_channels(),
                    [("spatial", 1), ("spatial", 3)],
                    output_dir,
                    grid_points=2,
                    margin_bohr=3.0,
                )

            self.assertEqual(generator.calls, [])
            self.assertEqual(list(output_dir.iterdir()), [])

    def test_failed_batch_publishes_no_final_cube_and_inventories_partials(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            generator = FakeCubeGenerator([[0.1] * 8, [0.2] * 7])

            with self.assertRaisesRegex(ValueError, "payload"):
                MO.generate_cube_artifacts(
                    FakeMolecule(),
                    generator,
                    spatial_channels(),
                    [("spatial", 1), ("spatial", 2)],
                    output_dir,
                    grid_points=2,
                    margin_bohr=3.0,
                )

            self.assertEqual(list(output_dir.glob("*.cube")), [])
            partials = sorted(output_dir.glob("*.partial"))
            self.assertEqual(len(partials), 2)
            inventory = MO.inventory_artifacts(output_dir)
            self.assertEqual(
                {record["status"] for record in inventory}, {"rejected_partial"}
            )
            self.assertTrue(all(record["sha256"] for record in inventory))

    def test_successful_batch_atomically_publishes_validated_cubes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            generator = FakeCubeGenerator([[0.1] * 8, [0.2] * 8])

            records = MO.generate_cube_artifacts(
                FakeMolecule(),
                generator,
                spatial_channels(),
                [("spatial", 1), ("spatial", 2)],
                output_dir,
                grid_points=2,
                margin_bohr=3.0,
            )

            self.assertEqual(len(records), 2)
            self.assertEqual(list(output_dir.glob("*.partial")), [])
            self.assertEqual(len(list(output_dir.glob("*.cube"))), 2)
            self.assertTrue(all(record["value_count"] == 8 for record in records))

    def test_failed_main_manifest_lists_every_rejected_partial(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            xyz = root / "molecule.xyz"
            xyz.write_text("1\nfixture\nO 0.0 0.0 0.0\n", encoding="utf-8")
            output_dir = root / "run"
            generator = FakeCubeGenerator([[0.1] * 8, [0.2] * 7])
            args = SimpleNamespace(
                xyz=xyz,
                output_dir=output_dir,
                unit="Angstrom",
                charge=0,
                spin=0,
                basis="fixture-basis",
                method="rks",
                xc="fixture-xc",
                max_cycle=20,
                grid_points=2,
                margin_bohr=3.0,
                orbital=[("spatial", 1), ("spatial", 2)],
            )

            fake_pyscf, fake_gto, fake_tools = fake_pyscf_modules(generator)

            with (
                mock.patch.object(MO, "parse_args", return_value=args),
                mock.patch.object(
                    MO, "build_mean_field", return_value=FakeMeanField()
                ),
                mock.patch.dict(
                    sys.modules,
                    {
                        "pyscf": fake_pyscf,
                        "pyscf.gto": fake_gto,
                        "pyscf.tools": fake_tools,
                    },
                ),
            ):
                with self.assertRaisesRegex(ValueError, "payload"):
                    MO.main()

            manifest = json.loads((output_dir / "run.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["failure_phase"], "generate_cube_artifacts")
            rejected = [
                record
                for record in manifest["artifacts"]
                if record["status"] == "rejected_partial"
            ]
            self.assertEqual(len(rejected), 2)
            self.assertEqual(list(output_dir.glob("*.cube")), [])


class RunSnapshotAndArtifactTests(unittest.TestCase):
    def run_with_fakes(
        self,
        args,
        mean_field,
        *,
        log_content: bytes | None = b"PySCF fixture log\n",
        stale_log: bool = False,
        captured: dict | None = None,
    ):
        generator = FakeCubeGenerator([])
        fake_pyscf, fake_gto, fake_tools = fake_pyscf_modules(
            generator,
            log_content=log_content,
            stale_log=stale_log,
            captured=captured,
        )
        with (
            mock.patch.object(MO, "parse_args", return_value=args),
            mock.patch.object(MO, "build_mean_field", return_value=mean_field),
            mock.patch.dict(
                sys.modules,
                {
                    "pyscf": fake_pyscf,
                    "pyscf.gto": fake_gto,
                    "pyscf.tools": fake_tools,
                },
            ),
        ):
            return MO.main()

    def test_dft_method_without_explicit_xc_fails_before_output_publication(self) -> None:
        for method in ("rks", "uks"):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "molecule.xyz"
                source.write_text(
                    "1\nfixture\nO 0.0 0.0 0.0\n", encoding="utf-8"
                )
                args = run_args(source, root / "run")
                args.method = method
                args.spin = 0 if method == "rks" else 2
                args.xc = None

                with mock.patch.object(MO, "parse_args", return_value=args):
                    with self.assertRaisesRegex(ValueError, "explicit --xc"):
                        MO.main()

                self.assertFalse(args.output_dir.exists())

    def test_xc_has_no_implicit_parser_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            argv = [
                str(SCRIPT),
                str(root / "molecule.xyz"),
                "--output-dir",
                str(root / "run"),
                "--unit",
                "Angstrom",
                "--charge",
                "0",
                "--spin",
                "0",
                "--basis",
                "sto-3g",
                "--method",
                "rks",
            ]
            with mock.patch.object(sys, "argv", argv):
                args = MO.parse_args()

            self.assertIsNone(args.xc)

    def test_source_rewrite_during_scf_cannot_change_snapshot_or_manifest_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "molecule.xyz"
            original = b"1\noriginal oxygen\nO 0.0 0.0 0.0\n"
            replacement = b"1\nreplacement hydrogen\nH 9.0 9.0 9.0\n"
            source.write_bytes(original)
            output_dir = root / "run"
            captured: dict = {}
            mean_field = FakeMeanField(on_kernel=lambda: source.write_bytes(replacement))

            result = self.run_with_fakes(
                run_args(source, output_dir), mean_field, captured=captured
            )

            self.assertEqual(result, 0)
            self.assertEqual(source.read_bytes(), replacement)
            snapshot = output_dir / "input.xyz"
            self.assertEqual(snapshot.read_bytes(), original)
            self.assertFalse((output_dir / "input.xyz.partial").exists())
            self.assertEqual(captured["atom"][0][0], "O")

            manifest = json.loads((output_dir / "run.json").read_text())
            self.assertEqual(manifest["status"], "accepted")
            self.assertEqual(manifest["input"]["path"], "input.xyz")
            self.assertEqual(manifest["input"]["source_path"], str(source))
            self.assertEqual(manifest["input"]["sha256"], sha256_bytes(original))
            self.assertEqual(
                manifest["runtime_artifacts"]["checkpoint"]["sha256"],
                sha256_bytes(b"chk"),
            )
            self.assertEqual(
                manifest["runtime_artifacts"]["pyscf_log"]["sha256"],
                sha256_bytes(b"PySCF fixture log\n"),
            )

    def test_missing_checkpoint_rejects_converged_scf(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "molecule.xyz"
            source.write_text("1\nfixture\nO 0.0 0.0 0.0\n", encoding="utf-8")
            output_dir = root / "run"

            with self.assertRaisesRegex(RuntimeError, "checkpoint"):
                self.run_with_fakes(
                    run_args(source, output_dir),
                    FakeMeanField(checkpoint_content=None),
                )

            manifest = json.loads((output_dir / "run.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["failure_phase"], "validate_runtime_artifacts")
            self.assertEqual(
                manifest["runtime_artifacts"]["checkpoint"]["status"], "invalid"
            )
            self.assertEqual(
                manifest["runtime_artifacts"]["pyscf_log"]["status"], "valid"
            )

    def test_missing_log_rejects_converged_scf(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "molecule.xyz"
            source.write_text("1\nfixture\nO 0.0 0.0 0.0\n", encoding="utf-8")
            output_dir = root / "run"

            with self.assertRaisesRegex(RuntimeError, "PySCF log"):
                self.run_with_fakes(
                    run_args(source, output_dir),
                    FakeMeanField(),
                    log_content=None,
                )

            manifest = json.loads((output_dir / "run.json").read_text())
            self.assertEqual(manifest["status"], "failed")
            self.assertEqual(manifest["failure_phase"], "validate_runtime_artifacts")
            self.assertEqual(
                manifest["runtime_artifacts"]["checkpoint"]["status"], "valid"
            )
            self.assertEqual(
                manifest["runtime_artifacts"]["pyscf_log"]["status"], "invalid"
            )

    def test_empty_or_stale_runtime_artifact_is_rejected(self) -> None:
        cases = (
            ("empty log", FakeMeanField(), b"", False, "non-empty"),
            ("stale log", FakeMeanField(), b"log\n", True, "fresh"),
            (
                "stale checkpoint",
                FakeMeanField(stale_checkpoint=True),
                b"log\n",
                False,
                "fresh",
            ),
        )
        for label, mean_field, log_content, stale_log, message in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source = root / "molecule.xyz"
                source.write_text(
                    "1\nfixture\nO 0.0 0.0 0.0\n", encoding="utf-8"
                )
                with self.assertRaisesRegex(RuntimeError, message):
                    self.run_with_fakes(
                        run_args(source, root / "run"),
                        mean_field,
                        log_content=log_content,
                        stale_log=stale_log,
                    )


if __name__ == "__main__":
    unittest.main()
