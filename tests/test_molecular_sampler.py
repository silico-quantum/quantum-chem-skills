from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "molecular-sampler" / "molecular_sampler.py"


def load_sampler():
    spec = importlib.util.spec_from_file_location("molecular_sampler", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load molecular_sampler.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MolecularSamplerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.sampler = load_sampler()

    def test_xyz_atom_count_is_strict(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xyz = Path(tmp) / "bad.xyz"
            xyz.write_text(
                "3\ntruncated\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8"
            )

            with self.assertRaises(ValueError):
                self.sampler.parse_xyz_file(xyz)

    def test_gaussian_bohr_geometry_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            gjf = Path(tmp) / "bohr.gjf"
            gjf.write_text(
                "#p hf/sto-3g units=bohr\n\n"
                "Bohr fixture\n\n"
                "0 1\n"
                "H 0 0 0\n"
                "H 0 0 1.4\n\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "bohr/atomic units"):
                self.sampler.parse_gjf_file(gjf)

    def test_gaussian_link1_input_is_rejected_instead_of_taking_first_geometry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            gjf = Path(tmp) / "multi_job.gjf"
            gjf.write_text(
                "%chk=first.chk\n#p hf/sto-3g\n\nfirst\n\n0 1\nH 0 0 0\n\n"
                "--Link1--\n"
                "%chk=second.chk\n#p hf/sto-3g\n\nsecond\n\n0 1\nH 1 0 0\n\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "exactly one Gaussian job"):
                self.sampler.parse_gjf_file(gjf)

    def test_charged_oniom_geometry_uses_trailing_layer_marker(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            gjf = Path(tmp) / "charged.gjf"
            gjf.write_text(
                "%chk=charged.chk\n"
                "# oniom(test:test)\n\n"
                "Charged ONIOM fixture\n\n"
                "-1 2 -1 2 0 1\n"
                "C-CT--0.100 0 0.000 0.000 0.000 H\n"
                "H-HC-0.100 0 0.000 0.000 1.090 L\n\n"
                "B 1 2 F\n",
                encoding="utf-8",
            )

            atoms = self.sampler.parse_gjf_file(gjf)

            self.assertEqual([atom["element"] for atom in atoms], ["C", "H"])
            self.assertEqual([atom["layer"] for atom in atoms], ["H", "L"])

    def test_fragment_detection_retains_small_components_and_all_atoms(self) -> None:
        atoms = [
            {"element": "H", "layer": "L", "x": 0.0, "y": 0.0, "z": 0.0, "source_index": 0},
            {"element": "H", "layer": "L", "x": 0.0, "y": 0.0, "z": 0.74, "source_index": 1},
            {"element": "H", "layer": "L", "x": 5.0, "y": 0.0, "z": 0.0, "source_index": 2},
            {"element": "H", "layer": "L", "x": 5.0, "y": 0.0, "z": 0.74, "source_index": 3},
        ]

        molecules = self.sampler.identify_molecules(atoms, layer_filter="all")

        self.assertEqual(len(molecules), 2)
        retained = sorted(
            atom["source_index"]
            for molecule in molecules
            for atom in molecule["atoms"]
        )
        self.assertEqual(retained, [0, 1, 2, 3])

    def test_sampling_is_unique_and_never_writes_undersized_complexes(self) -> None:
        molecules = [
            {
                "id": index,
                "atoms": [],
                "elements": {},
                "count": 0,
                "center": (float(index), 0.0, 0.0),
            }
            for index in range(4)
        ]

        samples = self.sampler.sample_molecules(molecules, n_samples=20)

        for size, name in (
            (2, "dimers"),
            (3, "trimers"),
            (4, "tetramers"),
            (5, "pentamers"),
        ):
            member_sets = [
                tuple(sorted(molecule["id"] for molecule in sample))
                for sample in samples[name]
            ]
            self.assertTrue(all(len(sample) == size for sample in samples[name]))
            self.assertEqual(len(member_sets), len(set(member_sets)))

        self.assertEqual(samples["pentamers"], [])

    def test_output_directory_must_be_new(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "samples"
            output.mkdir()
            (output / "stale.xyz").write_text("stale\n", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                self.sampler.create_output_directory(output)

    def test_xyz_cli_requires_an_explicit_unit_declaration(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xyz = Path(tmp) / "hydrogen.xyz"
            xyz.write_text(
                "2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "--xyz-units angstrom"):
                self.sampler.main(
                    [
                        str(xyz),
                        "--expected-fragments",
                        "1",
                        "--output-dir",
                        str(Path(tmp) / "samples"),
                    ]
                )

    def test_expected_fragment_count_fails_before_output_creation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xyz = Path(tmp) / "hydrogen.xyz"
            xyz.write_text(
                "2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8"
            )
            output = Path(tmp) / "samples"

            with self.assertRaisesRegex(ValueError, "Fragment-count mismatch"):
                self.sampler.main(
                    [
                        str(xyz),
                        "--xyz-units",
                        "angstrom",
                        "--expected-fragments",
                        "2",
                        "--output-dir",
                        str(output),
                    ]
                )
            self.assertFalse(output.exists())

    def test_expected_fragment_count_is_required_for_a_validated_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xyz = Path(tmp) / "hydrogen.xyz"
            xyz.write_text(
                "2\nhydrogen\nH 0 0 0\nH 0 0 0.74\n", encoding="utf-8"
            )
            output = Path(tmp) / "samples"

            with self.assertRaises(SystemExit):
                self.sampler.main(
                    [
                        str(xyz),
                        "--xyz-units",
                        "angstrom",
                        "--output-dir",
                        str(output),
                    ]
                )
            self.assertFalse(output.exists())

    def test_end_to_end_run_writes_a_validated_unique_dimer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xyz = Path(tmp) / "two_hydrogen_molecules.xyz"
            xyz.write_text(
                "4\ntwo H2 fragments; angstrom\n"
                "H 0 0 0\nH 0 0 0.74\n"
                "H 5 0 0\nH 5 0 0.74\n",
                encoding="utf-8",
            )
            output = Path(tmp) / "samples"

            exit_code = self.sampler.main(
                [
                    str(xyz),
                    "--xyz-units",
                    "angstrom",
                    "--layer",
                    "all",
                    "--expected-fragments",
                    "2",
                    "--output-dir",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            manifest = json.loads(
                (output / "sampling_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["input"]["selected_atom_count"], 4)
            self.assertEqual(manifest["status"], "internal_validation_passed")
            self.assertEqual(
                manifest["scientific_status"],
                "pending_independent_fragment_review",
            )
            self.assertEqual(
                manifest["fragment_summary"]["detected_before_size_filter"], 2
            )
            self.assertEqual(len(manifest["fragments"]), 2)
            dimers = [
                record for record in manifest["samples"] if record["type"] == "dimer"
            ]
            self.assertEqual(len(dimers), 1)
            self.assertEqual(dimers[0]["member_ids"], [0, 1])
            self.assertTrue((output / dimers[0]["path"]).is_file())
            self.assertFalse((output / "RUN_INCOMPLETE").exists())


if __name__ == "__main__":
    unittest.main()
