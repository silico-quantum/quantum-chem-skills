from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ScientificRegressionTests(unittest.TestCase):
    def test_dft_example_converts_hartree_to_kcal_once(self) -> None:
        source = (ROOT / "pyscf" / "scripts" / "dft_calculation.py").read_text(
            encoding="utf-8"
        )

        self.assertIn("(e - e_ref) * 627.509", source)
        self.assertNotIn("27.2114 * 627.509", source)

    def test_stage4_report_does_not_invent_method_provenance(self) -> None:
        module_path = ROOT / "tadf-screening" / "stage4_momap.py"
        spec = importlib.util.spec_from_file_location("stage4_momap", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "report.md"
            module.generate_report([], output)
            report = output.read_text(encoding="utf-8")

        self.assertNotIn("B3LYP/6-31G*", report)
        self.assertIn("verify software and calculation provenance", report)

    def test_excited_state_pes_uses_total_energy(self) -> None:
        guide = (
            ROOT
            / "pyscf"
            / "references"
            / "practice"
            / "2d-potential-energy-surface.md"
        ).read_text(encoding="utf-8")

        self.assertIn("e_s1 = e_gs + td.e[0]", guide)
        self.assertIn("excitation_energies_ev", guide)
        self.assertNotIn("C-C = 1.8000 Å", guide)

    def test_jk_integrals_are_labeled_coulomb_and_exchange(self) -> None:
        reference = (
            ROOT
            / "pyscf"
            / "references"
            / "theory"
            / "pyscf-api-reference.md"
        ).read_text(encoding="utf-8")

        self.assertIn("veff = mf.get_veff", reference)
        self.assertNotIn("vxc = mf.get_veff", reference)
        self.assertIn("# Coulomb matrix", reference)
        self.assertIn("# Exchange matrix", reference)
        self.assertNotIn("# Correlation\n", reference)

    def test_benzene_nto_summary_matches_logged_weights(self) -> None:
        examples = (ROOT / "examples" / "README.md").read_text(encoding="utf-8")

        self.assertIn("| S1 | 5.4951 | 225.7 | 0.0000 |", examples)
        self.assertIn("| S3 | 7.2776 | 170.4 | 0.5737 |", examples)
        self.assertIn("| S4 | 7.2781 | 170.4 | 0.5733 |", examples)
        self.assertIn("leading weights 0.4819 and 0.4805", examples)
        self.assertNotIn("single-pair transition", examples)

    def test_tddft_example_uses_callable_properties_and_one_based_nto_state(self) -> None:
        source = (ROOT / "pyscf" / "scripts" / "dft_calculation.py").read_text(
            encoding="utf-8"
        )

        self.assertIn("oscillator_strengths = td.oscillator_strength()", source)
        self.assertNotIn("td.oscillator_strength[", source)
        self.assertIn("weights, nto_coeff = td.get_nto(state=1)", source)
        self.assertNotIn("td.get_nto(state=0)", source)

    def test_tddft_references_match_pyscf_property_and_nto_interfaces(self) -> None:
        api_reference = (
            ROOT
            / "pyscf"
            / "references"
            / "theory"
            / "pyscf-api-reference.md"
        ).read_text(encoding="utf-8")
        advanced = (
            ROOT / "pyscf" / "references" / "theory" / "pyscf-advanced.md"
        ).read_text(encoding="utf-8")

        self.assertIn("oscillator_strengths = td.oscillator_strength()", api_reference)
        self.assertIn("transition_dipoles = td.transition_dipole()", api_reference)
        self.assertNotIn("td.oscillator_strength):", api_reference)
        self.assertNotIn("td.transition_dipole[", api_reference)
        self.assertIn("weights, nto_coeff = td.get_nto(state=1)", api_reference)
        self.assertNotIn("td.get_nto(state=0)", api_reference)
        self.assertIn("from pyscf import dft, scf", api_reference)
        self.assertNotIn("scf.RKS", api_reference)
        self.assertNotIn("scf.UKS", api_reference)

        self.assertIn("td.get_nto(state=n + 1)", advanced)
        self.assertNotIn("td.get_nto(state=n)", advanced)

    def test_emission_guides_use_state_specific_gradients_and_total_s1_energy(self) -> None:
        guide = (
            ROOT
            / "pyscf"
            / "references"
            / "practice"
            / "emission-spectrum-guide.md"
        ).read_text(encoding="utf-8")
        workflow = (
            ROOT
            / "pyscf"
            / "references"
            / "practice"
            / "emission-spectrum-workflow.md"
        ).read_text(encoding="utf-8")

        for document in (guide, workflow):
            self.assertIn("nuc_grad_method().as_scanner(state=1)", document)
            self.assertNotIn("optimize(td_s1, state=0)", document)

        self.assertIn("e_s1 = e_gs + td.e[0]", workflow)
        self.assertIn("s1_total_energies", workflow)
        self.assertNotIn("s1_energies.append(td.e[0]", workflow)

    def test_pyscf_examples_use_atomic_charges_and_valid_checkpoint_calls(self) -> None:
        source = (ROOT / "pyscf" / "scripts" / "dft_calculation.py").read_text(
            encoding="utf-8"
        )
        api_reference = (
            ROOT
            / "pyscf"
            / "references"
            / "theory"
            / "pyscf-api-reference.md"
        ).read_text(encoding="utf-8")

        self.assertIn("_, atomic_charges = mf.mulliken_pop", source)
        self.assertIn("charge = atomic_charges[i]", source)
        self.assertIn("mf.dump_chk(mf.chkfile)", source)
        self.assertNotIn("mf.dump_chk()", source)

        self.assertIn("mf.dump_chk(mf.chkfile)", api_reference)
        self.assertNotIn("mf.dump_chk()", api_reference)
        self.assertNotIn("mf.dump_chk(mf.chkfile, mf.mo_coeff", api_reference)
        self.assertIn("scf.chkfile.load_scf('calculation.chk')", api_reference)
        self.assertIn("scf_data['mo_coeff']", api_reference)

    def test_advanced_tddft_reference_avoids_nonexistent_density_helpers(self) -> None:
        advanced = (
            ROOT / "pyscf" / "references" / "theory" / "pyscf-advanced.md"
        ).read_text(encoding="utf-8")

        self.assertNotIn("td.generate_density", advanced)
        self.assertNotIn("td.get_rho", advanced)
        self.assertNotIn("iao.mulliken_pop", advanced)

    def test_rdkit_nci_uses_computed_charges_and_converts_dipole_units(self) -> None:
        source = (
            ROOT / "rdkit-chemistry" / "examples" / "nci_visualization.py"
        ).read_text(encoding="utf-8")

        self.assertNotIn("for atom in mol.GetAtoms():", source)
        self.assertIn("E_ANGSTROM_TO_DEBYE = 4.8032047", source)
        self.assertIn(
            "dipole_magnitude = np.linalg.norm(dipole) * E_ANGSTROM_TO_DEBYE",
            source,
        )

    def test_quantum_example_removes_placeholder_fukui_results_and_tracks_provenance(self) -> None:
        source = (
            ROOT / "rdkit-chemistry" / "examples" / "advanced_quantum_calc.py"
        ).read_text(encoding="utf-8")

        self.assertNotIn("f_plus = mulliken_charges", source)
        self.assertNotIn("Fukui Functions", source)
        self.assertIn("Method: DFT (B3LYP/6-31G)", source)
        self.assertIn('"benzene_HOMO.cube"', source)
        self.assertNotIn("DMAC-TRZ_HOMO.cube", source)
        self.assertIn("if not mf.converged:", source)
        self.assertIn("if not mf_plus.converged:", source)
        self.assertIn("if not mf_minus.converged:", source)


if __name__ == "__main__":
    unittest.main()
