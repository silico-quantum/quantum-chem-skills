from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

SKILLS = (
    "gaussian",
    "molecular-orbital-analysis",
    "molecular-sampler",
    "momap",
    "multiwfn",
    "pyscf",
    "rdkit-chemistry",
    "xtb-cluster-md",
    "xyzrender",
)

OPERATIONAL_SECTIONS = (
    "## Prerequisites",
    "## Input contract",
    "## Workflow",
    "## Validation and acceptance",
    "## Failure handling",
    "## Output and reporting",
    "## References",
)

DOMAIN_REQUIREMENTS = {
    "gaussian": (
        "charge and multiplicity",
        "normal termination",
        "imaginary frequencies",
        "checkpoint",
    ),
    "molecular-orbital-analysis": (
        "charge and spin",
        "orbital index",
        "cube",
        "isosurface",
    ),
    "molecular-sampler": (
        "atom count",
        "distance cutoff",
        "deterministic",
        "sampling_summary.txt",
    ),
    "momap": (
        "not_computed",
        "rate units",
        "spectrum",
        "gaussian",
    ),
    "multiwfn": (
        "wavefunction provenance",
        "electron count",
        "batch",
        "output file",
    ),
    "pyscf": (
        "charge and spin",
        "converged",
        "pyscf version",
        "units",
    ),
    "rdkit-chemistry": (
        "sanitization",
        "embedding return code",
        "random seed",
        "force-field parameters",
    ),
    "xtb-cluster-md": (
        "charge and unpaired electrons",
        "time step",
        "trajectory",
        "normal termination",
    ),
    "xyzrender": (
        "atom count",
        "input format",
        "non-empty",
        "xyzrender version",
    ),
}


class SkillCompletenessTests(unittest.TestCase):
    def assert_skill_complete(self, name: str) -> None:
        skill_file = ROOT / name / "SKILL.md"
        text = skill_file.read_text(encoding="utf-8")
        lowered = text.lower()

        for heading in OPERATIONAL_SECTIONS:
            with self.subTest(skill=name, heading=heading):
                self.assertTrue(
                    heading in text, f"{name}/SKILL.md is missing {heading}"
                )

        for phrase in DOMAIN_REQUIREMENTS[name]:
            with self.subTest(skill=name, phrase=phrase):
                self.assertTrue(
                    phrase in lowered,
                    f"{name}/SKILL.md is missing domain requirement: {phrase}",
                )

        metadata_file = ROOT / name / "agents" / "openai.yaml"
        with self.subTest(skill=name, metadata="agents/openai.yaml"):
            self.assertTrue(metadata_file.is_file())
        if not metadata_file.is_file():
            return

        metadata = metadata_file.read_text(encoding="utf-8")
        self.assertRegex(metadata, r'(?m)^interface:\s*$')
        self.assertRegex(metadata, r'(?m)^  display_name: "[^"]+"$')
        short_match = re.search(
            r'(?m)^  short_description: "([^"]+)"$', metadata
        )
        self.assertIsNotNone(short_match)
        if short_match:
            self.assertGreaterEqual(len(short_match.group(1)), 25)
            self.assertLessEqual(len(short_match.group(1)), 64)
        prompt_match = re.search(r'(?m)^  default_prompt: "([^"]+)"$', metadata)
        self.assertIsNotNone(prompt_match)
        if prompt_match:
            self.assertIn(f"${name}", prompt_match.group(1))

    def test_gaussian(self) -> None:
        self.assert_skill_complete("gaussian")

    def test_molecular_orbital_analysis(self) -> None:
        self.assert_skill_complete("molecular-orbital-analysis")

    def test_molecular_sampler(self) -> None:
        self.assert_skill_complete("molecular-sampler")

    def test_momap(self) -> None:
        self.assert_skill_complete("momap")

    def test_multiwfn(self) -> None:
        self.assert_skill_complete("multiwfn")

    def test_pyscf(self) -> None:
        self.assert_skill_complete("pyscf")

    def test_rdkit_chemistry(self) -> None:
        self.assert_skill_complete("rdkit-chemistry")

    def test_xtb_cluster_md(self) -> None:
        self.assert_skill_complete("xtb-cluster-md")

    def test_xyzrender(self) -> None:
        self.assert_skill_complete("xyzrender")


if __name__ == "__main__":
    unittest.main()
