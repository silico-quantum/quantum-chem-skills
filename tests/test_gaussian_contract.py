from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class GaussianContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.skill = (ROOT / "gaussian" / "SKILL.md").read_text(encoding="utf-8")
        cls.workflows = (ROOT / "gaussian" / "WORKFLOWS.md").read_text(
            encoding="utf-8"
        )

    def test_emission_chain_uses_distinct_checkpoint_stages(self) -> None:
        for checkpoint in (
            "%Chk=<s0-optfreq-run-id>.chk",
            "%OldChk=<s0-optfreq-run-id>.chk",
            "%Chk=<vertical-run-id>.chk",
            "%OldChk=<vertical-run-id>.chk",
            "%Chk=<s1-opt-run-id>.chk",
            "%OldChk=<s1-opt-run-id>.chk",
            "%Chk=<s1-freq-run-id>.chk",
            "%Chk=<emission-run-id>.chk",
        ):
            with self.subTest(checkpoint=checkpoint):
                self.assertIn(checkpoint, self.workflows)

    def test_emission_observable_and_unverified_minimum_are_explicit(self) -> None:
        self.assertIn("printed TD excitation energy", self.workflows)
        self.assertIn("Do not subtract raw SCF", self.workflows)
        self.assertIn("minimum_not_verified", self.workflows)
        self.assertIn("oscillator strength", self.skill)
        self.assertIn("state-continuity", self.skill)

    def test_checkpoint_conversion_and_scan_syntax_are_directional(self) -> None:
        self.assertIn("formchk molecule.chk molecule.fchk", self.workflows)
        self.assertIn(
            "unfchk molecule.fchk molecule-restored.chk", self.workflows
        )
        self.assertNotIn("unfchk molecule.fchk molecule.chk", self.workflows)
        self.assertIn("D 1 2 3 4 S 20 5.0", self.workflows)
        self.assertNotIn("D 1 2 3 4 0.0 S 20 5.0", self.workflows)


if __name__ == "__main__":
    unittest.main()
