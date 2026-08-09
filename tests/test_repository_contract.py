from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "scripts" / "validate_repository.py"


def run_validator(path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(path)],
        check=False,
        capture_output=True,
        text=True,
    )


class RepositoryContractTests(unittest.TestCase):
    def test_saved_logs_do_not_expose_local_build_identity(self) -> None:
        for name in ("xtb_md.log", "xtb_run.log"):
            log = (
                ROOT / "examples" / "xtb-cluster-md" / name
            ).read_text(encoding="utf-8")
            with self.subTest(name=name):
                self.assertNotIn("runner@Mac-", log)

        pyscf_log = (ROOT / "examples" / "pyscf" / "benzene_output.log").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("MolBotdeMac-mini-2.local", pyscf_log)

    def test_tadf_adapter_has_agent_neutral_momap_tool_discovery(self) -> None:
        source = (ROOT / "tadf-screening" / "stage4_momap.py").read_text(
            encoding="utf-8"
        )

        self.assertIn("MOMAP_TOOLS_DIR", source)
        self.assertNotIn("~/.openclaw/skills", source)

    def test_ci_runs_repository_contract_checks(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "validate.yml").read_text(
            encoding="utf-8"
        )

        self.assertIn("python3 -m unittest discover -s tests -v", workflow)
        self.assertIn("python3 scripts/validate_repository.py .", workflow)

    def test_node_metadata_does_not_model_the_python_runtime(self) -> None:
        package_path = ROOT / "molecular-sampler" / "package.json"
        package = json.loads(package_path.read_text(encoding="utf-8"))

        self.assertNotIn("python", package.get("dependencies", {}))
        self.assertFalse(str(package.get("main", "")).endswith(".py"))

    def test_readme_documents_supported_agent_install_locations(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        expected_fragments = (
            "OpenAI Codex",
            "Claude Code",
            "OpenClaw",
            "GitHub Copilot",
            ".agents/skills",
            "~/.agents/skills",
            ".claude/skills",
            "~/.claude/skills",
            "<workspace>/skills",
            "~/.openclaw/skills",
            ".github/skills",
            "~/.copilot/skills",
        )
        for fragment in expected_fragments:
            with self.subTest(fragment=fragment):
                self.assertIn(fragment, readme)

    def test_openclaw_workspace_install_avoids_external_symlinks(self) -> None:
        install = (ROOT / "INSTALL.md").read_text(encoding="utf-8")
        openclaw_section = install.split("### OpenClaw", 1)[1].split(
            "### GitHub Copilot", 1
        )[0]

        self.assertIn("openclaw skills install ./pyscf --as pyscf", openclaw_section)
        self.assertIn("cp -R pyscf /path/to/workspace/skills/", openclaw_section)
        self.assertNotIn(
            'ln -s "$PWD/pyscf" /path/to/workspace/skills/pyscf',
            openclaw_section,
        )

        self.assertIn("<state-dir>/skills", openclaw_section)
        self.assertIn("default state", openclaw_section)

    def test_copilot_cloud_install_uses_repository_content(self) -> None:
        install = (ROOT / "INSTALL.md").read_text(encoding="utf-8")
        copilot_section = install.split("### GitHub Copilot", 1)[1].split(
            "## 4. Install Scientific Dependencies Selectively", 1
        )[0]

        self.assertIn("local IDE or CLI sessions", copilot_section)
        self.assertIn("cloud-hosted Copilot", copilot_section)
        self.assertIn("cp -R pyscf /path/to/project/.github/skills/", copilot_section)

    def test_valid_skill_package_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "example-skill"
            skill.mkdir()
            (skill / "SKILL.md").write_text(
                textwrap.dedent(
                    """\
                    ---
                    name: example-skill
                    description: Use when a user needs an example skill.
                    license: MIT
                    compatibility: Requires Python 3.10 or newer.
                    ---

                    # Example Skill

                    Follow the documented workflow.
                    """
                ),
                encoding="utf-8",
            )
            (root / "README.md").write_text("# Example\n", encoding="utf-8")
            (root / "LICENSE").write_text(
                "MIT License\nPermission is hereby granted, free of charge, to any person obtaining a copy\n"
                "THE SOFTWARE IS PROVIDED \"AS IS\"\n",
                encoding="utf-8",
            )

            result = run_validator(root)

            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_invalid_skill_package_reports_actionable_errors(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "wrong-directory"
            skill.mkdir()
            cjk_text = "".join(chr(code) for code in (0x4E2D, 0x6587))
            (skill / "SKILL.md").write_text(
                "---\nname: different-name\ndescription: Vague description\n"
                f"version: 1.0\n---\n{cjk_text}\n",
                encoding="utf-8",
            )
            (root / "README.md").write_text(
                "[missing](missing.md)\n", encoding="utf-8"
            )
            (root / "LICENSE").write_text("MIT License\n", encoding="utf-8")

            result = run_validator(root)

            self.assertNotEqual(result.returncode, 0)
            output = result.stdout + result.stderr
            self.assertIn("name must match its parent directory", output)
            self.assertIn("unsupported frontmatter field", output)
            self.assertIn("description must start with 'Use when'", output)
            self.assertIn("contains CJK text", output)
            self.assertIn("broken local link", output)
            self.assertIn("LICENSE is not a complete MIT license", output)

    def test_public_python_source_must_be_english(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "example-skill"
            skill.mkdir()
            (skill / "SKILL.md").write_text(
                "---\n"
                "name: example-skill\n"
                "description: Use when testing public source language.\n"
                "license: MIT\n"
                "---\n",
                encoding="utf-8",
            )
            cjk_output = "".join(
                chr(code) for code in (0x4E2D, 0x6587, 0x8F93, 0x51FA)
            )
            (root / "example.py").write_text(
                f'print("{cjk_output}")\n', encoding="utf-8"
            )
            (root / "LICENSE").write_text(
                "MIT License\nPermission is hereby granted, free of charge\n"
                'THE SOFTWARE IS PROVIDED "AS IS"\n',
                encoding="utf-8",
            )

            result = run_validator(root)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "example.py: contains CJK text; public source must be English",
                result.stdout + result.stderr,
            )

    def test_skill_reference_must_be_linked_from_skill_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "example-skill"
            references = skill / "references"
            references.mkdir(parents=True)
            (skill / "SKILL.md").write_text(
                "---\n"
                "name: example-skill\n"
                "description: Use when testing reference discovery.\n"
                "license: MIT\n"
                "---\n",
                encoding="utf-8",
            )
            (references / "details.md").write_text(
                "# Details\n", encoding="utf-8"
            )
            (root / "LICENSE").write_text(
                "MIT License\nPermission is hereby granted, free of charge\n"
                'THE SOFTWARE IS PROVIDED "AS IS"\n',
                encoding="utf-8",
            )

            result = run_validator(root)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "example-skill/references/details.md: not linked from example-skill/SKILL.md",
                result.stdout + result.stderr,
            )

    def test_skill_entrypoint_must_stay_within_context_budget(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "example-skill"
            skill.mkdir()
            body = "\n".join(f"Guidance line {index}." for index in range(501))
            (skill / "SKILL.md").write_text(
                "---\n"
                "name: example-skill\n"
                "description: Use when testing progressive disclosure.\n"
                "license: MIT\n"
                "---\n\n"
                f"{body}\n",
                encoding="utf-8",
            )
            (root / "LICENSE").write_text(
                "MIT License\nPermission is hereby granted, free of charge\n"
                'THE SOFTWARE IS PROVIDED "AS IS"\n',
                encoding="utf-8",
            )

            result = run_validator(root)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("SKILL.md exceeds 500 lines", result.stdout + result.stderr)

    def test_markdown_python_fence_must_parse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            skill = root / "example-skill"
            skill.mkdir()
            (skill / "SKILL.md").write_text(
                "---\n"
                "name: example-skill\n"
                "description: Use when testing fenced Python examples.\n"
                "license: MIT\n"
                "---\n\n"
                "```python\nnot valid Python\n```\n",
                encoding="utf-8",
            )
            (root / "LICENSE").write_text(
                "MIT License\nPermission is hereby granted, free of charge\n"
                'THE SOFTWARE IS PROVIDED "AS IS"\n',
                encoding="utf-8",
            )

            result = run_validator(root)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("invalid Python code block", result.stdout + result.stderr)

    def test_repository_satisfies_contract(self) -> None:
        result = run_validator(ROOT)

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
