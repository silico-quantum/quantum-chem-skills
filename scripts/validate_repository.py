#!/usr/bin/env python3
"""Validate the repository's portable Agent Skills contract."""

from __future__ import annotations

import argparse
import ast
import re
import sys
from pathlib import Path
from urllib.parse import unquote


ALLOWED_FRONTMATTER_FIELDS = {
    "name",
    "description",
    "license",
    "compatibility",
    "metadata",
    "allowed-tools",
}
NAME_PATTERN = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
CJK_PATTERN = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")
MARKDOWN_LINK_PATTERN = re.compile(r"\[[^\]]*\]\(([^)\s]+)(?:\s+[^)]*)?\)")
PYTHON_FENCE_PATTERN = re.compile(r"```python[ \t]*\n(.*?)```", re.DOTALL)
MAX_SKILL_LINES = 500


def relative(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def parse_frontmatter(path: Path, root: Path) -> tuple[dict[str, str], list[str]]:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    errors: list[str] = []
    label = relative(path, root)

    if not lines or lines[0] != "---":
        return {}, [f"{label}: missing YAML frontmatter"]

    try:
        closing = lines.index("---", 1)
    except ValueError:
        return {}, [f"{label}: unclosed YAML frontmatter"]

    fields: dict[str, str] = {}
    for line_number, line in enumerate(lines[1:closing], start=2):
        if not line or line[0].isspace():
            continue
        if ":" not in line:
            errors.append(f"{label}:{line_number}: malformed frontmatter entry")
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        fields[key] = value
        if key not in ALLOWED_FRONTMATTER_FIELDS:
            errors.append(
                f"{label}:{line_number}: unsupported frontmatter field '{key}'"
            )

    name = fields.get("name", "")
    description = fields.get("description", "")
    if not name:
        errors.append(f"{label}: frontmatter requires 'name'")
    elif not NAME_PATTERN.fullmatch(name) or len(name) > 64:
        errors.append(f"{label}: invalid skill name '{name}'")
    elif name != path.parent.name:
        errors.append(
            f"{label}: name must match its parent directory ({path.parent.name})"
        )

    if not description:
        errors.append(f"{label}: frontmatter requires 'description'")
    elif len(description) > 1024:
        errors.append(f"{label}: description exceeds 1024 characters")
    elif not description.startswith("Use when"):
        errors.append(f"{label}: description must start with 'Use when'")

    compatibility = fields.get("compatibility")
    if compatibility is not None and not 1 <= len(compatibility) <= 500:
        errors.append(f"{label}: compatibility must contain 1-500 characters")

    return fields, errors


def check_markdown(path: Path, root: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    label = relative(path, root)
    errors: list[str] = []

    if CJK_PATTERN.search(text):
        errors.append(f"{label}: contains CJK text; public Markdown must be English")

    for match in PYTHON_FENCE_PATTERN.finditer(text):
        try:
            ast.parse(match.group(1), filename=str(path))
        except SyntaxError as exc:
            block_line = text.count("\n", 0, match.start(1)) + 1
            error_line = block_line + (exc.lineno or 1) - 1
            errors.append(
                f"{label}:{error_line}: invalid Python code block: {exc.msg}"
            )

    for match in MARKDOWN_LINK_PATTERN.finditer(text):
        raw_target = match.group(1)
        if raw_target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target_text = unquote(raw_target.split("#", 1)[0])
        if not target_text:
            continue
        target = (path.parent / target_text).resolve()
        if not target.exists():
            errors.append(f"{label}: broken local link '{raw_target}'")

    return errors


def check_reference_discovery(skill_file: Path, root: Path) -> list[str]:
    """Require supplemental Markdown to be directly discoverable from SKILL.md."""
    skill_text = skill_file.read_text(encoding="utf-8")
    linked_paths: set[Path] = set()
    for match in MARKDOWN_LINK_PATTERN.finditer(skill_text):
        raw_target = match.group(1)
        if raw_target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        target_text = unquote(raw_target.split("#", 1)[0])
        if target_text:
            linked_paths.add((skill_file.parent / target_text).resolve())

    errors: list[str] = []
    for reference in sorted(skill_file.parent.rglob("*.md")):
        if reference.name in {"SKILL.md", "README.md"}:
            continue
        if reference.resolve() not in linked_paths:
            errors.append(
                f"{relative(reference, root)}: not linked from "
                f"{relative(skill_file, root)}"
            )
    return errors


def check_skill_size(skill_file: Path, root: Path) -> list[str]:
    """Keep the always-loaded skill entrypoint within a small context budget."""
    line_count = len(skill_file.read_text(encoding="utf-8").splitlines())
    if line_count > MAX_SKILL_LINES:
        return [
            f"{relative(skill_file, root)}: SKILL.md exceeds "
            f"{MAX_SKILL_LINES} lines ({line_count}); move details to linked references"
        ]
    return []


def check_license(root: Path) -> list[str]:
    path = root / "LICENSE"
    if not path.is_file():
        return ["LICENSE: missing repository license"]
    text = path.read_text(encoding="utf-8")
    required = (
        "MIT License",
        "Permission is hereby granted, free of charge",
        'THE SOFTWARE IS PROVIDED "AS IS"',
    )
    if not all(fragment in text for fragment in required):
        return ["LICENSE is not a complete MIT license"]
    return []


def check_python_source(root: Path) -> list[str]:
    errors: list[str] = []
    for path in sorted(root.rglob("*.py")):
        if ".git" in path.parts:
            continue
        text = path.read_text(encoding="utf-8")
        relative_parts = path.relative_to(root).parts
        if "tests" not in relative_parts and CJK_PATTERN.search(text):
            errors.append(
                f"{relative(path, root)}: contains CJK text; public source must be English"
            )
        try:
            ast.parse(text, filename=str(path))
        except (SyntaxError, UnicodeDecodeError) as exc:
            location = f":{exc.lineno}" if isinstance(exc, SyntaxError) else ""
            errors.append(
                f"{relative(path, root)}{location}: Python syntax error: {exc.msg if isinstance(exc, SyntaxError) else exc}"
            )
    return errors


def validate_repository(root: Path) -> list[str]:
    root = root.resolve()
    errors: list[str] = []

    skill_files = sorted(root.rglob("SKILL.md"))
    if not skill_files:
        errors.append("repository contains no SKILL.md files")
    for skill_file in skill_files:
        _, frontmatter_errors = parse_frontmatter(skill_file, root)
        errors.extend(frontmatter_errors)
        errors.extend(check_skill_size(skill_file, root))
        errors.extend(check_reference_discovery(skill_file, root))

    for markdown_file in sorted(root.rglob("*.md")):
        if ".git" not in markdown_file.parts:
            errors.extend(check_markdown(markdown_file, root))

    errors.extend(check_license(root))
    errors.extend(check_python_source(root))
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate Agent Skills metadata, English Markdown, local links, license, and Python syntax."
    )
    parser.add_argument("root", nargs="?", default=".", type=Path)
    args = parser.parse_args(argv)

    errors = validate_repository(args.root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        print(f"Validation failed with {len(errors)} error(s).")
        return 1

    print("Repository validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
