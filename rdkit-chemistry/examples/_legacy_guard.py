#!/usr/bin/env python3
"""Fail-closed guard for unvalidated historical RDKit examples."""

from __future__ import annotations

from pathlib import Path


def refuse_legacy_example(script: str | Path) -> None:
    """Stop before optional imports or scientific work can begin."""
    name = Path(script).name
    raise SystemExit(
        f"QUARANTINED LEGACY EXAMPLE: {name} is retained for static review only. "
        "It does not satisfy the current embedding, force-field, provenance, "
        "fresh-output, or scientific-acceptance contract. Follow SKILL.md and "
        "write a new fail-closed workflow in a fresh run directory."
    )
