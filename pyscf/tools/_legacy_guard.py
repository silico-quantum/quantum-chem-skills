"""Shared fail-closed guard for quarantined historical command-line tools."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import NoReturn


def refuse_direct_execution(script: str) -> NoReturn:
    name = Path(script).name
    print(
        f"QUARANTINED: {name} is a historical prototype and direct execution "
        "is disabled. Use ../scripts/run_safe_dft_tda.py for the supported "
        "closed-shell RKS/TDA workflow; inspect tools/README.md for scope.",
        file=sys.stderr,
    )
    raise SystemExit(64)
