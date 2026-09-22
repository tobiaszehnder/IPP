#!/usr/bin/env python3
"""Orchestrates the pairwise-alignment pipeline used to build .pwaln files.

For each species pair, builds a small dependency graph of steps (fasta ->
faToTwoBit -> lastdb/lastal -> maf-convert -> axtChain -> chainMergeSort ->
chainPreNet -> collect_pwalns), skipping steps whose output is already
up to date, and runs independent species pairs concurrently.

Thin CLI entrypoint; the actual logic lives in the sibling orchestrator_*
modules (model, commands, dag, scheduler, core, cli).
"""

from __future__ import annotations

import sys
from pathlib import Path


# Allow running this file directly (`python compute_alignments/orchestrator.py`),
# where the repo root isn't on sys.path yet, as well as importing it as
# `compute_alignments.orchestrator` (e.g. from tests).
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from compute_alignments.orchestrator_cli import (
    _parse_list_arg,
    _resolve_targets,
    build_arg_parser,
    main,
    parse_force,
)
from compute_alignments.orchestrator_commands import (
    LASTAL_MEM_MB,
    _detect_total_mem_mb,
    _remove_species_prefix,
)
from compute_alignments.orchestrator_core import Orchestrator
from compute_alignments.orchestrator_model import (
    STEP_NAMES,
    PipelineConfig,
    PipelineStepError,
)


# Re-exported for backwards-compatible imports (e.g. `from compute_alignments
# import orchestrator as orch; orch.PipelineConfig(...)`).
__all__ = [
    "LASTAL_MEM_MB",
    "STEP_NAMES",
    "Orchestrator",
    "PipelineConfig",
    "PipelineStepError",
    "_detect_total_mem_mb",
    "_parse_list_arg",
    "_remove_species_prefix",
    "_resolve_targets",
    "build_arg_parser",
    "main",
    "parse_force",
]


if __name__ == "__main__":
    raise SystemExit(main())
