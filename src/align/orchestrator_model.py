"""Core data model for the alignment orchestrator: directory layout, pipeline
errors, and the command/node types used to describe a build step.
"""

from __future__ import annotations

import shlex
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path


# Valid --force step names (also used to identify a failing step in errors).
STEP_NAMES = frozenset(
    {
        "get_fasta",
        "fa_to_twobit",
        "twobit_to_sizes",
        "lastdb",
        "lastal",
        "maf_to_psl",
        "axt_chain",
        "all_chain",
        "all_pre_chain",
        "collect_and_binarize_pwalns",
    }
)


class PipelineStepError(RuntimeError):
    """Raised when a pipeline step fails; carries enough context to report clearly."""

    def __init__(self, step: str, label: str, message: str, stderr: str = ""):
        super().__init__(message)
        self.step = step
        self.label = label
        self.stderr = stderr


@dataclass
class PipelineConfig:
    data_dir: Path

    @property
    def assembly_dir(self) -> Path:
        return self.data_dir / "assembly"

    @property
    def alignment_dir(self) -> Path:
        return self.data_dir / "alignment"

    @property
    def fasta_dir(self) -> Path:
        return self.data_dir / "fasta"

    @property
    def lastdb_dir(self) -> Path:
        return self.alignment_dir / "lastdb"

    @property
    def chain_dir(self) -> Path:
        return self.alignment_dir / "chain"

    @property
    def maf_dir(self) -> Path:
        return self.alignment_dir / "maf"

    @property
    def psl_dir(self) -> Path:
        return self.alignment_dir / "psl"

    @property
    def bed_dir(self) -> Path:
        return self.alignment_dir / "bed"

    @property
    def pwaln_dir(self) -> Path:
        return self.alignment_dir / "pwaln"

    def make_dirs(self) -> None:
        for d in (
            self.assembly_dir,
            self.alignment_dir,
            self.fasta_dir,
            self.lastdb_dir,
            self.chain_dir,
            self.maf_dir,
            self.psl_dir,
            self.bed_dir,
            self.pwaln_dir,
        ):
            d.mkdir(parents=True, exist_ok=True)


@dataclass
class SubprocessCommand:
    argv: list[str]
    stdout: Path | None = None  # redirect stdout to this file
    log: Path | None = None  # redirect stdout+stderr to this file (e.g. axtChain)
    shell: str | None = None  # a literal shell pipeline instead of argv

    def describe(self) -> str:
        if self.shell is not None:
            return self.shell
        text = " ".join(shlex.quote(a) for a in self.argv)
        if self.stdout is not None:
            text += f" > {self.stdout}"
        if self.log is not None:
            text += f" &> {self.log}"
        return text


@dataclass
class PyCommand:
    description: str
    fn: Callable[[], None]


Command = SubprocessCommand | PyCommand


@dataclass
class Node:
    name: str  # step/rule name, used for --force and error messages
    output: Path
    inputs: tuple[Path, ...]
    ancient_inputs: frozenset[Path]
    label: str  # human-readable description, e.g. "mm39 vs galGal6"
    commands: list[Command]
    threads: int = 1
    needs_mem_slot: bool = False
