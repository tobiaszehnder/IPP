"""The `Orchestrator`: combines DAG construction and the execution scheduler
into the single stateful object the CLI drives.
"""

from __future__ import annotations

import threading
from pathlib import Path

from compute_alignments.orchestrator_commands import LASTAL_MEM_MB, _detect_total_mem_mb
from compute_alignments.orchestrator_dag import DagBuilderMixin
from compute_alignments.orchestrator_model import Node, PipelineConfig
from compute_alignments.orchestrator_scheduler import SchedulerMixin


class Orchestrator(DagBuilderMixin, SchedulerMixin):
    def __init__(
        self,
        cfg: PipelineConfig,
        *,
        nthreads: int = 1,
        dry_run: bool = False,
        force: str | frozenset[str] = "none",
        species_list: list[str] | None = None,
    ):
        self.cfg = cfg
        self.nthreads = max(1, nthreads)
        self.dry_run = dry_run
        self.force = force
        self.species_list = species_list
        self._nodes: dict[Path, Node] = {}
        self._mem_semaphore = threading.Semaphore(
            max(1, _detect_total_mem_mb() // LASTAL_MEM_MB)
        )
        self._log_lock = threading.Lock()

    def log(self, msg: str) -> None:
        with self._log_lock:
            print(msg, flush=True)
