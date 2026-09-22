"""Execution/scheduling loop: staleness checking, concurrent execution, and
error propagation. Mixed into `Orchestrator`.
"""

from __future__ import annotations

import threading
from collections.abc import Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from pathlib import Path

from align.commands import _run_subprocess
from align.model import Node, PipelineStepError, PyCommand


class SchedulerMixin:
    """Staleness checking and concurrent execution, mixed into `Orchestrator`.

    Expects `self.force`, `self.dry_run`, `self.nthreads`, `self.log`,
    `self._mem_semaphore: threading.Semaphore` and `self._node_for`
    (from `DagBuilderMixin`).
    """

    force: str | frozenset[str]
    dry_run: bool
    nthreads: int
    _mem_semaphore: threading.Semaphore
    _nodes: dict[Path, Node]

    def _needs_rebuild(self, node: Node) -> bool:
        if self.force == "all":
            return True
        if isinstance(self.force, frozenset) and node.name in self.force:
            return True
        if not node.output.exists():
            return True
        out_mtime = node.output.stat().st_mtime
        for inp in node.inputs:
            if inp in node.ancient_inputs:
                continue
            if inp.exists() and inp.stat().st_mtime > out_mtime:
                return True
        return False

    def _execute_node(self, node: Node) -> None:
        if not self._needs_rebuild(node):
            self.log(f"[skip]     {node.name:<28} {node.label} (up to date)")
            return
        if self.dry_run:
            self.log(f"[dry-run]  {node.name:<28} {node.label}")
            for cmd in node.commands:
                if isinstance(cmd, PyCommand):
                    self.log(f"           $ {cmd.description}")
                else:
                    self.log(f"           $ {cmd.describe()}")
            return

        self.log(f"[run]      {node.name:<28} {node.label}")
        node.output.parent.mkdir(parents=True, exist_ok=True)
        mem_slot = self._mem_semaphore if node.needs_mem_slot else None
        if mem_slot is not None:
            mem_slot.acquire()
        try:
            for cmd in node.commands:
                if isinstance(cmd, PyCommand):
                    try:
                        cmd.fn()
                    except Exception as exc:
                        raise PipelineStepError(
                            node.name,
                            node.label,
                            f"step '{node.name}' failed for {node.label}: {exc}",
                        ) from exc
                else:
                    _run_subprocess(cmd, step=node.name, label=node.label)
        finally:
            if mem_slot is not None:
                mem_slot.release()

    def run(self, targets: Sequence[Path]) -> None:
        order: list[Path] = []
        seen: set[Path] = set()

        def visit(path: Path) -> None:
            if path in seen:
                return
            seen.add(path)
            node = self._node_for(path)
            for dep in node.inputs:
                visit(dep)
            order.append(path)

        for t in targets:
            visit(t)

        futures: dict[Path, Future[None]] = {}
        with ThreadPoolExecutor(max_workers=self.nthreads) as executor:
            for path in order:
                node = self._nodes[path]

                def task(node=node, futures=futures) -> None:
                    for dep in node.inputs:
                        futures[dep].result()
                    self._execute_node(node)

                futures[path] = executor.submit(task)
            # Propagate the first failure (if any) once every target is done.
            first_error: BaseException | None = None
            for t in targets:
                try:
                    futures[t].result()
                except BaseException as exc:  # noqa: BLE001 - re-raised below
                    if first_error is None:
                        first_error = exc
            if first_error is not None:
                raise first_error
