"""Command execution primitives: running a pipeline step's subprocess/shell
command and detecting the machine's memory budget. Independent of DAG logic.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

from align.model import PipelineStepError, SubprocessCommand


# lastal is memory-hungry; cap how many may run concurrently regardless of
# --nthreads so we don't exceed the machine's memory budget.
LASTAL_MEM_MB = 40000


def _run_subprocess(cmd: SubprocessCommand, *, step: str, label: str) -> None:
    if cmd.shell is not None:
        proc = subprocess.run(
            ["bash", "-c", f"set -o pipefail; {cmd.shell}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stderr_text = proc.stderr.decode(errors="replace") if proc.stderr else ""
    elif cmd.log is not None:
        with open(cmd.log, "wb") as log_f:
            proc = subprocess.run(cmd.argv, stdout=log_f, stderr=subprocess.STDOUT)
        stderr_text = ""
        if proc.returncode != 0:
            stderr_text = cmd.log.read_text(errors="replace")[-4000:]
    elif cmd.stdout is not None:
        with open(cmd.stdout, "wb") as out_f:
            proc = subprocess.run(cmd.argv, stdout=out_f, stderr=subprocess.PIPE)
        stderr_text = proc.stderr.decode(errors="replace") if proc.stderr else ""
    else:
        proc = subprocess.run(cmd.argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stderr_text = proc.stderr.decode(errors="replace") if proc.stderr else ""

    if proc.returncode != 0:
        raise PipelineStepError(
            step,
            label,
            f"step '{step}' failed for {label} (exit code {proc.returncode}): "
            f"{cmd.describe()}",
            stderr_text,
        )


def _remove_species_prefix(assembly_dir: Path, maf_file: Path) -> None:
    """Undo lastal's occasional genome-build-name prefixing of chromosome names.

    Pure-Python re-implementation of the old `remove_species_prefix` shell
    script (sed + GNU parallel): same substitution logic, streamed line by
    line so large real .maf files don't need to fit in memory.
    """
    s1, s2 = maf_file.name.removesuffix(".maf").split(".")
    subs: list[tuple[str, str]] = []
    for sp in (s1, s2):
        sizes_file = assembly_dir / f"{sp}.sizes"
        contains_sp = sizes_file.is_file() and sp in sizes_file.read_text()
        subs.append((sp + sp, sp) if contains_sp else (sp, ""))

    tmp_file = maf_file.with_name(maf_file.name + ".tmp")
    with open(maf_file) as src, open(tmp_file, "w") as dst:
        for line in src:
            for old, new in subs:
                line = line.replace(old, new)
            dst.write(line)
    tmp_file.replace(maf_file)


def _detect_total_mem_mb() -> int:
    if "SLURM_MEM_PER_NODE" in os.environ:
        return int(os.environ["SLURM_MEM_PER_NODE"])
    try:
        out = subprocess.run(
            ["free", "-m"], capture_output=True, text=True, check=True
        ).stdout
        for line in out.splitlines():
            if line.startswith("Mem:"):
                return int(line.split()[1])
    except (FileNotFoundError, subprocess.CalledProcessError):
        pass
    try:
        out = subprocess.run(
            ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, check=True
        ).stdout
        return int(out.strip()) // 1024 // 1024
    except (FileNotFoundError, subprocess.CalledProcessError, ValueError):
        pass
    return 8192
