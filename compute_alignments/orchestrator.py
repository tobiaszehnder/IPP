#!/usr/bin/env python3
"""Orchestrates the pairwise-alignment pipeline used to build .pwaln files.

For each species pair, builds a small dependency graph of steps (fasta ->
faToTwoBit -> lastdb/lastal -> maf-convert -> axtChain -> chainMergeSort ->
chainPreNet -> collect_pwalns), skipping steps whose output is already
up to date, and runs independent species pairs concurrently.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
import threading
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path


# lastal is memory-hungry; cap how many may run concurrently regardless of
# --nthreads so we don't exceed the machine's memory budget.
LASTAL_MEM_MB = 40000

COMPUTE_ALIGNMENTS_DIR = Path(__file__).resolve().parent
COLLECT_PWALNS_SCRIPT = COMPUTE_ALIGNMENTS_DIR / "collect_pwalns.py"

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


def _parse_list_arg(value: str) -> list[str]:
    """Comma-separated list, or a file with one entry in its first column per line."""
    path = Path(value)
    if path.is_file():
        items = []
        for raw_line in path.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            items.append(line.split()[0])
        return items
    return [v.strip() for v in value.split(",") if v.strip()]


def parse_force(value: str) -> str | frozenset[str]:
    if value in ("none", "all"):
        return value
    names = frozenset(value.split())
    unknown = names - STEP_NAMES
    if unknown:
        raise ValueError(
            f"unknown step name(s) in --force: {', '.join(sorted(unknown))}. "
            f"Valid step names: {', '.join(sorted(STEP_NAMES))}"
        )
    return names


class Orchestrator:
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

    # ---- node construction (dispatch on directory + filename) -------------

    def _node_for(self, path: Path) -> Node:
        node = self._nodes.get(path)
        if node is not None:
            return node
        node = self._build_node(path)
        self._nodes[path] = node
        return node

    def _build_node(self, path: Path) -> Node:
        cfg = self.cfg
        name = path.name
        if path.parent == cfg.fasta_dir and name.endswith(".fa"):
            return self._node_get_fasta(path)
        if path.parent == cfg.assembly_dir and name.endswith(".2bit"):
            return self._node_fa_to_twobit(path)
        if path.parent == cfg.assembly_dir and name.endswith(".sizes"):
            return self._node_twobit_to_sizes(path)
        if path.parent == cfg.lastdb_dir and name.endswith(".prj"):
            return self._node_lastdb(path)
        if path.parent == cfg.maf_dir and name.endswith(".maf"):
            return self._node_lastal(path)
        if path.parent == cfg.psl_dir and name.endswith(".psl"):
            return self._node_maf_to_psl(path)
        if path.parent == cfg.chain_dir and name.endswith(".all.pre.chain"):
            return self._node_all_pre_chain(path)
        if path.parent == cfg.chain_dir and name.endswith(".all.chain"):
            return self._node_all_chain(path)
        if path.parent == cfg.chain_dir and name.endswith(".axtChain"):
            return self._node_axt_chain(path)
        if path.parent == cfg.pwaln_dir and name.endswith(".pwaln.bin"):
            return self._node_pwaln(path)
        raise ValueError(f"don't know how to build target: {path}")

    @staticmethod
    def _split_pair(path: Path, suffix: str) -> tuple[str, str]:
        s1, s2 = path.name.removesuffix(suffix).split(".")
        if s1 == s2:
            raise ValueError(f"self-pairs are not supported: {path}")
        return s1, s2

    def _node_get_fasta(self, path: Path) -> Node:
        species = path.name.removesuffix(".fa")
        link = f"http://hgdownload.cse.ucsc.edu/goldenpath/{species}/bigZips/{species}.fa.gz"
        cmd = SubprocessCommand(
            argv=[], shell=f"wget -nc {link} -O /dev/stdout | gzip -d > {path}"
        )
        return Node(
            name="get_fasta",
            output=path,
            inputs=(),
            ancient_inputs=frozenset(),
            label=species,
            commands=[cmd],
        )

    def _node_fa_to_twobit(self, path: Path) -> Node:
        species = path.name.removesuffix(".2bit")
        fa = self.cfg.fasta_dir / f"{species}.fa"
        cmd = SubprocessCommand(argv=["faToTwoBit", str(fa), str(path)])
        return Node(
            name="fa_to_twobit",
            output=path,
            inputs=(fa,),
            ancient_inputs=frozenset({fa}),
            label=species,
            commands=[cmd],
        )

    def _node_twobit_to_sizes(self, path: Path) -> Node:
        species = path.name.removesuffix(".sizes")
        twobit = self.cfg.assembly_dir / f"{species}.2bit"
        cmd = SubprocessCommand(
            argv=[], shell=f"twoBitInfo {twobit} stdout | sort -k2rn > {path}"
        )
        return Node(
            name="twobit_to_sizes",
            output=path,
            inputs=(twobit,),
            ancient_inputs=frozenset({twobit}),
            label=species,
            commands=[cmd],
        )

    def _node_lastdb(self, path: Path) -> Node:
        species = path.name.removesuffix(".prj")
        fa = self.cfg.fasta_dir / f"{species}.fa"
        prefix = self.cfg.lastdb_dir / species
        threads = min(self.nthreads, 2)
        cmd = SubprocessCommand(
            argv=["lastdb", "-v", "-P", str(threads), "-c", str(prefix), str(fa)]
        )
        return Node(
            name="lastdb",
            output=path,
            inputs=(fa,),
            ancient_inputs=frozenset({fa}),
            label=species,
            commands=[cmd],
            threads=threads,
        )

    def _node_lastal(self, path: Path) -> Node:
        s1, s2 = self._split_pair(path, ".maf")
        prj = self.cfg.lastdb_dir / f"{s1}.prj"
        fa = self.cfg.fasta_dir / f"{s2}.fa"
        prefix = self.cfg.lastdb_dir / s1
        threads = min(self.nthreads, 2)
        align_cmd = SubprocessCommand(
            argv=["lastal", "-v", "-P", str(threads), str(prefix), str(fa)],
            stdout=path,
        )
        rename_cmd = PyCommand(
            description=f"remove_species_prefix {self.cfg.assembly_dir} {path}",
            fn=lambda: _remove_species_prefix(self.cfg.assembly_dir, path),
        )
        return Node(
            name="lastal",
            output=path,
            inputs=(prj, fa),
            ancient_inputs=frozenset({prj, fa}),
            label=f"{s1} vs {s2}",
            commands=[align_cmd, rename_cmd],
            threads=threads,
            needs_mem_slot=True,
        )

    def _node_maf_to_psl(self, path: Path) -> Node:
        s1, s2 = self._split_pair(path, ".psl")
        maf = self.cfg.maf_dir / f"{s1}.{s2}.maf"
        cmd = SubprocessCommand(argv=["maf-convert", "psl", str(maf)], stdout=path)
        return Node(
            name="maf_to_psl",
            output=path,
            inputs=(maf,),
            ancient_inputs=frozenset({maf}),
            label=f"{s1} vs {s2}",
            commands=[cmd],
        )

    def _node_axt_chain(self, path: Path) -> Node:
        s1, s2 = self._split_pair(path, ".axtChain")
        psl = self.cfg.psl_dir / f"{s1}.{s2}.psl"
        twobit_s1 = self.cfg.assembly_dir / f"{s1}.2bit"
        twobit_s2 = self.cfg.assembly_dir / f"{s2}.2bit"
        log_path = path.with_name(path.name + ".log")
        cmd = SubprocessCommand(
            argv=[
                "axtChain",
                "-psl",
                "-linearGap=loose",
                str(psl),
                str(twobit_s1),
                str(twobit_s2),
                str(path),
            ],
            log=log_path,
        )
        return Node(
            name="axt_chain",
            output=path,
            inputs=(psl, twobit_s1, twobit_s2),
            ancient_inputs=frozenset({psl, twobit_s1, twobit_s2}),
            label=f"{s1} vs {s2}",
            commands=[cmd],
        )

    def _node_all_chain(self, path: Path) -> Node:
        s1, s2 = self._split_pair(path, ".all.chain")
        axt_chain = self.cfg.chain_dir / f"{s1}.{s2}.axtChain"
        cmd = SubprocessCommand(argv=["chainMergeSort", str(axt_chain)], stdout=path)
        return Node(
            name="all_chain",
            output=path,
            inputs=(axt_chain,),
            ancient_inputs=frozenset({axt_chain}),
            label=f"{s1} vs {s2}",
            commands=[cmd],
        )

    def _node_all_pre_chain(self, path: Path) -> Node:
        s1, s2 = self._split_pair(path, ".all.pre.chain")
        all_chain = self.cfg.chain_dir / f"{s1}.{s2}.all.chain"
        sizes_s1 = self.cfg.assembly_dir / f"{s1}.sizes"
        sizes_s2 = self.cfg.assembly_dir / f"{s2}.sizes"
        cmd = SubprocessCommand(
            argv=[
                "chainPreNet",
                str(all_chain),
                str(sizes_s1),
                str(sizes_s2),
                str(path),
            ]
        )
        return Node(
            name="all_pre_chain",
            output=path,
            inputs=(all_chain, sizes_s1, sizes_s2),
            # Both sizes files are static per-species outputs reused across
            # many pairs; neither should trigger a rebuild on its own.
            ancient_inputs=frozenset({all_chain, sizes_s1, sizes_s2}),
            label=f"{s1} vs {s2}",
            commands=[cmd],
        )

    def _node_pwaln(self, path: Path) -> Node:
        ref, qry = self._split_pair(path, ".pwaln.bin")
        if self.species_list is None:
            raise ValueError(
                f"cannot build {path}: a pwaln collection requires the full "
                "species list (pass --species), even when using --targets"
            )
        chain_inputs = tuple(
            self.cfg.chain_dir / f"{sp1}.{sp2}.all.pre.chain"
            for sp1 in self.species_list
            for sp2 in self.species_list
            if sp1 != sp2
        )
        species_csv = ",".join(self.species_list)
        cmd = SubprocessCommand(
            argv=[
                sys.executable,
                str(COLLECT_PWALNS_SCRIPT),
                str(self.cfg.chain_dir),
                str(self.cfg.assembly_dir),
                species_csv,
                str(path),
            ]
        )
        return Node(
            name="collect_and_binarize_pwalns",
            output=path,
            inputs=chain_inputs,
            ancient_inputs=frozenset(),
            label=f"{ref} vs {qry} (collection over {species_csv})",
            commands=[cmd],
        )

    # ---- staleness / force ---------------------------------------------

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

    # ---- execution -------------------------------------------------------

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


def _resolve_targets(
    args: argparse.Namespace, cfg: PipelineConfig
) -> tuple[list[Path], list[str] | None]:
    if args.species:
        species_list = _parse_list_arg(args.species)
        targets = [
            cfg.chain_dir / f"{sp1}.{sp2}.all.pre.chain"
            for sp1 in species_list
            for sp2 in species_list
            if sp1 != sp2
        ]
        if args.create_pwaln_collection:
            if args.qry is None:
                targets += [
                    cfg.pwaln_dir / f"{args.ref}.{s}.pwaln.bin"
                    for s in species_list
                    if s != args.ref
                ]
            else:
                targets += [cfg.pwaln_dir / f"{args.ref}.{args.qry}.pwaln.bin"]
        return targets, species_list
    targets = [Path(t) for t in _parse_list_arg(args.targets)]
    return targets, None


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute pairwise genome alignments for IPP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--species",
        default=None,
        help=(
            "Comma-separated list or file (first column) with species (genome "
            "builds) for which to compute pairwise alignments"
        ),
    )
    parser.add_argument(
        "-t",
        "--targets",
        default=None,
        help="Comma-separated list or file with target files to be produced",
    )
    parser.add_argument(
        "-c",
        "--create_pwaln_collection",
        action="store_true",
        help="Create a pairwise alignment collection file. Used as input for IPP",
    )
    parser.add_argument(
        "-r", "--ref", default=None, help="Reference species (required with -c)"
    )
    parser.add_argument(
        "-q",
        "--qry",
        default=None,
        help="Query species. If omitted, generate for every other species",
    )
    parser.add_argument(
        "-d",
        "--data_dir",
        default="./ipp-data",
        help="Output directory. A predefined folder structure is created inside it",
    )
    parser.add_argument(
        "-f",
        "--force",
        default="none",
        help="'none' | 'all' | space-separated step names to force re-run",
    )
    parser.add_argument(
        "-@", "--nthreads", type=int, default=1, help="Number of parallel threads"
    )
    parser.add_argument(
        "-n",
        "--dry_run",
        action="store_true",
        help="Only print planned commands instead of executing them",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_arg_parser().parse_args(list(argv) if argv is not None else None)

    if not args.species and not args.targets:
        print("error: one of --species or --targets is required", file=sys.stderr)
        return 2
    if args.create_pwaln_collection and not args.ref:
        print(
            "error: --ref is required with --create_pwaln_collection", file=sys.stderr
        )
        return 2

    try:
        force = parse_force(args.force)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    cfg = PipelineConfig(data_dir=Path(args.data_dir).resolve())
    cfg.make_dirs()
    targets, species_list = _resolve_targets(args, cfg)

    orch = Orchestrator(
        cfg,
        nthreads=args.nthreads,
        dry_run=args.dry_run,
        force=force,
        species_list=species_list,
    )
    try:
        orch.run(targets)
    except PipelineStepError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        if exc.stderr:
            print("--- stderr ---", file=sys.stderr)
            print(exc.stderr, file=sys.stderr)
        return 1
    except ValueError as exc:
        # DAG-construction errors: unknown target pattern, self-pairs, a
        # pwaln target requested without --species, etc.
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
