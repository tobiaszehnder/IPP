"""CLI argument parsing and the `main` entrypoint for the alignment orchestrator."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

from align.core import Orchestrator
from align.model import (
    STEP_NAMES,
    PipelineConfig,
    PipelineStepError,
)


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
