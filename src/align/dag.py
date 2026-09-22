"""DAG construction: maps an output path to the Node that builds it.

One `_node_*` method per pipeline step (get_fasta, fa_to_twobit,
twobit_to_sizes, lastdb, lastal, maf_to_psl, axt_chain, all_chain,
all_pre_chain, pwaln). Mixed into `Orchestrator`, which supplies `cfg`,
`nthreads`, `species_list` and the `_nodes` cache.
"""

from __future__ import annotations

import shlex
import sys
from pathlib import Path

from align.commands import _remove_species_prefix
from align.model import (
    Node,
    PipelineConfig,
    PyCommand,
    SubprocessCommand,
)


COMPUTE_ALIGNMENTS_DIR = Path(__file__).resolve().parent
COLLECT_PWALNS_SCRIPT = COMPUTE_ALIGNMENTS_DIR / "collect_pwalns.py"


class DagBuilderMixin:
    """Node-construction methods, mixed into `Orchestrator`.

    Expects `self.cfg: PipelineConfig`, `self.nthreads: int`,
    `self.species_list: list[str] | None` and `self._nodes: dict[Path, Node]`.
    """

    cfg: PipelineConfig
    nthreads: int
    species_list: list[str] | None
    _nodes: dict[Path, Node]

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
            argv=[],
            shell=(
                f"wget -nc {shlex.quote(link)} -O /dev/stdout"
                f" | gzip -d > {shlex.quote(str(path))}"
            ),
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
            argv=[],
            shell=(
                f"twoBitInfo {shlex.quote(str(twobit))} stdout"
                f" | sort -k2rn > {shlex.quote(str(path))}"
            ),
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
