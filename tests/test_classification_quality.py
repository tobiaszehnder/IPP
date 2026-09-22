"""Classification-quality end-to-end test using a real LAST toolchain.

Skipped unless lastal/lastdb/maf-convert are on PATH -- most environments
won't have them, so this only runs where it can. Builds real alignments for
the synthetic test genomes (data/test_fasta), converts the resulting PSL
alignment to chain format directly (a single already-aligned PSL hit per
chromosome pair needs no further UCSC axtChain/chainMergeSort/chainPreNet
joining or netting), and checks that the resulting DC/IC/NC classification
matches data/test_regions.bed's designed regions (see commit c9ae530, which
placed 3 DC, 1 IC, and 1 NC region by scanning real chain files).

tests/test_end_to_end.py covers the align<->ipp wiring with fake tools;
this test covers actual classification correctness, which fake placeholder
alignment data can't exercise.
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
COLLECT_PWALNS_SCRIPT = REPO_ROOT / "src" / "align" / "collect_pwalns.py"
PROJECT_SCRIPT = REPO_ROOT / "src" / "ipp" / "project.py"
TEST_FASTA_DIR = REPO_ROOT / "data" / "test_fasta"
REGIONS_FILE = REPO_ROOT / "data" / "test_regions.bed"
SPECIES = ["sp1", "sp2", "sp3"]

EXPECTED_CLASSIFICATION = {
    "region1_dc": "DC",
    "region2_dc": "DC",
    "region3_ic": "IC",
    "region4_dc": "DC",
    "region5_nc": "NC",
}

requires_last = pytest.mark.skipif(
    not all(shutil.which(tool) for tool in ("lastal", "lastdb", "maf-convert")),
    reason="requires a real LAST toolchain (lastal/lastdb/maf-convert) on PATH",
)


def _write_sizes_from_fasta(fasta_path: Path, out_path: Path) -> None:
    sizes = {}
    chrom = None
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                sizes[chrom] = 0
            elif line:
                sizes[chrom] += len(line)
    with open(out_path, "w") as f:
        for chrom, length in sizes.items():
            f.write(f"{chrom}\t{length}\n")


def _psl_to_chain(psl_path: Path, out_path: Path) -> None:
    """Faithful PSL -> UCSC chain format conversion for a single
    already-aligned hit. LAST produces one coherent gapped alignment per
    chromosome pair, so there's nothing left to chain/net -- this is a
    direct field mapping, not new alignment logic. PSL column layout:
    https://genome.ucsc.edu/FAQ/FAQformat.html#format2
    """
    with open(psl_path) as f, open(out_path, "w") as out:
        chain_id = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith(("psLayout", "match", "---")):
                continue
            fields = line.split("\t")
            if len(fields) < 21:
                continue
            matches = int(fields[0])
            strand = fields[8]
            q_name, q_size, q_start, q_end = (
                fields[9],
                int(fields[10]),
                int(fields[11]),
                int(fields[12]),
            )
            t_name, t_size, t_start, t_end = (
                fields[13],
                int(fields[14]),
                int(fields[15]),
                int(fields[16]),
            )
            block_sizes = [int(x) for x in fields[18].rstrip(",").split(",")]
            q_starts = [int(x) for x in fields[19].rstrip(",").split(",")]
            t_starts = [int(x) for x in fields[20].rstrip(",").split(",")]

            chain_id += 1
            out.write(
                f"chain {matches} {t_name} {t_size} + {t_start} {t_end} "
                f"{q_name} {q_size} {strand} {q_start} {q_end} {chain_id}\n"
            )
            n = len(block_sizes)
            for i in range(n):
                if i < n - 1:
                    dt = t_starts[i + 1] - (t_starts[i] + block_sizes[i])
                    dq = q_starts[i + 1] - (q_starts[i] + block_sizes[i])
                    out.write(f"{block_sizes[i]}\t{dt}\t{dq}\n")
                else:
                    out.write(f"{block_sizes[i]}\n")
            out.write("\n")


@requires_last
def test_real_alignment_classification_matches_designed_regions(tmp_path):
    fasta_dir = tmp_path / "fasta"
    assembly_dir = tmp_path / "assembly"
    lastdb_dir = tmp_path / "lastdb"
    chain_dir = tmp_path / "chain"
    for d in (fasta_dir, assembly_dir, lastdb_dir, chain_dir):
        d.mkdir()

    for sp in SPECIES:
        shutil.copy(TEST_FASTA_DIR / f"{sp}.fa", fasta_dir / f"{sp}.fa")
        _write_sizes_from_fasta(fasta_dir / f"{sp}.fa", assembly_dir / f"{sp}.sizes")
        subprocess.run(
            ["lastdb", "-c", str(lastdb_dir / sp), str(fasta_dir / f"{sp}.fa")],
            check=True,
            capture_output=True,
        )

    for s1 in SPECIES:
        for s2 in SPECIES:
            if s1 == s2:
                continue
            maf = tmp_path / f"{s1}.{s2}.maf"
            psl = tmp_path / f"{s1}.{s2}.psl"
            with open(maf, "w") as f:
                subprocess.run(
                    ["lastal", str(lastdb_dir / s1), str(fasta_dir / f"{s2}.fa")],
                    check=True,
                    stdout=f,
                )
            with open(psl, "w") as f:
                subprocess.run(["maf-convert", "psl", str(maf)], check=True, stdout=f)
            _psl_to_chain(psl, chain_dir / f"{s1}.{s2}.all.pre.chain")

    pwaln = tmp_path / "sp1.sp2.pwaln.bin"
    collect_result = subprocess.run(
        [
            sys.executable,
            str(COLLECT_PWALNS_SCRIPT),
            str(chain_dir),
            str(assembly_dir),
            ",".join(SPECIES),
            str(pwaln),
        ],
        capture_output=True,
        text=True,
    )
    assert collect_result.returncode == 0, collect_result.stderr

    out_dir = tmp_path / "projection"
    out_dir.mkdir()
    project_result = subprocess.run(
        [
            sys.executable,
            str(PROJECT_SCRIPT),
            "-o",
            str(out_dir),
            "-n",
            "1",
            "-ddc",
            "0.0001",
            "-dic",
            "0.0001",
            str(REGIONS_FILE),
            "sp1",
            "sp2",
            str(pwaln),
        ],
        capture_output=True,
        text=True,
    )
    assert project_result.returncode == 0, project_result.stderr

    proj_file = out_dir / "test_regions.sp1-sp2.proj"
    assert proj_file.is_file(), "expected at least one region to be projected"
    with open(proj_file, newline="") as f:
        rows = {row["id"]: row for row in csv.DictReader(f, delimiter="\t")}

    for name, expected_class in EXPECTED_CLASSIFICATION.items():
        assert name in rows, f"{name} was not projected at all"
        actual_class = rows[name]["sequence_conservation"]
        assert actual_class == expected_class, (
            f"{name}: expected {expected_class}, got {actual_class}"
        )
