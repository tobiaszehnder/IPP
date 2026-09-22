"""End-to-end test wiring the align and ipp packages together: build
alignments with the orchestrator, then project regions through them with
project.py, exactly like a real user would.

Uses the fake external-tool stand-ins from tests/align/fake_bin, so it runs
without lastal/UCSC tools installed. The fake lastal always emits the same
tiny placeholder alignment block (chr1:0-10), so none of data/test_regions.bed
falls inside it and every region ends up unmapped -- this test verifies the
two packages are wired together correctly (the .pwaln.bin format align
writes is what ipp reads, and the CLI hand-off works end to end), not
alignment quality or classification. tests/test_alignments.sh covers that
with real tools.
"""

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
ORCHESTRATOR_SCRIPT = REPO_ROOT / "src" / "align" / "orchestrator.py"
PROJECT_SCRIPT = REPO_ROOT / "src" / "ipp" / "project.py"
FAKE_BIN = REPO_ROOT / "tests" / "align" / "fake_bin"
TEST_FASTA_DIR = REPO_ROOT / "data" / "test_fasta"
SPECIES_FILE = REPO_ROOT / "data" / "species_test.txt"
REGIONS_FILE = REPO_ROOT / "data" / "test_regions.bed"
SPECIES = ["sp1", "sp2", "sp3"]


def _fake_path_env() -> dict:
    env = os.environ.copy()
    env["PATH"] = f"{FAKE_BIN}{os.pathsep}{env.get('PATH', '')}"
    return env


def test_align_then_project_round_trip(tmp_path):
    data_dir = tmp_path / "ipp-data"
    fasta_dir = data_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    for sp in SPECIES:
        shutil.copy(TEST_FASTA_DIR / f"{sp}.fa", fasta_dir / f"{sp}.fa")

    env = _fake_path_env()

    align_result = subprocess.run(
        [
            sys.executable,
            str(ORCHESTRATOR_SCRIPT),
            "-s",
            str(SPECIES_FILE),
            "-r",
            "sp1",
            "-q",
            "sp2",
            "-c",
            "-d",
            str(data_dir),
            "-@",
            "2",
        ],
        capture_output=True,
        text=True,
        env=env,
    )
    assert align_result.returncode == 0, align_result.stderr

    pwaln = data_dir / "alignment" / "pwaln" / "sp1.sp2.pwaln.bin"
    assert pwaln.is_file()
    assert pwaln.stat().st_size > 0

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

    basename = REGIONS_FILE.stem
    proj_file = out_dir / f"{basename}.sp1-sp2.proj"
    unmapped_file = out_dir / f"{basename}.sp1-sp2.unmapped"
    assert unmapped_file.is_file()

    proj_rows = []
    if proj_file.is_file():
        with open(proj_file, newline="") as f:
            proj_rows = list(csv.DictReader(f, delimiter="\t"))
    with open(unmapped_file) as f:
        unmapped_names = [line.strip() for line in f if line.strip()]

    with open(REGIONS_FILE) as f:
        input_names = [line.split("\t")[3].strip() for line in f if line.strip()]

    # Every input region is accounted for exactly once, either projected or
    # explicitly reported unmapped -- no region silently disappears.
    assert len(proj_rows) + len(unmapped_names) == len(input_names)
    # The fake alignment only covers chr1:0-10, so nothing in
    # data/test_regions.bed falls inside it -- every region is unmapped.
    # This asserts the current, deterministic fake-tool behavior; if the
    # fake fixtures ever grow richer alignment data, this will need updating.
    assert set(unmapped_names) == set(input_names)
