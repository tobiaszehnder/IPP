import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from align import orchestrator as orch


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
ORCHESTRATOR_SCRIPT = REPO_ROOT / "src" / "align" / "orchestrator.py"
FAKE_BIN = Path(__file__).resolve().parent / "fake_bin"
TEST_FASTA_DIR = REPO_ROOT / "data" / "test_fasta"
SPECIES_FILE = REPO_ROOT / "data" / "species_test.txt"
SPECIES = ["sp1", "sp2", "sp3"]


def _fake_path_env(extra_bin: Path | None = None) -> dict:
    env = os.environ.copy()
    bindirs = [str(extra_bin)] if extra_bin else []
    bindirs.append(str(FAKE_BIN))
    env["PATH"] = f"{os.pathsep.join(bindirs)}{os.pathsep}{env.get('PATH', '')}"
    return env


def _make_data_dir(tmp_path: Path) -> Path:
    data_dir = tmp_path / "ipp-data"
    fasta_dir = data_dir / "fasta"
    fasta_dir.mkdir(parents=True)
    for sp in SPECIES:
        shutil.copy(TEST_FASTA_DIR / f"{sp}.fa", fasta_dir / f"{sp}.fa")
    return data_dir


def run_cli(args: list[str], env: dict | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(ORCHESTRATOR_SCRIPT), *args]
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


# --------------------------------------------------------------------------
# Pure logic tests: no subprocess, no fake tools needed.
# --------------------------------------------------------------------------


def test_parse_list_arg_comma_separated():
    assert orch._parse_list_arg("a,b, c") == ["a", "b", "c"]


def test_parse_list_arg_from_file(tmp_path):
    f = tmp_path / "species.txt"
    f.write_text("# a comment\nmm39 some_other_column\n\nhg38\n")
    assert orch._parse_list_arg(str(f)) == ["mm39", "hg38"]


def test_parse_force_none_and_all_pass_through():
    assert orch.parse_force("none") == "none"
    assert orch.parse_force("all") == "all"


def test_parse_force_space_separated_step_names():
    assert orch.parse_force("lastal lastdb") == frozenset({"lastal", "lastdb"})


def test_parse_force_rejects_unknown_step_name():
    with pytest.raises(ValueError):
        orch.parse_force("not_a_real_step")


def _resolve(
    species=None, targets=None, create_pwaln_collection=False, ref=None, qry=None
):
    cfg = orch.PipelineConfig(data_dir=Path("/unused"))
    args = SimpleNamespace(
        species=species,
        targets=targets,
        create_pwaln_collection=create_pwaln_collection,
        ref=ref,
        qry=qry,
    )
    return orch._resolve_targets(args, cfg)


def test_species_mode_builds_ordered_cross_product_excluding_self_pairs():
    targets, species_list = _resolve(species="sp1,sp2,sp3")
    names = {t.name for t in targets}
    assert species_list == ["sp1", "sp2", "sp3"]
    assert len(targets) == 6  # 3*2 ordered, non-self pairs
    assert "sp1.sp1.all.pre.chain" not in names
    # both directions of a pair are distinct, independent targets
    assert "sp1.sp2.all.pre.chain" in names
    assert "sp2.sp1.all.pre.chain" in names


def test_pwaln_targets_generated_for_every_other_species_when_qry_omitted():
    targets, _ = _resolve(
        species="sp1,sp2,sp3", create_pwaln_collection=True, ref="sp1", qry=None
    )
    pwaln_names = {t.name for t in targets if t.name.endswith(".pwaln.bin")}
    assert pwaln_names == {"sp1.sp2.pwaln.bin", "sp1.sp3.pwaln.bin"}


def test_pwaln_target_single_when_qry_given():
    targets, _ = _resolve(
        species="sp1,sp2,sp3", create_pwaln_collection=True, ref="sp1", qry="sp2"
    )
    pwaln_names = [t.name for t in targets if t.name.endswith(".pwaln.bin")]
    assert pwaln_names == ["sp1.sp2.pwaln.bin"]


def test_targets_mode_bypasses_species_cross_product():
    targets, species_list = _resolve(targets="alignment/chain/sp1.sp2.all.pre.chain")
    assert species_list is None
    assert targets == [Path("alignment/chain/sp1.sp2.all.pre.chain")]


def test_split_pair_rejects_self_pairs(tmp_path):
    orc = orch.Orchestrator(orch.PipelineConfig(data_dir=tmp_path))
    with pytest.raises(ValueError):
        orc._split_pair(Path("/x/mm39.mm39.maf"), ".maf")


def test_pwaln_node_without_species_list_raises_clear_error(tmp_path):
    cfg = orch.PipelineConfig(data_dir=tmp_path)
    cfg.make_dirs()
    orc = orch.Orchestrator(cfg, species_list=None)
    with pytest.raises(ValueError, match="species list"):
        orc._node_for(cfg.pwaln_dir / "sp1.sp2.pwaln.bin")


def test_lastal_thread_hint_capped_at_two(tmp_path):
    cfg = orch.PipelineConfig(data_dir=tmp_path)
    cfg.make_dirs()
    node_one = orch.Orchestrator(cfg, nthreads=1)._node_lastal(
        cfg.maf_dir / "sp1.sp2.maf"
    )
    node_many = orch.Orchestrator(cfg, nthreads=8)._node_lastal(
        cfg.maf_dir / "sp3.sp2.maf"
    )
    assert node_one.commands[0].argv[3] == "1"
    assert node_many.commands[0].argv[3] == "2"  # min(nthreads, 2)


def test_node_get_fasta_shell_command_quotes_path_with_space(tmp_path):
    data_dir = tmp_path / "dir with space"
    cfg = orch.PipelineConfig(data_dir=data_dir)
    cfg.make_dirs()
    node = orch.Orchestrator(cfg)._node_get_fasta(cfg.fasta_dir / "sp1.fa")

    shell = node.commands[0].shell
    assert shlex.split(shell)[-1] == str(cfg.fasta_dir / "sp1.fa")


def test_node_twobit_to_sizes_shell_command_quotes_path_with_space(tmp_path):
    data_dir = tmp_path / "dir with space"
    cfg = orch.PipelineConfig(data_dir=data_dir)
    cfg.make_dirs()
    twobit = cfg.assembly_dir / "sp1.2bit"
    node = orch.Orchestrator(cfg)._node_twobit_to_sizes(cfg.assembly_dir / "sp1.sizes")

    shell = node.commands[0].shell
    parts = shlex.split(shell)
    assert str(twobit) in parts
    assert parts[-1] == str(cfg.assembly_dir / "sp1.sizes")


def test_mem_semaphore_matches_detected_memory_budget(tmp_path):
    cfg = orch.PipelineConfig(data_dir=tmp_path)
    orc = orch.Orchestrator(cfg, nthreads=8)
    expected = max(1, orch._detect_total_mem_mb() // orch.LASTAL_MEM_MB)
    assert orc._mem_semaphore._value == expected


def test_remove_species_prefix_collapses_doubled_prefix(tmp_path):
    assembly_dir = tmp_path / "assembly"
    assembly_dir.mkdir()
    (assembly_dir / "spA.sizes").write_text("spAchr1\t1000\n")
    (assembly_dir / "spB.sizes").write_text("chr1\t1000\n")
    maf_file = tmp_path / "spA.spB.maf"
    maf_file.write_text("s spAspAchr1 0 10 + 1000 ACGT\ns spBchr1 0 10 + 1000 ACGT\n")

    orch._remove_species_prefix(assembly_dir, maf_file)

    text = maf_file.read_text()
    # spA's real chrom names embed "spA" -> doubled prefix collapses to one
    assert "spAchr1" in text
    assert "spAspAchr1" not in text
    # spB's real chrom names don't embed "spB" -> prefix is stripped entirely
    assert "spBchr1" not in text
    assert "chr1" in text


# --------------------------------------------------------------------------
# Integration tests: real orchestrator CLI, fake external tools on PATH.
# --------------------------------------------------------------------------


def test_full_pipeline_creates_expected_layout_and_pwaln_collection(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    result = run_cli(
        [
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
        env=_fake_path_env(),
    )
    assert result.returncode == 0, result.stderr

    for rel in (
        "assembly",
        "alignment/lastdb",
        "alignment/chain",
        "alignment/maf",
        "alignment/psl",
        "alignment/bed",
        "alignment/pwaln",
        "fasta",
    ):
        assert (data_dir / rel).is_dir(), f"missing directory: {rel}"

    pwaln = data_dir / "alignment" / "pwaln" / "sp1.sp2.pwaln.bin"
    assert pwaln.is_file()
    assert pwaln.stat().st_size > 0

    for sp1 in SPECIES:
        for sp2 in SPECIES:
            if sp1 == sp2:
                continue
            assert (
                data_dir / "alignment" / "chain" / f"{sp1}.{sp2}.all.pre.chain"
            ).is_file()


def test_pwaln_collection_built_once_and_shared_across_qry_targets(tmp_path):
    # Omitting -q requests one .pwaln.bin per other species (sp1.sp2.pwaln.bin,
    # sp1.sp3.pwaln.bin), but they all contain the identical full-species-list
    # collection -- collect_pwalns.py must run exactly once, not once per target.
    data_dir = _make_data_dir(tmp_path)
    result = run_cli(
        [
            "-s",
            str(SPECIES_FILE),
            "-r",
            "sp1",
            "-c",
            "-d",
            str(data_dir),
            "-@",
            "2",
        ],
        env=_fake_path_env(),
    )
    assert result.returncode == 0, result.stderr

    run_lines = [
        line for line in result.stdout.splitlines() if line.startswith("[run]")
    ]
    collection_builds = [line for line in run_lines if "collection over" in line]
    collection_copies = [
        line for line in run_lines if "copy of shared collection" in line
    ]
    assert len(collection_builds) == 1, (
        f"expected collect_pwalns.py to run exactly once, got: {collection_builds}"
    )
    assert len(collection_copies) == 2  # sp1.sp2 and sp1.sp3

    pwaln_dir = data_dir / "alignment" / "pwaln"
    sp2_pwaln = pwaln_dir / "sp1.sp2.pwaln.bin"
    sp3_pwaln = pwaln_dir / "sp1.sp3.pwaln.bin"
    assert sp2_pwaln.is_file()
    assert sp3_pwaln.is_file()
    assert sp2_pwaln.read_bytes() == sp3_pwaln.read_bytes()


def test_rerun_without_changes_skips_every_step(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    env = _fake_path_env()
    args = [
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
    ]

    first = run_cli(args, env=env)
    assert first.returncode == 0, first.stderr

    second = run_cli(args, env=env)
    assert second.returncode == 0, second.stderr
    assert "[run]" not in second.stdout
    assert "[skip]" in second.stdout


def test_force_all_reruns_every_step(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    env = _fake_path_env()
    args = [
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
    ]

    assert run_cli(args, env=env).returncode == 0

    result = run_cli([*args, "-f", "all"], env=env)
    assert result.returncode == 0, result.stderr
    assert "[skip]" not in result.stdout
    # every step name shows up in a forced full rebuild
    for step in orch.STEP_NAMES:
        assert step in result.stdout


def test_force_single_step_only_reruns_that_step(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    env = _fake_path_env()
    args = [
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
    ]

    assert run_cli(args, env=env).returncode == 0

    result = run_cli([*args, "-f", "lastdb"], env=env)
    assert result.returncode == 0, result.stderr
    run_lines = [
        line for line in result.stdout.splitlines() if line.startswith("[run]")
    ]
    assert len(run_lines) == len(SPECIES)  # one lastdb call per species
    assert all("lastdb" in line for line in run_lines)


def test_dry_run_plans_without_executing(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    result = run_cli(
        [
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
            "-n",
        ],
        env=_fake_path_env(),
    )
    assert result.returncode == 0, result.stderr
    assert "[dry-run]" in result.stdout
    assert not (data_dir / "alignment" / "lastdb" / "sp1.prj").exists()
    assert not (data_dir / "alignment" / "pwaln" / "sp1.sp2.pwaln.bin").exists()


def test_error_from_failing_step_is_surfaced_clearly_and_stops_that_branch(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    failing_bin = tmp_path / "failing_bin"
    failing_bin.mkdir()
    (failing_bin / "lastal").write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        'if "sp3.fa" in sys.argv[-1]:\n'
        '    sys.stderr.write("lastal: boom\\n")\n'
        "    sys.exit(3)\n"
        "print('##maf version=1')\n"
        "print('a score=100')\n"
        "print('s chr1 0 10 + 1000 ACGTACGTAC')\n"
        "print('s chr1 0 10 + 1000 ACGTACGTAC')\n"
    )
    (failing_bin / "lastal").chmod(0o755)

    result = run_cli(
        [
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
            "3",
        ],
        env=_fake_path_env(extra_bin=failing_bin),
    )

    assert result.returncode == 1
    assert "lastal" in result.stderr
    assert "boom" in result.stderr
    # a target whose DAG never touched the failing pair still completes...
    assert (data_dir / "alignment" / "chain" / "sp2.sp1.all.pre.chain").is_file()
    # ...but the pwaln collection (which needs every pair) must not appear
    assert not (data_dir / "alignment" / "pwaln" / "sp1.sp2.pwaln.bin").exists()


def test_independent_branches_run_in_parallel(tmp_path):
    data_dir = _make_data_dir(tmp_path)
    slow_bin = tmp_path / "slow_bin"
    slow_bin.mkdir()
    (slow_bin / "axtChain").write_text(
        "#!/usr/bin/env python3\n"
        "import sys, time\n"
        "time.sleep(0.4)\n"
        "out_path = sys.argv[-1]\n"
        "with open(out_path, 'w') as f:\n"
        "    f.write('chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1\\n100\\n\\n')\n"
    )
    (slow_bin / "axtChain").chmod(0o755)
    env = _fake_path_env(extra_bin=slow_bin)
    args = ["-s", str(SPECIES_FILE), "-d", str(data_dir), "-@"]

    start = time.monotonic()
    serial = run_cli([*args, "1"], env=env)
    serial_elapsed = time.monotonic() - start
    assert serial.returncode == 0, serial.stderr

    shutil.rmtree(data_dir)
    data_dir = _make_data_dir(tmp_path)
    env = _fake_path_env(extra_bin=slow_bin)
    start = time.monotonic()
    parallel = run_cli([*args, "6"], env=env)
    parallel_elapsed = time.monotonic() - start
    assert parallel.returncode == 0, parallel.stderr

    assert parallel_elapsed < serial_elapsed * 0.6


def test_help_flag_runs_without_error():
    result = run_cli(["--help"])
    assert result.returncode == 0
    assert "--species" in result.stdout
