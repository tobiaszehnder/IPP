import os
import shutil
import subprocess
import tempfile


SNAKEFILE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "compute_alignments",
    "Snakefile",
)


def _dry_run_job_count(data_dir):
    result = subprocess.run(
        [
            "uv",
            "run",
            "snakemake",
            "--snakefile",
            SNAKEFILE,
            "--config",
            "species_list=sp1,sp2",
            "create_pwaln_collection=False",
            "ref=None",
            "qry=None",
            f"data_dir={data_dir}",
            "-n",
            "--cores",
            "1",
            os.path.join(data_dir, "alignment", "chain", "sp1.sp2.all.pre.chain"),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def test_chainprenet_ignores_stale_sizes_s2_mtime():
    # sizes_s1 and sizes_s2 are static per-species files reused across many
    # species pairs. If sizes_s2 isn't marked ancient(), touching it (e.g. by
    # recomputing an unrelated pair involving that species) spuriously
    # triggers a chainPreNet rerun for every pair. Both sizes inputs must be
    # ancient() for the same reason allchain and sizes_s1 already are.
    data_dir = tempfile.mkdtemp()
    try:
        assembly_dir = os.path.join(data_dir, "assembly")
        chain_dir = os.path.join(data_dir, "alignment", "chain")
        os.makedirs(assembly_dir)
        os.makedirs(chain_dir)

        allchain = os.path.join(chain_dir, "sp1.sp2.all.chain")
        sizes_s1 = os.path.join(assembly_dir, "sp1.sizes")
        sizes_s2 = os.path.join(assembly_dir, "sp2.sizes")
        pre_chain = os.path.join(chain_dir, "sp1.sp2.all.pre.chain")

        for path, mtime in [
            (allchain, 1000),
            (sizes_s1, 1000),
            (sizes_s2, 1000),
            (pre_chain, 2000),
        ]:
            open(path, "w").close()
            os.utime(path, (mtime, mtime))

        stdout = _dry_run_job_count(data_dir)
        assert "Nothing to be done" in stdout

        # Touch sizes_s2 so it's newer than the already-built output.
        os.utime(sizes_s2, (3000, 3000))

        stdout = _dry_run_job_count(data_dir)
        assert "Nothing to be done" in stdout, (
            "chainPreNet was rebuilt just because sizes_s2's mtime changed; "
            "sizes_s2 must be wrapped in ancient() like sizes_s1"
        )
    finally:
        shutil.rmtree(data_dir)
