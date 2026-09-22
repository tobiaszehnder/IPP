import subprocess
import sys
from pathlib import Path

from align import collect_pwalns


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
COLLECT_PWALNS_SCRIPT = REPO_ROOT / "src" / "align" / "collect_pwalns.py"


def _reset_chrom_globals():
    collect_pwalns.chroms = []
    collect_pwalns.chroms_map = {}


def _write_chain(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return path


# --------------------------------------------------------------------------
# read_alignment_blocks_from_chain: chain-file parsing and coordinate math.
# --------------------------------------------------------------------------


def test_gapped_blocks_accumulate_ref_and_qry_start(tmp_path):
    _reset_chrom_globals()
    chain = _write_chain(
        tmp_path,
        "sp1.sp2.all.pre.chain",
        "chain 100 chr1 1000 + 0 40 chr2 2000 + 0 40 1\n10 5 5\n10 5 5\n10\n",
    )

    arr = collect_pwalns.read_alignment_blocks_from_chain(str(chain))

    # 3 blocks of width 10, each followed by a 5bp gap on both ref and qry.
    assert arr[:, 1].tolist() == [0, 15, 30]  # ref_start
    assert arr[:, 3].tolist() == [0, 15, 30]  # qry_start (forward strand)
    assert arr[:, 4].tolist() == [10, 10, 10]  # length, MSB unset (forward)


def test_reverse_strand_flips_qry_coordinate_and_sets_strand_bit(tmp_path):
    _reset_chrom_globals()
    chain = _write_chain(
        tmp_path,
        "sp1.sp2.all.pre.chain",
        "chain 100 chr1 1000 + 100 110 chr2 2000 - 50 60 1\n10\n",
    )

    arr = collect_pwalns.read_alignment_blocks_from_chain(str(chain))

    assert arr.shape == (1, 5)
    ref_chrom, ref_start, qry_chrom, qry_start, length_and_strand = arr[0].tolist()
    assert ref_start == 100
    # qry_chrom_size - qry_start - 1 == 2000 - 50 - 1
    assert qry_start == 1949
    # MSB (1 << 15) marks the reverse strand; low 15 bits are the length.
    assert length_and_strand == (10 | (1 << 15))
    assert length_and_strand & 0x7FFF == 10


def test_chromosome_names_map_to_stable_indices(tmp_path):
    _reset_chrom_globals()
    chain = _write_chain(
        tmp_path,
        "sp1.sp2.all.pre.chain",
        "chain 100 chr1 1000 + 0 10 chr2 2000 + 0 10 1\n10\n",
    )

    collect_pwalns.read_alignment_blocks_from_chain(str(chain))

    assert collect_pwalns.chroms == ["chr1", "chr2"]
    assert collect_pwalns.chroms_map == {"chr1": 0, "chr2": 1}


def test_duplicate_alignment_blocks_are_deduplicated(tmp_path):
    _reset_chrom_globals()
    # Two chains with an identical single block each; np.unique must collapse
    # them to one row rather than keeping both.
    chain = _write_chain(
        tmp_path,
        "sp1.sp2.all.pre.chain",
        "chain 100 chr1 1000 + 0 10 chr2 2000 + 0 10 1\n10\n"
        "chain 100 chr1 1000 + 0 10 chr2 2000 + 0 10 2\n10\n",
    )

    arr = collect_pwalns.read_alignment_blocks_from_chain(str(chain))

    assert arr.shape == (1, 5)


# --------------------------------------------------------------------------
# main(): full binary round-trip, decoded byte-by-byte against the format
# documented at the top of collect_pwalns.py.
# --------------------------------------------------------------------------


def _read_cstr(f):
    b = bytearray()
    while True:
        c = f.read(1)
        if c in (b"", b"\x00"):
            break
        b += c
    return b.decode()


def _read_int(f, length):
    return int.from_bytes(f.read(length), byteorder=sys.byteorder)


def _decode_pwalns(path):
    with open(path, "rb") as f:
        version = _read_int(f, 1)
        magic = _read_int(f, 2)
        n_sp1 = _read_int(f, 1)
        species = {}
        for _ in range(n_sp1):
            sp1_name = _read_cstr(f)
            genome_size = _read_int(f, 8)
            n_sp2 = _read_int(f, 1)
            pairs = {}
            for _ in range(n_sp2):
                sp2_name = _read_cstr(f)
                n_ref_chrom_entries = _read_int(f, 4)
                entries = []
                for _ in range(n_ref_chrom_entries):
                    ref_chrom = _read_int(f, 4)
                    n_blocks = _read_int(f, 4)
                    blocks = [
                        (
                            _read_int(f, 4),  # ref_start
                            _read_int(f, 4),  # qry_start
                            _read_int(f, 4),  # qry_chrom
                            _read_int(f, 2),  # length_and_strand
                        )
                        for _ in range(n_blocks)
                    ]
                    entries.append((ref_chrom, blocks))
                pairs[sp2_name] = entries
            species[sp1_name] = {"genome_size": genome_size, "pairs": pairs}
        n_chroms = _read_int(f, 4)
        chrom_names = [_read_cstr(f) for _ in range(n_chroms)]
    return {
        "version": version,
        "magic": magic,
        "species": species,
        "chrom_names": chrom_names,
    }


def test_main_writes_correctly_decodable_binary_pwaln_file(tmp_path):
    assembly_dir = tmp_path / "assembly"
    chain_dir = tmp_path / "chain"
    assembly_dir.mkdir()
    chain_dir.mkdir()

    (assembly_dir / "sp1.sizes").write_text("chr1\t1000\n")
    (assembly_dir / "sp2.sizes").write_text("chr2\t2000\n")

    _write_chain(
        chain_dir,
        "sp1.sp2.all.pre.chain",
        "chain 100 chr1 1000 + 0 10 chr2 2000 + 0 10 1\n10\n",
    )
    _write_chain(
        chain_dir,
        "sp2.sp1.all.pre.chain",
        "chain 100 chr2 2000 + 5 15 chr1 1000 - 2 12 2\n10\n",
    )

    outfile = tmp_path / "out.pwaln.bin"
    result = subprocess.run(
        [
            sys.executable,
            str(COLLECT_PWALNS_SCRIPT),
            str(chain_dir),
            str(assembly_dir),
            "sp1,sp2",
            str(outfile),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    decoded = _decode_pwalns(outfile)
    assert decoded["version"] == collect_pwalns.FORMAT_VERSION
    assert decoded["magic"] == 0xAFFE
    assert decoded["chrom_names"] == ["chr1", "chr2"]

    sp1 = decoded["species"]["sp1"]
    assert sp1["genome_size"] == 1000
    ((ref_chrom, blocks),) = sp1["pairs"]["sp2"]
    assert ref_chrom == 0  # chr1
    assert blocks == [(0, 0, 1, 10)]  # ref_start, qry_start, qry_chrom(chr2), length

    sp2 = decoded["species"]["sp2"]
    assert sp2["genome_size"] == 2000
    ((ref_chrom, blocks),) = sp2["pairs"]["sp1"]
    assert ref_chrom == 1  # chr2
    # reverse strand: qry_start = qry_chrom_size(1000) - qStart(2) - 1 = 997,
    # and the strand bit is set on length_and_strand.
    assert blocks == [(5, 997, 0, 10 | (1 << 15))]
