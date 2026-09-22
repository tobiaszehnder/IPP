"""main()'s exit-code/error-reporting contract: usage errors (bad input) exit
2 with a clean "error: ..." message, unexpected runtime errors exit 1 with a
clean "ERROR: ..." message (full traceback only under --verbose), and success
exits 0. Mirrors the same contract as align.cli.main -- see that module.
"""

import pytest

from ipp import project


DUMMY_REGIONS = "data/dummy_regions.bed"
DUMMY_PWALN = "data/dummy.pwaln.bin"


def test_missing_regions_file_exits_2_with_clean_error(tmp_path, capsys):
    out_dir = tmp_path / "out"
    exit_code = project.main(
        [
            "-o",
            str(out_dir),
            str(tmp_path / "does_not_exist.bed"),
            "mm39",
            "hg38",
            DUMMY_PWALN,
        ]
    )
    assert exit_code == 2
    captured = capsys.readouterr()
    assert captured.err.startswith("error:")
    assert "not found" in captured.err


def test_invalid_pwaln_file_exits_1_with_clean_error_not_traceback(tmp_path, capsys):
    bad_pwaln = tmp_path / "bad.pwaln.bin"
    bad_pwaln.write_text("not a real pwaln file")
    out_dir = tmp_path / "out"

    exit_code = project.main(
        ["-o", str(out_dir), DUMMY_REGIONS, "mm39", "hg38", str(bad_pwaln)]
    )

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "ERROR:" in captured.err
    assert "Traceback" not in captured.err


def test_invalid_pwaln_file_reraises_full_traceback_when_verbose(tmp_path):
    bad_pwaln = tmp_path / "bad.pwaln.bin"
    bad_pwaln.write_text("not a real pwaln file")
    out_dir = tmp_path / "out"

    with pytest.raises(RuntimeError):
        project.main(
            [
                "-v",
                "-o",
                str(out_dir),
                DUMMY_REGIONS,
                "mm39",
                "hg38",
                str(bad_pwaln),
            ]
        )


def test_successful_run_exits_0(tmp_path):
    out_dir = tmp_path / "out"
    exit_code = project.main(
        ["-o", str(out_dir), "-n", "1", DUMMY_REGIONS, "mm39", "hg38", DUMMY_PWALN]
    )
    assert exit_code == 0
