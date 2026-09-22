import ipp_cpp  # type: ignore
import pandas as pd
import pytest

from ipp.io_utils import build_output_paths, convert_df_to_bed, write_unmapped


BED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
]


def _make_results_df():
    return pd.DataFrame(
        {
            "coords_ref": [ipp_cpp.Coords(("chr1", 100))],
            "coords_multi": [ipp_cpp.Coords(("chr2", 200))],
            "score_multi": [0.75],
            "sequence_conservation": ["IC"],
            "functional_conservation": ["+"],
        },
        index=["peak1"],
    )


def test_convert_df_to_bed_ref_has_bed9_columns_in_order():
    df = _make_results_df()
    bed = convert_df_to_bed(df, "ref")
    assert bed.columns.tolist() == BED_COLUMNS


def test_convert_df_to_bed_ref_uses_ref_coords_for_chrom_start():
    df = _make_results_df()
    bed = convert_df_to_bed(df, "ref")
    assert bed.loc["peak1", "chrom"] == "chr1"
    assert bed.loc["peak1", "start"] == 100
    assert bed.loc["peak1", "end"] == 101
    # name embeds the *opposite* (multi) coords
    assert bed.loc["peak1", "name"] == "peak1_chr2:200_IC+"


def test_convert_df_to_bed_qry_swaps_coordinate_space():
    df = _make_results_df()
    bed = convert_df_to_bed(df, "qry")
    assert bed.loc["peak1", "chrom"] == "chr2"
    assert bed.loc["peak1", "start"] == 200
    assert bed.loc["peak1", "end"] == 201
    # name embeds the *opposite* (ref) coords when in qry mode
    assert bed.loc["peak1", "name"] == "peak1_chr1:100_IC+"


@pytest.mark.parametrize(
    "sequence_conservation,functional_conservation,expected_rgb",
    [
        ("DC", "+", "255,0,0"),
        ("DC", "-", "30,30,30"),
        ("IC", "+", "255,0,0"),
        ("IC", "-", "30,30,30"),
        ("NC", "+", "141,153,174"),
        ("NC", "-", "141,153,174"),
    ],
)
def test_convert_df_to_bed_item_rgb_per_class(
    sequence_conservation, functional_conservation, expected_rgb
):
    df = _make_results_df()
    df["sequence_conservation"] = sequence_conservation
    df["functional_conservation"] = functional_conservation
    bed = convert_df_to_bed(df, "ref")
    assert bed.loc["peak1", "itemRgb"] == expected_rgb


def test_build_output_paths():
    paths = build_output_paths("/out", "/data/my_regions.bed", "mm39", "hg38")
    assert paths == {
        "unmapped": "/out/my_regions.mm39-hg38.unmapped",
        "table": "/out/my_regions.mm39-hg38.proj",
        "bed_ref": "/out/my_regions.mm39.bed",
        "bed_qry": "/out/my_regions.hg38.bed",
    }


def test_write_unmapped_empty_list_writes_empty_file(tmp_path):
    out_path = tmp_path / "empty.unmapped"
    write_unmapped([], str(out_path))
    assert out_path.read_bytes() == b""


def test_write_unmapped_writes_one_region_per_line(tmp_path):
    out_path = tmp_path / "regions.unmapped"
    write_unmapped(["peak1", "peak2"], str(out_path))
    assert out_path.read_text() == "peak1\npeak2\n"
