import ipp_cpp  # type: ignore
import pandas as pd
import pyranges as pr

from ipp.classify import classify_conservation


def _make_projections_df(rows):
    """rows: list of (name, score_direct, score_multi, chrom, loc)."""
    data = {
        "score_direct": [r[1] for r in rows],
        "score_multi": [r[2] for r in rows],
        "coords_multi": [ipp_cpp.Coords((r[3], r[4])) for r in rows],
    }
    return pd.DataFrame(data, index=[r[0] for r in rows])


def test_dc_when_score_direct_at_or_above_threshold():
    df = _make_projections_df([("peak1", 0.9, 0.5, "chr1", 100)])
    seq_cons, _ = classify_conservation(
        df, None, thresh_dc=0.9, thresh_ic=0.7, maxgap=0
    )
    assert seq_cons["peak1"] == "DC"


def test_ic_when_below_dc_but_score_multi_at_or_above_ic_threshold():
    df = _make_projections_df([("peak1", 0.5, 0.7, "chr1", 100)])
    seq_cons, _ = classify_conservation(
        df, None, thresh_dc=0.9, thresh_ic=0.7, maxgap=0
    )
    assert seq_cons["peak1"] == "IC"


def test_ic_boundary_is_inclusive():
    df = _make_projections_df([("peak1", 0.5, 0.7, "chr1", 100)])
    seq_cons, _ = classify_conservation(
        df, None, thresh_dc=0.9, thresh_ic=0.7, maxgap=0
    )
    assert seq_cons["peak1"] == "IC"


def test_nc_when_below_both_thresholds():
    df = _make_projections_df([("peak1", 0.5, 0.6999, "chr1", 100)])
    seq_cons, _ = classify_conservation(
        df, None, thresh_dc=0.9, thresh_ic=0.7, maxgap=0
    )
    assert seq_cons["peak1"] == "NC"


def test_functional_conservation_is_empty_string_without_target_regions():
    df = _make_projections_df([("peak1", 0.9, 0.9, "chr1", 100)])
    _, func_cons = classify_conservation(
        df, None, thresh_dc=0.9, thresh_ic=0.7, maxgap=500
    )
    assert list(func_cons) == [""]


def _target_regions(chrom, start, end):
    return pr.PyRanges(
        pd.DataFrame({"Chromosome": [chrom], "Start": [start], "End": [end]})
    )


def test_functional_conservation_plus_when_within_target_region():
    df = _make_projections_df([("peak1", 0.9, 0.9, "chr1", 5005)])
    target = _target_regions("chr1", 5000, 5010)
    _, func_cons = classify_conservation(
        df, target, thresh_dc=0.9, thresh_ic=0.7, maxgap=0
    )
    assert func_cons["peak1"] == "+"


def test_functional_conservation_minus_when_outside_target_region_and_maxgap():
    # "anchor" always overlaps, so pr_multi.overlap() isn't entirely empty
    # (pyranges drops the Name column on a fully-empty overlap result).
    df = _make_projections_df(
        [("anchor", 0.9, 0.9, "chr1", 5005), ("peak1", 0.9, 0.9, "chr1", 100)]
    )
    target = _target_regions("chr1", 5000, 5010)
    _, func_cons = classify_conservation(
        df, target, thresh_dc=0.9, thresh_ic=0.7, maxgap=500
    )
    assert func_cons["peak1"] == "-"


def test_functional_conservation_plus_just_inside_maxgap_expansion():
    # target region chr1:5000-5010, maxgap 500 -> expanded to chr1:4500-5510.
    # A hit at 4500 (the new, expanded start) falls just inside.
    df = _make_projections_df([("peak1", 0.9, 0.9, "chr1", 4500)])
    target = _target_regions("chr1", 5000, 5010)
    _, func_cons = classify_conservation(
        df, target, thresh_dc=0.9, thresh_ic=0.7, maxgap=500
    )
    assert func_cons["peak1"] == "+"


def test_functional_conservation_minus_just_outside_maxgap_expansion():
    # One bp further out than the expanded boundary (4499 < 4500) misses.
    # "anchor" always overlaps, so pr_multi.overlap() isn't entirely empty
    # (pyranges drops the Name column on a fully-empty overlap result).
    df = _make_projections_df(
        [("anchor", 0.9, 0.9, "chr1", 4500), ("peak1", 0.9, 0.9, "chr1", 4499)]
    )
    target = _target_regions("chr1", 5000, 5010)
    _, func_cons = classify_conservation(
        df, target, thresh_dc=0.9, thresh_ic=0.7, maxgap=500
    )
    assert func_cons["peak1"] == "-"
