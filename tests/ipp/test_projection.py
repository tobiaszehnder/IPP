from ipp import projection
from ipp.logging_utils import LOG_LEVEL_QUIET, ProjectionLogger


DUMMY_PWALN = "data/dummy.pwaln.bin"
DUMMY_REGIONS = "data/dummy_regions.bed"


def test_score_threshold_from_distance_at_half_life_distance_is_half():
    # genome_size_ref == GENOME_SIZE_BASIS collapses the formula to
    # 0.5 ** (distance / half_life_distance).
    score = projection.score_threshold_from_distance(
        distance=10000,
        genome_size_ref=projection.GENOME_SIZE_BASIS,
        half_life_distance=10000,
    )
    assert score == 0.5


def test_score_threshold_from_distance_is_rounded_to_3_decimals():
    score = projection.score_threshold_from_distance(
        distance=1,
        genome_size_ref=projection.GENOME_SIZE_BASIS,
        half_life_distance=10000,
    )
    assert score == round(score, 3)


def test_read_regions_centers_coordinates_and_maps_names():
    ref_coords, coord_names = projection.read_regions(DUMMY_REGIONS)
    assert len(ref_coords) == 4

    first = ref_coords[0]
    assert first.chrom == "chr1"
    assert first.loc == 125  # mean(90, 160)
    assert coord_names[first] == (0, "region1")

    last = ref_coords[-1]
    assert last.chrom == "chr2"
    assert last.loc == 410  # mean(400, 420)
    assert coord_names[last] == (3, "region4")


def test_run_projection_against_real_dummy_pwaln_fixture():
    ipp = projection.load_ipp(DUMMY_PWALN)
    ref_coords, coord_names = projection.read_regions(DUMMY_REGIONS)
    logger = ProjectionLogger(LOG_LEVEL_QUIET)

    results_df, unmapped_regions = projection.run_projection(
        ipp,
        "mm39",
        "hg38",
        ref_coords,
        coord_names,
        n_cores=1,
        include_anchors=True,
        simple_coords=False,
        logger=logger,
    )

    # The dummy pwaln fixture has no real alignment data for these regions,
    # so every region ends up unmapped -- this is the "all unmapped" case.
    assert sorted(unmapped_regions) == ["region1", "region2", "region3", "region4"]
    assert results_df.shape[0] == 0
    assert results_df.columns.tolist() == [
        c for c in projection.RESULT_COLUMNS if c not in ("sort_idx", "id")
    ]


def test_run_projection_excludes_anchor_columns_when_not_requested():
    ipp = projection.load_ipp(DUMMY_PWALN)
    ref_coords, coord_names = projection.read_regions(DUMMY_REGIONS)
    logger = ProjectionLogger(LOG_LEVEL_QUIET)

    results_df, _ = projection.run_projection(
        ipp,
        "mm39",
        "hg38",
        ref_coords,
        coord_names,
        n_cores=1,
        include_anchors=False,
        simple_coords=False,
        logger=logger,
    )

    assert results_df.columns.tolist() == projection.BASE_COLUMNS


def test_run_projection_with_zero_regions_returns_empty_frame():
    ipp = projection.load_ipp(DUMMY_PWALN)
    logger = ProjectionLogger(LOG_LEVEL_QUIET)

    results_df, unmapped_regions = projection.run_projection(
        ipp,
        "mm39",
        "hg38",
        [],
        {},
        n_cores=1,
        include_anchors=True,
        simple_coords=False,
        logger=logger,
    )

    assert results_df.shape[0] == 0
    assert unmapped_regions == []
