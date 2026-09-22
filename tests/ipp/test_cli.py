import pytest

from ipp import cli


def test_build_parser_defaults():
    parser = cli.build_parser()
    args = parser.parse_args(["regions.bed", "mm39", "hg38", "pwaln.bin"])
    assert args.score_DC == 0.98
    assert args.score_IC == 0.84
    assert args.distance_DC == 300
    assert args.distance_IC == 2500
    assert args.distance_FC == 500
    assert args.n_cores == 1
    assert args.target_bedfile is None
    assert args.quiet is False
    assert args.verbose is False
    assert args.simple_coords is False
    assert args.include_anchors is False


def test_validate_score_thresholds_passes_when_dc_at_least_ic():
    cli.validate_score_thresholds(0.98, 0.84)
    cli.validate_score_thresholds(0.9, 0.9)


def test_validate_score_thresholds_exits_when_dc_below_ic():
    with pytest.raises(SystemExit) as exc_info:
        cli.validate_score_thresholds(0.5, 0.8)
    assert "score_DC" in str(exc_info.value)


def test_validate_input_files_passes_for_existing_files(tmp_path):
    regions_file = tmp_path / "regions.bed"
    regions_file.write_text("chr1\t0\t10\tpeak1\n")
    target_file = tmp_path / "target.bed"
    target_file.write_text("chr1\t0\t10\n")
    # should not raise
    cli.validate_input_files(str(regions_file), str(target_file))


def test_validate_input_files_passes_when_target_bedfile_is_none(tmp_path):
    regions_file = tmp_path / "regions.bed"
    regions_file.write_text("chr1\t0\t10\tpeak1\n")
    cli.validate_input_files(str(regions_file), None)


def test_validate_input_files_exits_with_clear_message_for_missing_regions_file(
    tmp_path,
):
    missing = tmp_path / "does_not_exist.bed"
    with pytest.raises(SystemExit) as exc_info:
        cli.validate_input_files(str(missing), None)
    message = str(exc_info.value)
    assert "not found" in message
    assert str(missing) in message


def test_validate_input_files_exits_with_clear_message_for_missing_target_bedfile(
    tmp_path,
):
    regions_file = tmp_path / "regions.bed"
    regions_file.write_text("chr1\t0\t10\tpeak1\n")
    missing_target = tmp_path / "does_not_exist_target.bed"
    with pytest.raises(SystemExit) as exc_info:
        cli.validate_input_files(str(regions_file), str(missing_target))
    assert "not found" in str(exc_info.value)
