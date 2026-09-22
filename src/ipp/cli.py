"""CLI argument parser and validation for the IPP projection script."""

import argparse
import os


def build_parser() -> argparse.ArgumentParser:
    """Build the argparse parser for src/ipp/project.py."""
    parser = argparse.ArgumentParser(
        description=(
            "Interspecies Point Projections (IPP).\nA method for projecting genomic "
            "point coordinates between genomes with large evolutionary distances."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )  # RawTextHelpFormatter)
    parser.add_argument(
        "regions_file",
        help=(
            "Bed file containing genomic coordinates. "
            "Regions with width > 1 will be centered."
        ),
    )
    parser.add_argument("ref", help="Reference species")
    parser.add_argument("qry", help="Query species")
    parser.add_argument("path_pwaln")
    parser.add_argument(
        "-o", "--out_dir", default=os.getcwd(), help="Directory for output files"
    )
    parser.add_argument(
        "-sdc",
        "--score_DC",
        type=float,
        default=0.98,
        help="Score threshold for direct conservation detection",
    )
    parser.add_argument(
        "-sic",
        "--score_IC",
        type=float,
        default=0.84,
        help="Score threshold for indirect conservation detection",
    )
    parser.add_argument(
        "-ddc",
        "--distance_DC",
        type=float,
        default=300,
        help=(
            "Distance threshold for direct conservation detection.\n"
            "Regions up to this distance to the closest anchor will be "
            "considered as directly conserved.\nUsed instead of score_threshold."
        ),
    )
    parser.add_argument(
        "-dic",
        "--distance_IC",
        type=float,
        default=2500,
        help=(
            "Distance threshold for indirect conservation detection.\n"
            "Regions up to this distance to the closest anchor will be "
            "considered as indirectly conserved.\nUsed instead of score_threshold."
        ),
    )
    parser.add_argument(
        "-dfc",
        "--distance_FC",
        type=float,
        default=500,
        help=(
            "Distance threshold for functional conservation detection.\n"
            "Regions up to this distance to the closest region in the "
            "target_bedfile will be considered as functionally conserved."
        ),
    )
    parser.add_argument("-n", "--n_cores", type=int, default=1, help="Number of CPUs")
    parser.add_argument(
        "-t",
        "--target_bedfile",
        default=None,
        help=(
            "Functional regions in target species to check for overlap with "
            "projections for classification"
        ),
    )
    # parser.add_argument(
    #     '-dhl', '--half_life_distance', type=int, default=10000,
    #     help='distance to closest anchor point at which projection score is 0.5'
    # )
    parser.add_argument(
        "-q", "--quiet", action="store_true", help="Do not produce any log output"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Produce additional debugging output",
    )
    parser.add_argument(
        "-c",
        "--simple_coords",
        action="store_true",
        help="Make coord numbers in debug output as small as possible",
    )
    parser.add_argument(
        "-a",
        "--include_anchors",
        action="store_true",
        help="Include anchors in results table",
    )
    return parser


def validate_input_files(regions_file: str, target_bedfile: str | None) -> None:
    """Raise ValueError with a clear message if a required input file is missing."""
    for label, path in (
        ("regions file", regions_file),
        ("target bedfile", target_bedfile),
    ):
        if path is None:
            continue
        if not os.path.isfile(path):
            raise ValueError(f"{label} not found: {path}")


def validate_score_thresholds(score_dc: float, score_ic: float) -> None:
    """DC is the stricter class, so score_DC must not be lower than score_IC."""
    if score_dc < score_ic:
        raise ValueError("score_DC must not be lower than score_IC")
