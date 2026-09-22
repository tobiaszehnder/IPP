"""Core projection: score thresholds, loading pwalns, and running project_coords."""

import ipp_cpp  # type: ignore
import numpy as np
import pandas as pd
import tqdm

from ipp.debug_output import debug_shortest_path
from ipp.logging_utils import ProjectionLogger


GENOME_SIZE_BASIS = 2728222451  # mouse mm39 genome size, the scoring reference genome
HALF_LIFE_DISTANCE = 10000  # distance at which score is 0.5 (scaled to genome size)

ANCHOR_COLUMNS = [
    "ref_anchor_direct_left_start",
    "ref_anchor_direct_left_end",
    "ref_anchor_direct_right_start",
    "ref_anchor_direct_right_end",
    "ref_anchor_multi_left_start",
    "ref_anchor_multi_left_end",
    "ref_anchor_multi_right_start",
    "ref_anchor_multi_right_end",
    "qry_anchor_direct_left_start",
    "qry_anchor_direct_left_end",
    "qry_anchor_direct_right_start",
    "qry_anchor_direct_right_end",
    "qry_anchor_multi_left_start",
    "qry_anchor_multi_left_end",
    "qry_anchor_multi_right_start",
    "qry_anchor_multi_right_end",
]

RESULT_COLUMNS = [
    "sort_idx",
    "id",
    "coords_ref",
    "coords_direct",
    "coords_multi",
    "score_direct",
    "score_multi",
    "bridging_species",
    *ANCHOR_COLUMNS,
]

# Columns kept when include_anchors is False.
BASE_COLUMNS = [
    "coords_ref",
    "coords_direct",
    "coords_multi",
    "score_direct",
    "score_multi",
    "bridging_species",
]


def score_threshold_from_distance(
    distance: float,
    genome_size_ref: int,
    half_life_distance: float = HALF_LIFE_DISTANCE,
) -> float:
    """Convert a bp distance threshold to a score threshold, rounded to 3 decimals."""
    return round(
        0.5 ** (distance * GENOME_SIZE_BASIS / (genome_size_ref * half_life_distance)),
        3,
    )


def load_ipp(
    path_pwaln: str, half_life_distance: float = HALF_LIFE_DISTANCE
) -> ipp_cpp.Ipp:
    """Load a pwaln binary file into a new ipp_cpp.Ipp instance."""
    ipp = ipp_cpp.Ipp()
    ipp.load_pwalns(path_pwaln)
    ipp.set_half_life_distance(half_life_distance)
    return ipp


def read_regions(regions_file: str) -> tuple[list, dict]:
    """Parse a BED file into centered ipp_cpp.Coords and a coords -> (idx, name) map."""
    ref_coords = []
    coord_names = {}
    with open(regions_file) as f:
        for i, line in enumerate(f.readlines()):
            cols = line.strip().split("\t")
            name = cols[3]
            ref_chrom = cols[0]
            ref_loc = int(np.mean([int(cols[1]), int(cols[2])]))
            coords = ipp_cpp.Coords((ref_chrom, ref_loc))
            ref_coords.append(coords)
            coord_names[coords] = (i, name)
    return ref_coords, coord_names


def run_projection(
    ipp,
    ref: str,
    qry: str,
    ref_coords: list,
    coord_names: dict,
    n_cores: int,
    include_anchors: bool,
    simple_coords: bool,
    logger: ProjectionLogger,
) -> tuple[pd.DataFrame, list]:
    """Run project_coords, assembling the raw results DataFrame, indexed by name.

    Returns (results_df, unmapped_region_names). score_direct/score_multi are floored
    to 3 decimals (a score of 1.0 is reserved for regions overlapping alignments
    throughout the path).
    """
    pbar = tqdm.tqdm(total=len(ref_coords), leave=False)
    logger.pbar = pbar

    results = []
    unmapped_regions = []

    # Start the projection
    logger.log("Projecting regions from %s to %s" % (ref, qry))

    def on_job_done_callback(
        ref_coord,
        direct_score,
        direct_coords,
        direct_up_anchor,
        direct_down_anchor,
        multi_shortest_path,
    ):
        """Callback for one finished mapping job; appends a result or unmapped name."""
        pbar.update()
        coord_idx, coord_name = coord_names[ref_coord]

        logger.debug()
        if not len(multi_shortest_path):
            logger.debug(
                "no mapping found for {} ({}:{})".format(
                    coord_name, ref_coord.chrom, ref_coord.loc
                )
            )
            logger.debug()
            unmapped_regions.append(coord_name)
            return

        logger.debug("({})".format(coord_name))
        debug_shortest_path(logger, multi_shortest_path, simple_coords)

        direct_refs = (
            tuple(
                map(
                    str,
                    (
                        direct_up_anchor.ref_start,
                        direct_up_anchor.ref_end,
                        direct_down_anchor.ref_start,
                        direct_down_anchor.ref_end,
                    ),
                )
            )
            if direct_up_anchor
            else ("", "", "", "")
        )
        direct_qrys = (
            tuple(
                map(
                    str,
                    (
                        direct_up_anchor.qry_start,
                        direct_up_anchor.qry_end,
                        direct_down_anchor.qry_start,
                        direct_down_anchor.qry_end,
                    ),
                )
            )
            if direct_up_anchor
            else ("", "", "", "")
        )

        # The first intermediate species on the path.
        multi_first_entry = multi_shortest_path[1]

        # The qry species on the shortest path.
        multi_last_entry = multi_shortest_path[-1]
        assert multi_last_entry.species == qry
        multi_score = multi_last_entry.score
        multi_coords = multi_last_entry.coords

        multi_bridging_species = ",".join(
            [spe.species for spe in multi_shortest_path[1:-1]]
        )
        results.append(
            [
                coord_idx,
                coord_name,
                ref_coord,
                direct_coords,
                multi_coords,
                direct_score,
                multi_score,
                multi_bridging_species,
                *direct_refs,
                multi_first_entry.up_anchor.ref_start,
                multi_first_entry.up_anchor.ref_end,
                multi_first_entry.down_anchor.ref_start,
                multi_first_entry.down_anchor.ref_end,
                *direct_qrys,
                multi_last_entry.up_anchor.qry_start,
                multi_last_entry.up_anchor.qry_end,
                multi_last_entry.down_anchor.qry_start,
                multi_last_entry.down_anchor.qry_end,
            ]
        )

    ipp.project_coords(ref, qry, ref_coords, n_cores, on_job_done_callback)
    pbar.close()
    logger.pbar = None

    results_df = (
        pd.DataFrame(results, columns=RESULT_COLUMNS)
        .sort_values(by=["sort_idx"])
        .drop(columns=["sort_idx"])
        .set_index("id")
    )

    if not include_anchors:
        results_df = results_df.loc[:, BASE_COLUMNS]

    # Trim scores to third decimal (floor) such that 0.9999 becomes 0.999
    # instead of 1.000.
    results_df[["score_direct", "score_multi"]] = results_df[
        ["score_direct", "score_multi"]
    ].apply(lambda x: np.floor(x * 1000) / 1000)

    return results_df, unmapped_regions
