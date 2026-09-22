#!/usr/bin/env python

import argparse
import os
import sys
from collections.abc import Iterable

import pyranges as pr

from ipp import classify, cli, io_utils, projection
from ipp.logging_utils import (
    LOG_LEVEL_DEBUG,
    LOG_LEVEL_LOG,
    LOG_LEVEL_QUIET,
    ProjectionLogger,
)


def main(argv: Iterable[str] | None = None) -> int:
    parser = cli.build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    try:
        return _run(args)
    except ValueError as exc:
        # Usage-type errors: bad input files, inconsistent thresholds.
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        if args.verbose:
            raise
        return 1


def _run(args: argparse.Namespace) -> int:
    cli.validate_input_files(args.regions_file, args.target_bedfile)

    level = LOG_LEVEL_LOG
    if args.quiet:
        level = LOG_LEVEL_QUIET
    elif args.verbose:
        level = LOG_LEVEL_DEBUG
    logger = ProjectionLogger(level)

    # define variables and create output directory
    os.makedirs(args.out_dir, exist_ok=True)

    logger.log("Loading pairwise alignments")
    ipp = projection.load_ipp(args.path_pwaln)

    # compute score thresholds if distance thresholds were passed
    # score = 0.5^{minDist * genomeSizeBasis / (genomeSize * halfLifeDistance_)}
    genome_size_ref = ipp.get_genome_size(args.ref)
    score_DC = args.score_DC
    score_IC = args.score_IC
    if args.distance_DC is not None:
        score_DC = projection.score_threshold_from_distance(
            args.distance_DC, genome_size_ref
        )
        logger.log(
            f"DC distance threshold passed ({args.distance_DC}). "
            f"Corresponding DC score threshold set to {score_DC}"
        )
    if args.distance_IC is not None:
        score_IC = projection.score_threshold_from_distance(
            args.distance_IC, genome_size_ref
        )
        logger.log(
            f"IC distance threshold passed ({args.distance_IC}). "
            f"Corresponding IC score threshold set to {score_IC}"
        )
    cli.validate_score_thresholds(score_DC, score_IC)

    logger.log("Reading regions from %s" % (args.regions_file))
    ref_coords, coord_names = projection.read_regions(args.regions_file)

    results_df, unmapped_regions = projection.run_projection(
        ipp,
        args.ref,
        args.qry,
        ref_coords,
        coord_names,
        args.n_cores,
        args.include_anchors,
        args.simple_coords,
        logger,
    )

    if logger.is_debug():
        logger.debug(results_df.to_string())

    out_paths = io_utils.build_output_paths(
        args.out_dir, args.regions_file, args.ref, args.qry
    )

    # write list of coord names of unmapped regions to file
    logger.log("Writing unmapped regions to:\n\t%s" % out_paths["unmapped"])
    io_utils.write_unmapped(unmapped_regions, out_paths["unmapped"])

    # classify projections according to conservation of sequence (DC/IC/NC)
    # and function (+/-)
    logger.log("Classifying projections")
    target_regions = None
    if args.target_bedfile is not None:
        target_regions = pr.read_bed(args.target_bedfile)
    sequence_conservation, functional_conservation = classify.classify_conservation(
        results_df, target_regions, score_DC, score_IC, args.distance_FC
    )
    results_df.insert(5, "sequence_conservation", sequence_conservation)
    results_df.insert(6, "functional_conservation", functional_conservation)

    if results_df.shape[0] == 0:
        logger.log("No regions from the input bed files could be projected\nDone")
        return 0

    # write results table to file
    logger.log(
        "Writing table with successful projections to:\n\t%s" % out_paths["table"]
    )
    io_utils.write_proj_table(results_df, out_paths["table"])

    # write results to bed files (bed9 format for both reference and target species)
    logger.log(
        "Writing output bed files:\n\t%s\n\t%s"
        % (out_paths["bed_ref"], out_paths["bed_qry"])
    )
    io_utils.write_bed_files(results_df, out_paths["bed_ref"], out_paths["bed_qry"])

    logger.log(
        "NOTE: If you compare projections from different reference genomes to "
        "each other (e.g. mouse -> chicken and chicken -> mouse), always use "
        "the SAME SCORE thresholds.\nIPP adjusts the scores to the genome sizes. "
        "Projections with same distance thresholds from different genomes are "
        "not comparable.\n"
    )
    logger.log("Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# vim: tabstop=4 shiftwidth=4 expandtab
