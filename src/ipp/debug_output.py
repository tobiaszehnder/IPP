"""Debug printing of a multi-species shortest path."""

import tabulate

from ipp.logging_utils import ProjectionLogger


def debug_shortest_path(
    logger: ProjectionLogger, shortest_path: list, simple: bool
) -> None:
    """Print the shortest path as a table; no-op unless logger is in debug mode.

    "out anchors" are the ref coords of the anchors of the *next* entry; "in anchors"
    are the qry coords of the anchors of the *current* entry (up/down anchors collapsed
    to one when identical). If `simple`, coords are shifted towards zero for
    readability.
    """
    if not logger.is_debug() or not len(shortest_path):
        return

    headers = [
        "species",
        "score",
        "chrom",
        "in anchors (strand)",
        "loc",
        "out anchors (strand)",
    ]
    data = []
    for i in range(len(shortest_path)):
        spe = shortest_path[i]

        not_first = i > 0
        not_last = i < (len(shortest_path) - 1)
        nxt_spe = shortest_path[i + 1] if not_last else None
        in_anchors = (
            [
                spe.up_anchor.qry_start,
                spe.up_anchor.qry_end,
                spe.down_anchor.qry_start,
                spe.down_anchor.qry_end,
            ]
            if not_first
            else []
        )
        out_anchors = (
            [
                nxt_spe.up_anchor.ref_start,
                nxt_spe.up_anchor.ref_end,
                nxt_spe.down_anchor.ref_start,
                nxt_spe.down_anchor.ref_end,
            ]
            if not_last
            else []
        )

        min_coord = min(in_anchors + out_anchors) if simple else 0
        in_anchors = [i - min_coord for i in in_anchors]
        out_anchors = [i - min_coord for i in out_anchors]

        data.append(
            [
                spe.species,
                spe.score,
                spe.coords.chrom,
                _anchors_str(in_anchors),
                spe.coords.loc - min_coord,
                _anchors_str(out_anchors),
            ]
        )

    logger.debug(tabulate.tabulate(data, headers=headers, stralign="center"))
    logger.debug()


def _anchors_str(anchors: list) -> str:
    """Format an anchor pair as start-end (collapsed if identical), with strand."""
    if not len(anchors):
        # First in or last out
        return ""

    up_str = "{}-{}".format(anchors[0], anchors[1])
    down_str = "{}-{}".format(anchors[2], anchors[3])
    return "{} ({})".format(
        up_str if up_str == down_str else ", ".join([up_str, down_str]),
        "+" if anchors[0] < anchors[1] else "-",
    )
