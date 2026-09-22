"""Classification of projections by sequence (DC/IC/NC) and function (+/-)."""

import pandas as pd
import pyranges as pr


def classify_conservation(
    df_projections: pd.DataFrame,
    target_regions: pr.PyRanges | None,
    thresh_dc: float,
    thresh_ic: float,
    maxgap: float,
) -> tuple[pd.Series, pd.Series]:
    """Classify sequence conservation (DC/IC/NC) and functional conservation (+/-).

    Functional conservation is "" for every row when target_regions is None; otherwise
    it is "+" when the projected (multi) coordinate falls within `maxgap` of a target
    region, else "-". target_regions is expanded in place by maxgap on both ends.
    """
    sequence_conservation = df_projections.apply(
        lambda x: "DC"
        if x["score_direct"] >= thresh_dc
        else "IC"
        if x["score_multi"] >= thresh_ic
        else "NC",
        axis=1,
    )
    # create PyRanges object of projections
    df_multi = pd.DataFrame(
        {
            "Name": df_projections.index,
            "Chromosome": [x.chrom for x in df_projections["coords_multi"]],
            "Start": [x.loc for x in df_projections["coords_multi"]],
        }
    )
    df_multi["End"] = df_multi["Start"] + 1
    pr_multi = pr.PyRanges(df_multi)
    # determine functional conservation by checking for overlap with target regions
    functional_conservation = [""] * df_multi.shape[0]
    if target_regions is not None:
        # expand target regions by maxgap
        target_regions.Start -= maxgap
        target_regions.End += maxgap
        # save names of projected regions that overlap with a target region
        names_fc = pr_multi.overlap(target_regions).Name.values
        functional_conservation = pd.Series("-", index=df_projections.index)
        functional_conservation[names_fc] = "+"
    return sequence_conservation, functional_conservation
