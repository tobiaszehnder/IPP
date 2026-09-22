"""Writing projection outputs: the .proj table, .unmapped list, and BED9 files."""

import os

import pandas as pd


BED_COLORS = {
    "DC": "127,201,127",
    "IC": "253,180,98",
    "NC": "141,153,174",
    "DC+": "255,0,0",
    "IC+": "255,0,0",
    "DC-": "30,30,30",
    "IC-": "30,30,30",
    "NC+": "141,153,174",
    "NC-": "141,153,174",
}


def build_output_paths(out_dir: str, regions_file: str, ref: str, qry: str) -> dict:
    """Compute the unmapped list, proj table, and BED file output paths."""
    basename = os.path.splitext(os.path.basename(regions_file))[0]
    prefix = os.path.join(out_dir, basename)
    return {
        "unmapped": f"{prefix}.{ref}-{qry}.unmapped",
        "table": f"{prefix}.{ref}-{qry}.proj",
        "bed_ref": f"{prefix}.{ref}.bed",
        "bed_qry": f"{prefix}.{qry}.bed",
    }


def convert_df_to_bed(df: pd.DataFrame, which: str) -> pd.DataFrame:
    """Convert projections to BED9 columns; which is "ref" or "qry" (coord space).

    name column: "id_oppositeCoords_class", e.g. "peak_75_chr3:3880_IC+".
    """
    coords_col = "coords_ref"
    opposite_coords_col = "coords_multi"
    if which == "qry":
        coords_col, opposite_coords_col = opposite_coords_col, coords_col
    bed = pd.DataFrame(df[coords_col].tolist(), index=df.index).rename(
        columns={0: "chrom", 1: "start"}
    )
    bed["start"] = bed["start"].astype(int)
    bed["end"] = bed["start"] + 1
    bed["name"] = df.apply(
        lambda x: "_".join(
            [
                x.name,
                str(x[opposite_coords_col]),
                x["sequence_conservation"] + x["functional_conservation"],
            ]
        ),
        axis=1,
    )
    bed["score"] = df["score_multi"]
    bed["strand"] = "."
    bed["thickStart"] = bed["start"]
    bed["thickEnd"] = bed["end"]
    bed["itemRgb"] = (
        df["sequence_conservation"] + df["functional_conservation"]
    ).apply(lambda x: BED_COLORS[x])
    return bed


def write_unmapped(unmapped_regions: list, out_path: str) -> None:
    """Write one unmapped region name per line."""
    with open(out_path, "w") as f:
        f.write("\n".join(unmapped_regions) + "\n")


def write_proj_table(results_df: pd.DataFrame, out_path: str) -> None:
    """Write the classified results table, tab-separated with scores to 3 decimals."""
    results_df.to_csv(out_path, sep="\t", header=True, float_format="%.3f")


def write_bed_files(
    results_df: pd.DataFrame, out_path_ref: str, out_path_qry: str
) -> None:
    """Write BED9 files (no header) for both reference and query coordinate spaces."""
    bed_ref = convert_df_to_bed(results_df, "ref")
    bed_qry = convert_df_to_bed(results_df, "qry")
    bed_ref.to_csv(
        out_path_ref, sep="\t", index=False, header=None, float_format="%.3f"
    )
    bed_qry.to_csv(
        out_path_qry, sep="\t", index=False, header=None, float_format="%.3f"
    )
