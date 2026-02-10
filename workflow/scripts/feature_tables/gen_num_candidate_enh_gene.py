import os
import click
import pandas as pd


def _populate_enhancer_count_from_tss(df, enhancers, is_upstream):
    enh_indexes = enhancers.index
    if is_upstream:
        # start counting from the enhancer closest to TSS
        enh_indexes = reversed(enh_indexes)

    count_from_tss = 0
    for enh_idx in enh_indexes:
        count_from_tss += 1
        df.loc[enh_idx, "NumCandidateEnhGene"] = count_from_tss

def determine_num_candidate_enh_gene(pred_df, out_file):
    # This part is fine
    df = pred_df.copy() # Use a copy to avoid SettingWithCopyWarning
    df["midpoint"] = ((df["start"] + df["end"]) / 2).astype("int")

    # Sort enhancers relative to their gene's TSS for correct counting
    # Upstream: sort descending (enhancers closest to TSS come first)
    # Downstream: sort ascending (enhancers closest to TSS come first)
    df = df.sort_values(by=["TargetGene", "midpoint"], ascending=[True, True])
    
    # Create masks for upstream and downstream enhancers
    is_downstream = df["midpoint"] > df["TargetGeneTSS"]
    is_upstream = df["midpoint"] < df["TargetGeneTSS"]

    # Calculate cumulative counts in a vectorized way
    # For upstream, we sort descending by midpoint so cumcount counts away from TSS
    df.loc[is_upstream, "NumCandidateEnhGene"] = df[is_upstream].sort_values("midpoint", ascending=False).groupby("TargetGene").cumcount() + 1
    
    # For downstream, the default ascending sort is correct
    df.loc[is_downstream, "NumCandidateEnhGene"] = df[is_downstream].groupby("TargetGene").cumcount() + 1

    df = df.fillna(value=0) # Regions on top of the TSS will be NaN, fill them
    df["NumCandidateEnhGene"] = df["NumCandidateEnhGene"].astype("int")
    df[["name", "TargetGene", "NumCandidateEnhGene"]].to_csv(
        out_file,
        sep="\t",
        index=False,
    )
    print("Saved num candidate enhancers")

# def determine_num_candidate_enh_gene(pred_df, out_file):
#     # Need df to be sorted by midpoint for each chromosome
#     pred_df["midpoint"] = ((pred_df["start"] + pred_df["end"]) / 2).astype("int")
#     df = pred_df.sort_values(by=["chr", "midpoint"], ascending=True).reset_index(
#         drop=True
#     )

#     gene_groups = df.groupby(["TargetGene", "TargetGeneTSS"])
#     for (gene, tss), indexes in gene_groups.groups.items():
#         enhancers = df.loc[indexes]
#         upstream_enh = enhancers[enhancers["midpoint"] < tss]
#         downstream_enh = enhancers[enhancers["midpoint"] > tss]
#         _populate_enhancer_count_from_tss(df, upstream_enh, is_upstream=True)
#         _populate_enhancer_count_from_tss(df, downstream_enh, is_upstream=False)

#     df = df.fillna(value=0)
#     df["NumCandidateEnhGene"] = df["NumCandidateEnhGene"].astype("int")
#     df[["name", "TargetGene", "NumCandidateEnhGene"]].to_csv(
#         out_file,
#         sep="\t",
#         index=False,
#     )
#     print("Saved num candidate enhancers")


@click.command()
@click.option("--abc_predictions")
@click.option("--out_file")
def main(abc_predictions, out_file):
    pred_df = pd.read_csv(abc_predictions, sep="\t", compression="gzip")
    if len(pred_df) == 0:
        raise Exception("Did not find any enhancers in the Predictions file")

    determine_num_candidate_enh_gene(pred_df, out_file)


if __name__ == "__main__":
    main()
