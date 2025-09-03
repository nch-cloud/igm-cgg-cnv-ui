#!/usr/bin/env python

# import plotly.express as px
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
import argparse
import json
import math
import os
import time

# import scipy.stats as st
import traceback

import numpy as np
import pandas as pd
from plot_interactive_cnvloh_bokeh import generate_chrom_plot as cnvloh_bokeh_plot
from plot_interactive_cnvloh_bokeh import (
    generate_genome_plot as cnvloh_bokeh_plot_genome,
)

#### FUNCTIONS ####

### Accessory functions -- these are things that might get used more than once for various CNV-LOH processing purposes ####


def order_genome_info(
    genome_info_dict: dict, mains_only: bool = True, size_cutoff: int = None
):
    autosomes = {}
    allosomes = []
    other = []
    for i in genome_info_dict:
        contig_len = int(genome_info_dict[i]["length"])
        if not size_cutoff or contig_len >= size_cutoff:
            try:
                num_index = int(i.lower().replace("chr", ""))
                autosomes[num_index] = i
            except:
                if i.lower().replace("chr", "") in ["w", "x", "y", "z"]:
                    allosomes.append(i)
                else:
                    other.append(i)
    chroms_out = []
    for i in sorted(list(autosomes.keys())):
        chroms_out.append(autosomes[i])
    chroms_out = chroms_out + sorted(allosomes)
    if not mains_only:
        chroms_out = chroms_out + sorted(other)
    new_genome_info = {}
    for i in chroms_out:
        new_genome_info[i] = genome_info_dict[i]
    return new_genome_info


## Turns a set of values into a confidence interval
# def points_to_ci(values_list:list, alpha=0.95):
#     lower_ci,upper_ci = st.t.interval(alpha=alpha, df=len(values_list)-1, loc=np.mean(values_list), scale=st.sem(values_list))
#     return lower_ci, upper_ci


## Takes the gene 'list' df and converts it to x,y coordinates to plot as rectangles ##
def genes_to_rects(gene_list_df, to: float = 1.0):
    x = []
    y = []
    for index, row in gene_list_df.iterrows():
        x.extend([row["START"], row["START"], row["END"], row["END"], None])
        y.extend([0, to, to, 0, None])
    try:
        x.pop()
        y.pop()
    except IndexError:
        pass
    return x, y


## Takes a CNV dataframe and turns into line segments ##
def cnv_to_segs(cnv_df):
    x = []
    y = []
    for index, row in cnv_df.iterrows():
        x.extend([row["START"], row["END"], None])
        y.extend([row["LOG2"], row["LOG2"], None])
    try:
        x.pop()
        y.pop()
    except IndexError:
        pass
    return x, y


## Takes an LOH dataframe and turns into LOH rectangles ##
def loh_to_rects(seg_df):
    x = []
    y = []
    for index, row in seg_df.iterrows():
        x.extend([row["START"], row["START"], row["END"], row["END"], None])
        y.extend([row["LOH"], 1 - row["LOH"], 1 - row["LOH"], row["LOH"], None])
    try:
        x.pop()
        y.pop()
    except IndexError:
        pass
    return x, y


## Takes the list of genes to plot and turn it into plottable text ##
def get_label_data(gene_list_df, label_fields: list):
    list_of_lists = []
    for index, row in gene_list_df.iterrows():
        labels_list = []
        for label in label_fields:
            labels_list.append(row[label])
        list_of_lists.append(labels_list)
    labels_array = np.array(list_of_lists)
    return labels_array


### Data-Processing Functions -- These process specific data in specific ways to run the CNV-LOH plots ###


## Takes the genes list dataframe and adds columns that will make plotting easier downstream ##
def add_required_columns(gene_list_df):
    gene_list_df["POSITION"] = (gene_list_df["START"] + gene_list_df["END"]) / 2
    gene_list_df["ZERO"] = 0
    log2_text = []
    loh_text = []
    focal_text = []
    for index, row in gene_list_df.iterrows():
        # if row['Focal'] == 'YES':
        #    if row['focal_type'] == 'EXON':
        #        focal_text.append('YES (Exon)')
        #        log2_text.append(float(row['focal_log2_gene']))
        #    else:
        #        focal_text.append('YES (Gene)')
        #        log2_text.append(float(row['focal_log2_exon']))
        # else:
        #    focal_text.append('NO')
        log2_text.append(float(row["seg_mean"]))
        if str(row["LOH"]) not in ["LOH", "NEARBY"]:
            loh_text.append("NO")
        else:
            loh_text.append(str(row["LOH"]))
    gene_list_df["LOG2_TEXT"] = log2_text
    gene_list_df["LOH_TEXT"] = loh_text
    # gene_list_df['FOCAL_TEXT']=focal_text
    return gene_list_df


## Takes in a genes table (Pandas dataframe) and returns a dictionary that will be used by fig_add_helper_plot() ##
def process_genes_toplot(genes_table):
    loh_nearby_genes = genes_table["LOH"] == "NEARBY"
    loh_genes = genes_table["LOH"] == "YES"
    all_gain_genes = genes_table["seg_mean"] > 0.25
    all_gain_genes_df = genes_table[all_gain_genes]
    weak_gain_genes = all_gain_genes_df["seg_mean"] < 0.45
    strong_gain_genes = all_gain_genes_df["seg_mean"] >= 0.45
    # gain_sensitive_genes = ['A' in a for a in genes_table['CGC_MUTATION_TYPES'].fillna('')]
    all_loss_genes = genes_table["seg_mean"] < -0.25
    all_loss_genes_df = genes_table[all_loss_genes]
    weak_loss_genes = all_loss_genes_df["seg_mean"] > 0.45
    strong_loss_genes = all_loss_genes_df["seg_mean"] <= -0.45
    # loss_sensitive_genes = ['D' in a for a in genes_table['CGC_MUTATION_TYPES'].fillna('')]

    # all_focal_genes = genes_table['Focal'] == 'YES'
    # all_focal_genes_df = genes_table[all_focal_genes]
    # focal_gain_genes = all_focal_genes_df['LOG2_TEXT'] > 0
    # focal_loss_genes = all_focal_genes_df['LOG2_TEXT'] < 0

    plotted_genes = []
    for i in range(0, len(list(loh_nearby_genes))):
        if True in [
            list(loh_nearby_genes)[i],
            list(loh_genes)[i],
            list(all_gain_genes)[i],
            list(all_loss_genes)[i],
        ]:  # ,list(all_focal_genes)[i]]:
            plotted_genes.append(True)
        else:
            plotted_genes.append(False)

    plotted_genes_df = genes_table[plotted_genes]

    genes_to_plot_dict = {
        "genes_table": genes_table,
        "loh_nearby_genes": loh_nearby_genes,
        "loh_genes": loh_genes,
        "all_gain_genes": all_gain_genes,
        "all_gain_genes_df": all_gain_genes_df,
        "weak_gain_genes": weak_gain_genes,
        "strong_gain_genes": strong_gain_genes,
        "gain_sensitive_genes": gain_sensitive_genes,
        "all_loss_genes": all_loss_genes,
        "all_loss_genes_df": all_loss_genes_df,
        "weak_loss_genes": weak_loss_genes,
        "strong_loss_genes": strong_loss_genes,
        "loss_sensitive_genes": loss_sensitive_genes,
        #'all_focal_genes':all_focal_genes,
        #'all_focal_genes_df':all_focal_genes_df,
        #'focal_gain_genes':focal_gain_genes,
        #'focal_loss_genes':focal_loss_genes,
        "plotted_genes": plotted_genes,
        "plotted_genes_df": plotted_genes_df,
    }
    return genes_to_plot_dict


## Divide up segments into bins for coverage.
def segment_to_bins(
    segment_length: int, bin_size: int = (10**6), segment_start: int = 0
):
    bin_starts = []
    bin_ends = []
    bin_start = 0
    while bin_start + bin_size <= segment_length:
        bin_starts.append(bin_start)
        bin_ends.append(bin_start + bin_size - 1)
        bin_start += bin_size
    bin_starts.append(bin_start)
    bin_ends.append(segment_length)

    # combine first two bins
    if len(bin_starts) >= 2:
        bin_starts.pop(1)
        bin_ends.pop(0)
    # combine last two bins
    if len(bin_starts) >= 2:
        bin_starts.pop(-1)
        bin_ends.pop(-2)
    bin_starts = list(np.array(bin_starts) + segment_start)
    bin_ends = list(np.array(bin_ends) + segment_start)
    return bin_starts, bin_ends


def find_point_bin(
    point_pos: int,
    num_bins: int,
    segment_start: int,
    segment_end: int,
    bin_size: int = 1000000,
):
    segment_len = segment_end - segment_start
    bin_index = (
        max(
            min(
                num_bins,
                math.ceil(
                    (
                        ((point_pos - segment_start - bin_size) / segment_len)
                        * (num_bins + 1)
                    )
                ),
            ),
            1,
        )
        - 1
    )
    return bin_index


## Adjusts LOH calls when the GATK results do not reflect reality
def add_segment_avg_baf(segment: dict, min_len: int = 1000000):
    print(
        f"{segment['position']['chrom']}:{segment['position']['start']}-{segment['position']['end']}"
    )
    seg_len = segment["position"]["length"]
    seg_loh = segment["loh"]
    gatk_baf = seg_loh["loh_allele_fraction"]
    loh_vars = seg_loh["loh_supporting_variants"]
    n_vars = len(loh_vars)
    bafs = []
    for i in loh_vars:
        bafs.append(i["b_allele_freq"])
    if seg_len >= min_len or len(loh_vars) > 0:
        var_per_mb = n_vars / (seg_len / (10**6))
    else:
        var_per_mb = None
    if n_vars == 0:
        if seg_len >= min_len:
            avg_baf = 0
        else:
            avg_baf = None
    else:
        avg_baf = sum(bafs) / n_vars
    print(f"\tGATK-Calculated BAF: {gatk_baf}")
    print(f"\tBAF from Variants: {avg_baf}")
    print(f"\tVariants Per Megabase: {var_per_mb}")
    segment["loh"]["vars_per_mb"] = var_per_mb
    segment["loh"]["baf_from_vars"] = avg_baf
    return segment


def downweight_excessive_vars(
    segment: dict, highest_decile: float = 20.0, low_cutoff: float = -0.10
):
    allele_bafs = []
    allele_weights = []

    variants = segment["loh"]["loh_supporting_variants"]

    dup_vars = 0
    excl_vars = 0

    if len(variants) == 0:
        allele_bafs = [0]
        allele_weights = [1]
    elif len(variants) == 1:
        allele_bafs = [variants[0]["b_allele_freq"]]
        allele_weights = [1]
    else:
        prev_variant = variants.pop(0)
        distance = abs(prev_variant["position"] - segment["position"]["start"])
        print(distance)
        allele_bafs.append(prev_variant["b_allele_freq"])
        allele_weights.append(min(highest_decile * 100 / ((10**6) / distance), 1))
        for i in variants:
            distance = min(
                abs(prev_variant["position"] - i["position"]),
                abs(segment["position"]["end"] - i["position"]),
            )
            # print(distance)
            if i["b_allele_freq"] >= low_cutoff:
                allele_bafs.append(i["b_allele_freq"])
                allele_weights.append(
                    min(highest_decile * 100 / ((10**6) / distance), 1)
                )
                ## Double variants that are at 1 - cutoff to make up for the removal of the low-frequency alleles.
                # In normal non-degraded samples, the number below cutoff and above 1-cutoff should be nearly the same.
                if i["b_allele_freq"] >= (1 - low_cutoff):
                    allele_bafs.append(i["b_allele_freq"])
                    allele_weights.append(
                        min(highest_decile * 100 / ((10**6) / distance), 1)
                    )
                    dup_vars += 1
                prev_variant = i
            else:
                excl_vars += 1

    print(f"Dup. vars: {dup_vars} | Excl. vars: {excl_vars}")
    # print(allele_weights)

    return allele_bafs, allele_weights


def loh_rescue(segments_dict: dict, metadata_dict: dict = None, min_len: int = 5000000):
    print("Running LOH rescue pass.\n")
    for i in range(0, len(segments_dict)):
        segments_dict[i] = add_segment_avg_baf(segments_dict[i], min_len)

    vars_per_mb = sorted(
        [
            x["loh"]["vars_per_mb"]
            for x in segments_dict
            if x["loh"]["vars_per_mb"] is not None
        ]
    )
    avg_vars_per_mb = sum(vars_per_mb) / len(vars_per_mb)
    lowest_decile = vars_per_mb[int(len(vars_per_mb) * 0.1)] / 2
    highest_decile = vars_per_mb[int(len(vars_per_mb) * 0.9)]
    median_vars_per_mb = vars_per_mb[int(len(vars_per_mb) * 0.5)]

    print(
        f"Median variants per megabase: {median_vars_per_mb}\nLowest decile variants per megabase: {lowest_decile}"
    )

    if metadata_dict is not None:
        metadata_dict["median_vars_per_mb"] = median_vars_per_mb

    if abs(lowest_decile - median_vars_per_mb) >= (median_vars_per_mb * 0.33):
        fill_homoz = True
        print(
            "Lowest decile of variants/megabase deviates greatly from median. Filling in homozygotes."
        )
    else:
        fill_homoz = False
    corrections = []
    for i in range(0, len(segments_dict)):
        label = f"{segments_dict[i]['position']['chrom']}:{segments_dict[i]['position']['start']}-{segments_dict[i]['position']['end']}"

        # print(segments_dict[i]['loh']['loh_allele_fraction'] == 'nan')

        allele_bafs, allele_weights = downweight_excessive_vars(
            segments_dict[i], highest_decile
        )
        weighted_sum_baf = 0
        total_weight = sum(allele_weights)  # no lowest-decile value here
        for v in range(0, len(allele_bafs)):
            weighted_sum_baf += allele_bafs[v] * allele_weights[v]

        adj_vars_per_mb = total_weight / (
            segments_dict[i]["position"]["length"] / (10**6)
        )
        print(
            f"Variants per mb: {segments_dict[i]['loh']['vars_per_mb']} | Adjusted: {adj_vars_per_mb}"
        )

        if segments_dict[i]["loh"]["vars_per_mb"] is None:
            print(i)
            print(
                f"\tSegment {label} has no variants but was too short to fill ({segments_dict[i]['position']['length']})."
            )

        elif (
            (
                segments_dict[i]["loh"]["vars_per_mb"] < lowest_decile
                or adj_vars_per_mb < lowest_decile
            )
            and (segments_dict[i]["position"]["length"] > min_len)
            and fill_homoz is True
        ):
            print(label)
            # segments_dict[i]['loh']['adj_baf'] = (segments_dict[i]['loh']['baf_from_vars']*segments_dict[i]['loh']['vars_per_mb']) / lowest_decile

            segments_dict[i]["loh"]["adj_baf"] = weighted_sum_baf / (
                lowest_decile * (segments_dict[i]["position"]["length"] / (10**6))
            )

        elif (
            str(segments_dict[i]["loh"]["loh_allele_fraction"]).lower() == "nan"
            or segments_dict[i]["loh"]["loh_allele_fraction"] is np.nan
            or (
                segments_dict[i]["loh"]["loh_allele_fraction"]
                - (weighted_sum_baf / total_weight)
            )
            > 0.1
        ):  # Add an abs() to the second part if you want to be able to reduce LOH also.
            print(label)

            segments_dict[i]["loh"]["adj_baf"] = weighted_sum_baf / total_weight
            print(segments_dict[i]["loh"]["adj_baf"])

        if "adj_baf" in segments_dict[i]["loh"]:
            label = f"{label} | {segments_dict[i]['loh']['loh_allele_fraction']} -> {segments_dict[i]['loh']['adj_baf']}"
            corrections.append(label)
    print(f"Segments corrected: {len(corrections)}")
    for i in corrections:
        print("\t" + i)

    return segments_dict, metadata_dict


## Takes the CNVLOH json file path and creates 3 dataframes in a dictionary for plotting ##
def cnvloh_json2dfs(
    cnvloh_dict: str, chromosome: str = None, coverage_bin_size: int = 1000000
):
    segment_chrom = []
    segment_index = []
    segment_id = []
    segment_pos_start = []
    segment_pos_end = []
    segment_cnv = []
    segment_loh = []
    segment_cnv_marks = []
    segment_cnv_reads = []
    segment_cnv_stdev = []
    segment_loh_marks = []
    segment_loh_reads = []
    segment_num_bins = []

    subseg_chrom = []
    subseg_x0 = []
    subseg_x1 = []
    subseg_mid = []
    subseg_len = []
    subseg_half_len = []
    subseg_reads = []
    subseg_y = []
    subseg_bin = []
    subseg_parent_index = []

    var_chrom = []
    var_x = []
    var_y = []
    ref_base = []
    alt_base = []
    ref_reads = []
    alt_reads = []

    if chromosome:
        chroms_list = [chromosome]
    else:
        chroms_set = set()
        for i in cnvloh_dict:
            chroms_set.add(i["position"]["chrom"])
        chroms_list = list(chroms_set)

    for chromosome in chroms_list:
        for i in range(0, len(cnvloh_dict)):
            if cnvloh_dict[i]["position"]["chrom"] == chromosome:
                seg_index = cnvloh_dict[i]["position"]["seg_index"]
                seg_id = cnvloh_dict[i]["position"]["seg_id"]
                segment_chrom.append(cnvloh_dict[i]["position"]["chrom"])
                segment_index.append(seg_index)
                segment_id.append(seg_id)
                segment_pos_start.append(cnvloh_dict[i]["position"]["start"])
                segment_pos_end.append(cnvloh_dict[i]["position"]["end"])
                segment_cnv.append(cnvloh_dict[i]["cnv"]["log2_copy_ratio"])
                segment_cnv_marks.append(cnvloh_dict[i]["cnv"]["cnv_supporting_points"])
                segment_cnv_reads.append(cnvloh_dict[i]["cnv"]["cnv_supporting_reads"])
                segment_cnv_stdev.append(cnvloh_dict[i]["cnv"]["log2_stdev"])
                calc_seg_len = (
                    cnvloh_dict[i]["position"]["end"]
                    - cnvloh_dict[i]["position"]["start"]
                )
                segment_bin_starts, segment_bin_ends = segment_to_bins(
                    calc_seg_len, coverage_bin_size, cnvloh_dict[i]["position"]["start"]
                )
                calc_num_bins = len(segment_bin_starts)
                segment_num_bins.append(calc_num_bins)
                segment_points = []
                segment_reads = 0
                calc_bin = 0
                for j in range(
                    0, len(cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"])
                ):
                    calc_midpoint = round(
                        (
                            cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                                "start"
                            ]
                            + cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                                "end"
                            ]
                        )
                        / 2
                    )
                    calc_sub_len = (
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["end"]
                        - cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                            "start"
                        ]
                    )
                    calc_bin = find_point_bin(
                        calc_midpoint,
                        calc_num_bins,
                        cnvloh_dict[i]["position"]["start"],
                        cnvloh_dict[i]["position"]["end"],
                        coverage_bin_size,
                    )
                    subseg_parent_index.append(seg_index)
                    subseg_chrom.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["chrom"]
                    )
                    subseg_x0.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["start"]
                    )
                    subseg_x1.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["end"]
                    )
                    subseg_mid.append(calc_midpoint)
                    subseg_len.append(calc_sub_len)
                    subseg_half_len.append(
                        (
                            cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                                "end"
                            ]
                            - cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                                "start"
                            ]
                        )
                        / 2
                    )
                    subseg_y.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                            "log2_copy_ratio"
                        ]
                    )
                    subseg_reads.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["reads"]
                    )
                    if cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j]["reads"]:
                        segment_reads += cnvloh_dict[i]["cnv"][
                            "cnv_supporting_subsegments"
                        ][j]["reads"]
                    segment_points.append(
                        cnvloh_dict[i]["cnv"]["cnv_supporting_subsegments"][j][
                            "log2_copy_ratio"
                        ]
                    )
                    subseg_bin.append(calc_bin)

                if "loh" in cnvloh_dict[i]:
                    loh_supporting_points = cnvloh_dict[i]["loh"][
                        "loh_supporting_points"
                    ]
                    loh_supporting_reads = cnvloh_dict[i]["loh"]["loh_supporting_reads"]
                    loh_total_reads = cnvloh_dict[i]["loh"]["loh_supporting_reads"]
                    segment_loh_marks.append(loh_supporting_points)
                    segment_loh_reads.append(loh_supporting_reads)

                    if "adj_baf" in cnvloh_dict[i]["loh"]:
                        print("Adjusted segment BAF found. Substituting.")
                        segment_loh.append(cnvloh_dict[i]["loh"]["adj_baf"])
                    elif (
                        loh_supporting_points > 0
                        and "loh_allele_fraction" in cnvloh_dict[i]["loh"]
                    ):
                        segment_loh.append(cnvloh_dict[i]["loh"]["loh_allele_fraction"])
                    else:
                        # If there are no supporting points, we just insert a dummy exactly heterozygous value
                        # since we need the lists to be the same length for the upcoming zip()
                        segment_loh.append(0.5)

                    # Add supporting point data if available
                    if loh_supporting_points > 0:
                        for j in range(
                            0, len(cnvloh_dict[i]["loh"]["loh_supporting_variants"])
                        ):
                            var_chrom.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "chrom"
                                ]
                            )
                            var_x.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "position"
                                ]
                            )
                            var_y.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "alt_vaf"
                                ]
                            )
                            ref_reads.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "ref_count"
                                ]
                            )
                            ref_base.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "ref_allele"
                                ]
                            )
                            alt_reads.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "alt_count"
                                ]
                            )
                            alt_base.append(
                                cnvloh_dict[i]["loh"]["loh_supporting_variants"][j][
                                    "alt_allele"
                                ]
                            )
                else:
                    pass

    seg_df = pd.DataFrame(
        list(
            zip(
                segment_chrom,
                segment_index,
                segment_id,
                segment_pos_start,
                segment_pos_end,
                segment_cnv,
                segment_cnv_marks,
                segment_cnv_reads,
                segment_cnv_stdev,
                segment_loh,
                segment_loh_marks,
                segment_loh_reads,
                segment_num_bins,
            )
        ),
        columns=[
            "CHROM",
            "INDEX",
            "ID",
            "START",
            "END",
            "LOG2",
            "LOG2_MARKS",
            "LOG2_READS",
            "LOG2_STDEV",
            "LOH",
            "LOH_MARKS",
            "LOH_READS",
            "NUM_BINS",
        ],
    )
    seg_df = seg_df.sort_values("INDEX").reset_index(drop=True)
    cnv_df = pd.DataFrame(
        list(
            zip(
                subseg_chrom,
                subseg_parent_index,
                subseg_x0,
                subseg_x1,
                subseg_mid,
                subseg_len,
                subseg_half_len,
                subseg_reads,
                subseg_y,
                subseg_bin,
            )
        ),
        columns=[
            "CHROM",
            "PARENT_INDEX",
            "START",
            "END",
            "MIDPOINT",
            "LEN",
            "HALF_LEN",
            "READS",
            "LOG2",
            "BIN",
        ],
    )
    cnv_df = cnv_df.sort_values("PARENT_INDEX").reset_index(drop=True)
    try:
        loh_df = pd.DataFrame(
            list(
                zip(var_chrom, var_x, var_y, ref_base, alt_base, ref_reads, alt_reads)
            ),
            columns=["CHROM", "POS", "VAF", "REF", "ALT", "REF_READS", "ALT_READS"],
        )
    except:
        loh_df = None
    df_dict = {"seg_df": seg_df, "cnv_df": cnv_df, "loh_df": loh_df}
    df_dict["cnv_bins_df"] = calc_coverage_bins(seg_df, cnv_df, coverage_bin_size)
    return df_dict


def calc_coverage_bins(
    seg_df,
    cnv_df,
    bin_width: int = 1000000,
    min_support: int = 2,
    abs_positions: bool = False,
):
    bins_dict = {}
    all_chroms = []
    all_segs = []
    all_bins = []
    all_starts = []
    all_ends = []
    all_points = []
    all_means = []
    all_stdevs = []
    all_reads = []
    for i in set(seg_df["CHROM"]):
        chrom_df = seg_df[seg_df["CHROM"] == i]
        chrom_cnv_df = cnv_df[cnv_df["CHROM"] == i]
        bins_dict[i] = {}
        for j in set(chrom_df["INDEX"]):
            index_seg_df = chrom_df[chrom_df["INDEX"] == j]
            index_cnv_df = chrom_cnv_df[chrom_cnv_df["PARENT_INDEX"] == j]
            seg_start = min(index_seg_df["START"])
            seg_end = max(index_seg_df["END"])
            bins_dict[i][j] = {}
            segment_bin_starts, segment_bin_ends = segment_to_bins(
                seg_end - seg_start, bin_width, seg_start
            )
            bins = set(index_cnv_df[index_cnv_df["PARENT_INDEX"] == j]["BIN"])
            for k in bins:
                bins_df = index_cnv_df[index_cnv_df["BIN"] == k]
                log2s = list(bins_df["LOG2"])
                read_counts = list(bins_df["READS"])
                for l in range(0, len(log2s)):
                    point_log2 = log2s[l]
                    point_reads = read_counts[l]
                    if k not in bins_dict[i][j]:
                        bins_dict[i][j][k] = {
                            "start": int(segment_bin_starts[int(k)]),
                            "end": int(segment_bin_ends[int(k)]),
                            "log2_scores": [point_log2],
                            "reads": point_reads,
                        }
                    else:
                        bins_dict[i][j][k]["log2_scores"].append(point_log2)
                        if point_reads:
                            bins_dict[i][j][k]["reads"] += point_reads
                bins_dict[i][j][k]["points"] = len(bins_dict[i][j][k]["log2_scores"])
                bins_dict[i][j][k]["mean"] = sum(
                    bins_dict[i][j][k]["log2_scores"]
                ) / len(bins_dict[i][j][k]["log2_scores"])
                bins_dict[i][j][k]["stdev"] = np.nan_to_num(
                    np.std(bins_dict[i][j][k]["log2_scores"])
                )
                if bins_dict[i][j][k]["points"] > min_support:
                    all_chroms.append(i)
                    all_segs.append(j)
                    all_bins.append(k)
                    all_starts.append(bins_dict[i][j][k]["start"])
                    all_ends.append(bins_dict[i][j][k]["end"])
                    all_points.append(bins_dict[i][j][k]["points"])
                    all_means.append(bins_dict[i][j][k]["mean"])
                    all_stdevs.append(bins_dict[i][j][k]["stdev"])
                    all_reads.append(bins_dict[i][j][k]["reads"])
    cnv_bins_df = pd.DataFrame(
        list(
            zip(
                all_chroms,
                all_segs,
                all_bins,
                all_starts,
                all_ends,
                all_points,
                all_reads,
                all_means,
                all_stdevs,
            )
        ),
        columns=[
            "CHROM",
            "SEG",
            "BIN",
            "START",
            "END",
            "POINTS",
            "READS",
            "MEAN",
            "STDEV",
        ],
    )
    return cnv_bins_df


### Figure-Plotting -- These are used to actually make the CNV-LOH plots ###


## Sets up the three-subplot stacked canvass for CNV, Helper, and LOH ##
def setup_cnvloh_plot(genome_info, chromosome, include_helper=True):
    chrom_info = genome_info[chromosome]

    if include_helper:
        plot_rows = 3
    else:
        plot_rows = 2
    fig = make_subplots(
        rows=plot_rows, cols=1, shared_xaxes=True, vertical_spacing=0.005
    )

    chromosome_length = chrom_info["length"].astype(int)
    centromere_start = chrom_info["centromere_start"].astype(int)
    centromere_end = chrom_info["centromere_end"].astype(int)

    fig.add_vrect(
        x0=centromere_start,
        x1=centromere_end,
        fillcolor="Grey",
        opacity=0.5,
        layer="below",
        line_width=0,
        row=1,
        col=1,
    )
    boundary_df = pd.DataFrame(
        list(zip([0, chromosome_length], [-1, 1])), columns=["X", "Y"]
    )
    fig.add_trace(
        go.Scatter(
            x=list(boundary_df["X"]),
            y=list(boundary_df["Y"]),
            mode="markers",
            opacity=0,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(
            x=[0, chromosome_length],
            y=[0, 0],
            fill="toself",
            mode="lines",
            name="Zero Point",
        ),
        row=1,
        col=1,
    )
    fig["data"][-1]["line"]["color"] = "#000000"

    if include_helper:
        fig.add_vrect(
            x0=centromere_start,
            x1=centromere_end,
            fillcolor="Grey",
            opacity=0.5,
            layer="below",
            line_width=0,
            row=2,
            col=1,
        )
        boundary_df = pd.DataFrame(
            list(zip([0, chromosome_length], [-1, 1])), columns=["X", "Y"]
        )
        fig.add_trace(
            go.Scatter(
                x=list(boundary_df["X"]),
                y=list(boundary_df["Y"]),
                mode="markers",
                opacity=0,
            ),
            row=2,
            col=1,
        )

    fig.add_vrect(
        x0=centromere_start,
        x1=centromere_end,
        fillcolor="Grey",
        opacity=0.5,
        layer="below",
        line_width=0,
        row=plot_rows,
        col=1,
    )
    boundary_df = pd.DataFrame(
        list(zip([0, chromosome_length], [0, 1])), columns=["X", "Y"]
    )
    fig.add_trace(
        go.Scatter(
            x=list(boundary_df["X"]),
            y=list(boundary_df["Y"]),
            mode="markers",
            opacity=0,
        ),
        row=plot_rows,
        col=1,
    )

    fig.update_layout(
        height=400 * plot_rows, width=1000, title_text=f"CNV-LOH: {chromosome}"
    )

    return fig


## Adds the top CNV plot to the three-canvas 'fig' ##
def fig_add_cnv_plot(fig, seg_df, cnv_df, plot_row=1):
    plot_x, plot_y = cnv_to_segs(cnv_df)
    fig.add_trace(
        go.Scatter(
            x=plot_x,
            y=plot_y,
            fill="toself",
            mode="lines+markers",
            name="CNV Supporting Marks",
        ),
        row=plot_row,
        col=1,
    )
    fig["data"][-1]["line"]["color"] = "#329efb"

    plot_x, plot_y = cnv_to_segs(seg_df)
    fig.add_trace(
        go.Scatter(
            x=plot_x, y=plot_y, fill="toself", mode="lines", name="CNV Called Segments"
        ),
        row=plot_row,
        col=1,
    )
    fig["data"][-1]["line"]["color"] = "#f63729"
    return None


## Adds the middle Helper plot to the three-canvas 'fig' ##
def fig_add_helper_plot(fig, genes_to_plot_dict, plot_row=2):
    # Plot Gains #
    # Strong
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_gain_genes_df"][
            genes_to_plot_dict["strong_gain_genes"]
        ],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x, y=plot_y, fill="toself", mode="lines", name="Copy Gain"
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#329efb"
    # Weak
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_gain_genes_df"][genes_to_plot_dict["weak_gain_genes"]],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x, y=plot_y, fill="toself", mode="lines", name="Copy Gain (Weak)"
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#164c90"
    # Plot Losses #
    # Strong
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_loss_genes_df"][
            genes_to_plot_dict["strong_loss_genes"]
        ],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x, y=plot_y, fill="toself", mode="lines", name="Copy Loss"
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#f63729"
    # Weak
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_loss_genes_df"][genes_to_plot_dict["weak_loss_genes"]],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x, y=plot_y, fill="toself", mode="lines", name="Copy Loss (Weak)"
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#ab261d"
    # Plot Focal Events #
    # Gain
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_focal_genes_df"][
            genes_to_plot_dict["focal_gain_genes"]
        ],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x,
                y=plot_y,
                fill="toself",
                mode="lines",
                name="Copy Gain (Focal)",
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#0431f9"
    # Loss
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["all_focal_genes_df"][
            genes_to_plot_dict["focal_loss_genes"]
        ],
        to=1,
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x,
                y=plot_y,
                fill="toself",
                mode="lines",
                name="Copy Loss (Focal)",
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#f62400"
    # Plot LOH #
    # LOH
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["genes_table"][genes_to_plot_dict["loh_genes"]], to=-1
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(x=plot_x, y=plot_y, fill="toself", mode="lines", name="LOH"),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#fbc300"
    # Nearby LOH
    plot_x, plot_y = genes_to_rects(
        genes_to_plot_dict["genes_table"][genes_to_plot_dict["loh_nearby_genes"]], to=-1
    )
    if len(plot_x) > 0:
        fig.add_trace(
            go.Scatter(
                x=plot_x, y=plot_y, fill="toself", mode="lines", name="LOH (Nearby)"
            ),
            row=plot_row,
            col=1,
        )
        fig["data"][-1]["line"]["color"] = "#ca9c00"

    # If we've filtered all the genes out from being plotted, skip this block since it'll fail
    if not genes_to_plot_dict["plotted_genes_df"].empty:
        fig.add_trace(
            go.Scatter(
                x=list(genes_to_plot_dict["plotted_genes_df"]["POSITION"]),
                y=[
                    math.sin(i) * 0.9
                    for i in range(
                        0, len(list(genes_to_plot_dict["plotted_genes_df"]["POSITION"]))
                    )
                ],
                mode="markers+text",
                name="Affected Genes",
                text=list(genes_to_plot_dict["plotted_genes_df"]["SYMBOL"]),
                textposition="top center",
            ),
            row=plot_row,
            col=1,
        )

    # fig1.add_trace(px.scatter(plotted_genes_df,x="POSITION",y="ZERO",hover_data=['SYMBOL','START','END','LOG2_TEXT','LOH_TEXT','FOCAL_TEXT']))
    # fig1['data'][0]['hovertemplate']='<b>%{customdata[0]}</b><br>(%{customdata[1]}-%{customdata[2]})<br>Log2: %{customdata[3]}<br>LOH: %{customdata[4]}<br>Focal: %{customdata[5]}'
    fig["data"][-1]["marker"]["color"] = "#000000"
    return None


## Adds the bottom LOH plot to the three-plot 'fig' canvas
def fig_add_loh_plot(fig, seg_df, loh_df, plot_row=3):
    # LOH Plot
    plot_x, plot_y = loh_to_rects(seg_df)
    fig.add_trace(
        go.Scatter(x=plot_x, y=plot_y, fill="toself", mode="lines", name="LOH Regions"),
        row=plot_row,
        col=1,
    )
    fig["data"][-1]["line"]["color"] = "#fbc300"

    fig.add_trace(
        go.Scatter(
            x=list(loh_df["POS"]),
            y=list(loh_df["VAF"]),
            mode="markers",
            name="Supporting Variants",
        ),
        row=plot_row,
        col=1,
    )
    # fig2 = px.scatter(loh_df,x="POS",y="VAF",hover_data=['CHROM','POS','VAF','REF','ALT','REF_READS','ALT_READS'])
    # fig2['data'][0]['hovertemplate']='<b>%{customdata[0]}:%{x} %{customdata[1]}>%{customdata[2]}</b><br>VAF: %{y}<br>%{customdata[1]}: %{customdata[3]}<br>%{customdata[2]}: %{customdata[4]}'
    fig["data"][-1]["marker"]["color"] = "#000000"
    return None


## Make a CNV-LOH plot based on a chromosome, genome info, JSON, and genes table ##
# def plot_chrom_html_plotly(chromosome, genome_info, cnvloh_dict, out_path, genes_table=None, genes_list=None):
#     plot_path = os.path.join(out_path,f'{chromosome}_CNVLOH.html')

#     if genes_table is not None:
#         genes_table = add_required_columns(genes_table)

#     if not chromosome.startswith("chr"):
#         print("Prefixing chromosome with 'chr'")
#         chromosome = "chr" + chromosome

#     print(f'Plotting {chromosome} -> {plot_path}')

#     df_dict = cnvloh_json2dfs(cnvloh_dict,chromosome)
#     genes_table = subset_genes_table(genes_table,genes_list)
#     try:
#         genes_table_chrom = genes_table[genes_table['CHR'] == chromosome]
#         genes_to_plot_dict = process_genes_toplot(genes_table_chrom)
#         helper_plot = True
#     except:
#         helper_plot = False
#     fig = setup_cnvloh_plot(genome_info,chromosome, include_helper=helper_plot)
#     fig_add_cnv_plot(fig,df_dict['seg_df'],df_dict['cnv_df'])
#     if helper_plot:
#         fig_add_helper_plot(fig,genes_to_plot_dict)
#         loh_row = 3
#     else:
#         loh_row = 2
#     fig_add_loh_plot(fig,df_dict['seg_df'],df_dict['loh_df'], loh_row)
#     # fig.show()
#     fig.write_html(plot_path)
#     return None


def plot_chrom_html_bokeh(
    chromosome,
    genome_info,
    cnvloh_dict,
    out_path,
    genes_table=None,
    genes_list=None,
    focal_events=None,
    output_type: str = "json",
    cnvloh_dict_2: dict = None,
    focals_dict: dict = None,
    metadata: dict = None,
):
    plot_path = os.path.join(out_path, f"{chromosome}_CNVLOH.{output_type}")

    genome_info = order_genome_info(genome_info)

    if not chromosome.startswith("chr"):
        print("\tPrefixing chromosome with 'chr'")
        chromosome = "chr" + chromosome

    if (genes_table) is not None:
        genes_table_chrom = genes_table[genes_table["CHR"] == chromosome].copy()
        genes_table_chrom = add_required_columns(genes_table)
    else:
        genes_table_chrom = None

    print(f"Plotting {chromosome} -> {plot_path}")
    print("*** DAVE TEST plot_chrom_html_bokeh ***")
    print("*** GRANT TEST plot_chrom_html_bokeh ***")
    if genes_table is not None:
        print("genes_table columns:")
        print(",".join(genes_table.keys()))

    df_dict = cnvloh_json2dfs(cnvloh_dict, chromosome)
    if cnvloh_dict_2:
        df_dict_2 = cnvloh_json2dfs(cnvloh_dict_2, chromosome)
    else:
        df_dict_2 = None
    # genes_table = subset_genes_table(genes_table,genes_list)

    if metadata is not None and "sex" in metadata and metadata["sex"] is not None:
        if metadata["sex"].lower() == "male":
            print("Patient is male")
            is_male = True
        else:
            print("Patient is female")
            is_male = False
    else:
        print("Patient sex unknown. Defaulting to female.")
        is_male = False

    try:
        cnvloh_bokeh_plot(
            genes_table=genes_table_chrom,
            genes_list=genes_list,
            cnvloh_dna_dfs=df_dict,
            chrom_info_dict=genome_info[chromosome],
            chromosome=chromosome,
            outpath=plot_path,
            output_type=output_type,
            cnvloh_dna_dfs_2=df_dict_2,
            focals_dict=focals_dict,
            male_adjustment=is_male,
        )
    except Exception as e:
        print("Exception: " + str(e))
        traceback.print_exc()
    return None


def plot_genome_html_bokeh(
    genome_info,
    cnvloh_dict,
    out_path,
    genes_table=None,
    genes_list=None,
    focal_events=None,
    output_type: str = "json",
    cnvloh_dict_2: dict = None,
    focals_dict: dict = None,
    metadata: dict = None,
):
    plot_path = os.path.join(out_path, f"GenomeView_CNVLOH.{output_type}")

    genome_info = order_genome_info(genome_info)

    try:
        genes_table = add_required_columns(genes_table)
    except:
        genes_table = None

    print(f"Plotting Genome View -> {plot_path}")

    df_dict = cnvloh_json2dfs(cnvloh_dict, coverage_bin_size=5000000)
    # genes_table = subset_genes_table(genes_table,genes_list)

    if cnvloh_dict_2:
        df_dict_2 = cnvloh_json2dfs(cnvloh_dict_2, coverage_bin_size=5000000)
    else:
        df_dict_2 = None

    if metadata is not None and "sex" in metadata and metadata["sex"] is not None:
        if metadata["sex"].lower() == "male":
            print("Patient is male")
            is_male = True
        else:
            print("Patient is female")
            is_male = False
    else:
        print("Patient sex unknown. Defaulting to female.")
        is_male = False

    cnvloh_bokeh_plot_genome(
        genes_table=genes_table,
        cnvloh_dna_dfs=df_dict,
        genome_info=genome_info,
        genes_list=genes_list,
        outpath=plot_path,
        output_type=output_type,
        cnvloh_dna_dfs_2=df_dict_2,
        focals_dict=focals_dict,
        male_adjustment=is_male,
    )

    return None


## Go through and plot all the chromosomes in a genome info file ##
def plot_all_chrom_html(
    genome_info,
    cnvloh_dict,
    out_path,
    genes_table=None,
    genes_list=None,
    output_type: str = "json",
    cnvloh_dict_2: dict = None,
    focals_dict: dict = None,
    metadata: dict = None,
):
    # genes_table = subset_genes_table(genes_table,genes_list)
    time_start = time.time()
    plot_genome_html_bokeh(
        out_path=out_path,
        genome_info=genome_info,
        cnvloh_dict=cnvloh_dict,
        genes_table=genes_table,
        genes_list=genes_list,
        output_type=output_type,
        cnvloh_dict_2=cnvloh_dict_2,
        focals_dict=focals_dict,
        metadata=metadata,
    )
    time_elapsed = time.time() - time_start
    print(f"\tGenome plotted in: {time_elapsed:.2f} seconds.")
    for i in genome_info:
        time_start = time.time()
        if (genes_table) is not None:
            genes_table_chrom = genes_table[genes_table["CHR"] == i].copy()
        else:
            genes_table_chrom = None
        plot_chrom_html_bokeh(
            chromosome=i,
            out_path=out_path,
            genome_info=genome_info,
            cnvloh_dict=cnvloh_dict,
            genes_table=genes_table_chrom,
            genes_list=genes_list,
            output_type=output_type,
            cnvloh_dict_2=cnvloh_dict_2,
            focals_dict=focals_dict,
            metadata=metadata,
        )
        time_elapsed = time.time() - time_start
        print(f"\tChromosome {i} plotted in: {time_elapsed:.2f} seconds.")


#### Arg-Parser -- This will only be used for the 'name=main' direct call of the script ####


def run_cnvplot_argparser():
    parser = argparse.ArgumentParser(
        description="Creates interactive plots for CNV-LOH."
    )
    parser.add_argument(
        "--sample-name",
        help="(Optional) Name of CNV-LOH sample being processed.",
        default="cnvloh_sample",
    )
    parser.add_argument(
        "--chromosome",
        help="(Optional) Chromosome to plot. Or, enter 'genome' for genome-wide plot. If not specified, plots all chromosomes.",
    )
    parser.add_argument(
        "--genes-table",
        help="(Optional) Path to genes table file, for annotation plotting.",
    )
    parser.add_argument(
        "--gene-sublist",
        help="Path to text file with subset of genes to display. If not given, will plot all genes.",
    )
    parser.add_argument(
        "--genome-info",
        help="Path to file showing chromosome lengths and centromere positions.",
        required=True,
    )
    parser.add_argument(
        "--cnvloh-dir", help="Directory to output CNV-LOH plots.", required=True
    )
    parser.add_argument(
        "--cnvloh-json",
        help="Path to CNV-LOH JSON. If not specified, must instead use all three of Counts, Hets, and Modelsegs.",
    )
    parser.add_argument(
        "--focals-json",
        help="(Optional) Path to JSON containing gene and exon focal CNV events.",
        default=None,
    )
    parser.add_argument("--hets-file", help="Path to CNV-LOH *.hets.tsv File")
    parser.add_argument(
        "--counts-file", help="Path to CNV-LOH *.counts.tsv.denoisedCR.tsv"
    )
    parser.add_argument(
        "--modelsegs-file", help="Path to CNV-LOH *.modelFinal.seg file"
    )
    parser.add_argument(
        "--sample-type",
        help="(Optional) Sample type -- typically Normal, Tumor, or Tumor_Normal",
        default=None,
    )
    parser.add_argument(
        "--case-name",
        help="(Optional) Name of the case to which the sample belongs.",
        default=None,
    )
    parser.add_argument(
        "--out-format",
        help="(Optional) Output format.",
        default="html",
        choices=["html", "json"],
    )
    parser.add_argument(
        "--cnvloh-json-2",
        help="Path to CNV-LOH JSON, for optional comparator sample.",
        default=None,
    )
    return parser


#### Main Function ####
# May need to put some of this code in functions too, to make it easier to call from another script #

if __name__ == "__main__":
    time_start = time.time()
    args = run_cnvplot_argparser().parse_args()
    # Use the cnvloh-dir argument to figure out the relevant file paths and other info.
    output_dir = args.cnvloh_dir
    cnvloh_json = args.cnvloh_json
    cnvloh_json_2 = args.cnvloh_json_2
    out_format = args.out_format
    if not cnvloh_json:
        from .cnvloh2json import generate_cnvloh_json

        modeled_seg_file = args.modelsegs_file
        cnv_subseg_file = args.counts_file
        loh_calls_file = args.hets_file
        if not modeled_seg_file or not cnv_subseg_file or not loh_calls_file:
            raise Exception(
                "Insufficient data provided to begin plotting. Must specify either a CNV-LOH JSON, or all three of ModeledSeg, Counts, and Hets."
            )
        else:
            sample_name = args.sample_name
            case_name = args.case_name
            cnv_params = None
            data_type = "DNA"
            software = "GATK"
            software_version = "4.2.2.0"
            sample_type = args.sample_type
            cnvloh_json = os.path.join(
                output_dir, f"cnvloh_output_{data_type}-{software}.json"
            )
            generate_cnvloh_json(
                output_dir,
                sample_name,
                modeled_seg_file,
                cnv_subseg_file,
                loh_calls_file,
                data_type,
                case_name,
                cnv_params,
                software,
                software_version,
                sample_type,
            )

    # Get files needed for plotting #
    # Read in the genes table file as a df
    if not args.genes_table:
        genes_table = None
        gene_sublist = None
    else:
        genes_table_file = args.genes_table
        if os.path.exists(genes_table_file):
            genes_table = pd.read_csv(genes_table_file, sep=",", header=0)
            # Sublist file contains the subset of genes to plot
            gene_sublist_file = args.gene_sublist
            if gene_sublist_file:
                if os.path.exists(gene_sublist_file):
                    try:
                        gene_sublist = list(
                            pd.read_csv(gene_sublist_file, header=None)[0].astype(str)
                        )
                    except:
                        print("Gene list is an empty file. Disabling gene plotting.")
                        gene_sublist = []
                else:
                    raise Exception(
                        "A gene list file was specified but the file does not exist."
                    )
            else:
                gene_sublist = None
        else:
            print("No annotations table provided. Plotting in non-annotated mode.")
            genes_table = None
            gene_sublist = None

    # Genome info contains chromosomes, chrom lengths, and centromere positions
    genome_info_file = args.genome_info
    with open(genome_info_file, "r") as genome_info_data:
        genome_info_raw = json.loads(genome_info_data.read())
    genome_info = order_genome_info(genome_info_raw)

    # Make a directory for the HTML files #
    interactive_plots_dir = os.path.join(output_dir, "Interactive_Plots_HTML")
    if not os.path.exists(interactive_plots_dir):
        os.mkdir(interactive_plots_dir)

    with open(cnvloh_json, "r") as cnvloh_json_file:
        cnvloh_dict = json.loads(cnvloh_json_file.read())
        metadata_dict = None
        if "metadata" in cnvloh_dict:
            metadata_dict = cnvloh_dict["metadata"]
        print(metadata_dict)
        if "segments" in cnvloh_dict:
            segments_dict = cnvloh_dict["segments"]
        elif "rows" in cnvloh_dict:
            segments_dict = cnvloh_dict["rows"]
        else:
            segments_dict = cnvloh_dict
        segments_dict, metadata_dict = loh_rescue(segments_dict, metadata_dict)

    if cnvloh_json_2:
        with open(cnvloh_json_2, "r") as cnvloh_json_file_2:
            cnvloh_dict_2 = json.loads(cnvloh_json_file_2.read())
            segments_dict_2 = cnvloh_dict_2["segments"]
            metadata_dict_2 = None
            if "metadata" in cnvloh_dict_2:
                metadata_dict_2 = cnvloh_dict_2["metadata"]
            segments_dict_2, metadata_dict_2 = loh_rescue(
                cnvloh_dict_2, metadata_dict_2
            )
    else:
        segments_dict_2 = None
        metadata_dict_2 = None

    # print(f'Found {len(cnvloh_dict)} total segments')

    if args.focals_json:
        with open(args.focals_json, "r") as focals_json:
            focals_dict = json.loads(focals_json.read())["rows"]
    else:
        focals_dict = None

    # If specified, plots just one; otherwise, plots all
    chromosome = args.chromosome
    if chromosome:
        if chromosome.lower() == "genome":
            plot_genome_html_bokeh(
                out_path=interactive_plots_dir,
                genome_info=genome_info,
                cnvloh_dict=segments_dict,
                genes_table=genes_table,
                genes_list=gene_sublist,
                output_type=out_format,
                cnvloh_dict_2=segments_dict_2,
                focals_dict=focals_dict,
                metadata=metadata_dict,
            )
        else:
            plot_chrom_html_bokeh(
                chromosome=chromosome,
                genome_info=genome_info,
                cnvloh_dict=segments_dict,
                out_path=interactive_plots_dir,
                genes_table=genes_table,
                genes_list=gene_sublist,
                output_type=out_format,
                cnvloh_dict_2=segments_dict_2,
                focals_dict=focals_dict,
                metadata=metadata_dict,
            )
    else:
        plot_all_chrom_html(
            genome_info=genome_info,
            cnvloh_dict=segments_dict,
            out_path=interactive_plots_dir,
            genes_list=gene_sublist,
            genes_table=genes_table,
            output_type=out_format,
            cnvloh_dict_2=segments_dict_2,
            focals_dict=focals_dict,
            metadata=metadata_dict,
        )

    time_elapsed = time.time() - time_start
    print(f"Total runtime: {time_elapsed:.2f} seconds.")
