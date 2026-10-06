#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script: rt_qpcr_analysis.py

Description:
    Analyse RT-qPCR data using the Delta-Delta Ct method and generate
    publication-style figures for basal thionin expression and
    aphid-responsive expression.

    Technical replicates are averaged within each biological replicate
    before Delta Ct calculation. Statistical analyses are performed
    using biological-replicate Delta Ct values.

Analysis:
    1. Basal expression analysis
       - Untreated plants only.
       - Akashinriki untreated plants are used as the calibrator.
       - Four genotypes are compared using one-way ANOVA followed by
         Tukey's HSD test.

    2. Aphid-response analysis
       - Each genotype is analysed independently.
       - The corresponding clip-cage control is used as the calibrator.
       - Clip-cage control, M. persicae infestation, and R. padi
         infestation are compared using one-way ANOVA followed by
         Tukey's HSD test.
       - Only the two aphid treatments are displayed as bars.

Input:
    RT-qPCR Ct data in Excel format.

    Required columns:
        genotype
        CT
        reference_gene
        repeat
        treatment
        type

    The 'type' column should contain:
        sample
        reference_gene

Output:
    1. rt-qPCR_analysis_results.xlsx

       Includes:
           filtered raw data
           technical replicate QC
           mean sample Ct
           mean reference Ct
           biological-replicate Delta Ct
           calibrator values
           Delta-Delta Ct
           relative expression
           log2-transformed relative expression
           summary statistics
           ANOVA results
           Tukey HSD results

    2. figures/

       Basal-expression figures:
           Actin_Figure1_Untreated_ANOVA_Tukey.tif
           Ubiquitin_Figure1_Untreated_ANOVA_Tukey.tif

       Aphid-response figures:
           Actin_Figure2_Aphid_ANOVA_Tukey.tif
           Ubiquitin_Figure2_Aphid_ANOVA_Tukey.tif
           Figure2_Combined_AB.tif

Calculation:
    Delta Ct = Ct_target - Ct_reference

    Basal expression:
        DeltaDelta Ct =
            Delta Ct_sample -
            mean Delta Ct_Akashinriki untreated

    Aphid-response expression:
        DeltaDelta Ct =
            Delta Ct_sample -
            mean Delta Ct_genotype-specific clip-cage control

    Relative expression = 2^-DeltaDeltaCt

    log2 relative expression =
        log2(2^-DeltaDeltaCt)
        = -DeltaDeltaCt

Statistics:
    Statistics are performed on biological-replicate Delta Ct values.

    Basal expression:
        one-way ANOVA + Tukey HSD across genotypes

    Aphid response:
        one-way ANOVA + Tukey HSD across treatments
        within each genotype

Significance:
    **  P < 0.01
    *   P < 0.05
    ns  P >= 0.05

Usage:
    python rt_qpcr_analysis.py \
        --input_file input.xlsx \
        --output_dir output
"""


# =====================================================================
# Imports
# =====================================================================

import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.lines import Line2D
from scipy.stats import f_oneway, studentized_range


warnings.filterwarnings("ignore")


# =====================================================================
# Global plotting parameters
# =====================================================================

plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 14
plt.rcParams["axes.linewidth"] = 1.3


# =====================================================================
# Analysis parameters
# =====================================================================

GENOTYPE_ORDER = [
    "Akashinriki",
    "HOR10350",
    "HOR21599",
    "Morex",
]

REFERENCE_GENES = [
    "Actin",
    "Ubiquitin",
]


UNTREATED_TREATMENT = "Untreated control"
CLIP_CONTROL_TREATMENT = "Clip-cage control"
MP_TREATMENT = "M.persicae infestation"
RP_TREATMENT = "R.padi infestation"


ALL_TREATMENTS = [
    UNTREATED_TREATMENT,
    CLIP_CONTROL_TREATMENT,
    MP_TREATMENT,
    RP_TREATMENT,
]


# Clip-cage control is required for DeltaDeltaCt calculation
# and statistical analysis, but is not displayed as an aphid-response bar.

TREATMENTS_FOR_APHID_CALCULATION = [
    CLIP_CONTROL_TREATMENT,
    MP_TREATMENT,
    RP_TREATMENT,
]

TREATMENTS_FOR_APHID_PLOT = [
    MP_TREATMENT,
    RP_TREATMENT,
]


TREATMENT_DISPLAY = {
    CLIP_CONTROL_TREATMENT: "Clip-cage control",
    MP_TREATMENT: r"$M.\ persicae$",
    RP_TREATMENT: r"$R.\ padi$",
}


# Final publication colors

FIG1_BAR_COLOR = "#D9D9D9"

TREATMENT_COLORS = {
    MP_TREATMENT: "#EED1CC",
    RP_TREATMENT: "#C9DCC4",
}


# Biological replicate symbols

REPEAT_ORDER = [1, 2, 3]

REPEAT_MARKERS = {
    1: "o",
    2: "s",
    3: "^",
}

REPEAT_LABELS = {
    1: "Repeat 1",
    2: "Repeat 2",
    3: "Repeat 3",
}


# Significance thresholds

P_STAR = 0.05
P_DOUBLE_STAR = 0.01


# Figure fonts

TICK_FONT_SIZE = 14
AXIS_FONT_SIZE = 17
SIGNIFICANCE_FONT_SIZE = 16
LEGEND_FONT_SIZE = 14
PANEL_LABEL_FONT_SIZE = 24


# =====================================================================
# Argument parser
# =====================================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Analyse RT-qPCR data using the Delta-Delta Ct method, "
            "perform ANOVA/Tukey HSD tests on biological-replicate "
            "Delta Ct values, and generate publication-style figures."
        )
    )

    parser.add_argument(
        "--input_file",
        required=True,
        help="RT-qPCR Ct data in Excel format.",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Output directory.",
    )

    return parser.parse_args()


# =====================================================================
# Utility functions
# =====================================================================

def p_to_star(p):
    """
    Convert a P value to a significance label.
    """

    if pd.isna(p):
        return ""

    if p < P_DOUBLE_STAR:
        return "**"

    if p < P_STAR:
        return "*"

    return "ns"


def add_black_border(ax):
    """
    Add a black frame and publication-style ticks to an axis.
    """

    for spine_name in [
        "left",
        "right",
        "top",
        "bottom",
    ]:

        ax.spines[spine_name].set_visible(True)
        ax.spines[spine_name].set_color("black")
        ax.spines[spine_name].set_linewidth(1.3)

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        width=1.2,
        length=5,
        color="black",
        labelcolor="black",
        labelsize=TICK_FONT_SIZE,
    )


def get_repeat_legend_handles():
    """
    Generate legend handles for biological replicates.
    """

    handles = []

    for repeat_id in REPEAT_ORDER:

        handles.append(
            Line2D(
                [0],
                [0],
                marker=REPEAT_MARKERS[repeat_id],
                linestyle="None",
                markerfacecolor="black",
                markeredgecolor="black",
                markeredgewidth=0.7,
                markersize=7,
                label=REPEAT_LABELS[repeat_id],
            )
        )

    return handles


# =====================================================================
# Data loading and cleaning
# =====================================================================

def load_input_data(input_file):
    """
    Read the RT-qPCR Ct dataset.
    """

    print("Loading input data...")

    df = pd.read_excel(
        input_file,
        na_values=[
            "na",
            "NA",
            "Na",
            "n/a",
            "N/A",
            "",
        ],
    )

    print(
        f"Loaded {df.shape[0]} rows "
        f"and {df.shape[1]} columns."
    )

    return df


def clean_input_data(df):
    """
    Validate and clean the input dataset.
    """

    print("Cleaning input data...")

    required_columns = [
        "genotype",
        "CT",
        "reference_gene",
        "repeat",
        "treatment",
        "type",
    ]

    missing_columns = [
        col
        for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:

        raise ValueError(
            "Input file is missing required columns: "
            + ", ".join(missing_columns)
        )


    # -------------------------------------------------------------
    # Clean text fields
    # -------------------------------------------------------------

    for col in [
        "genotype",
        "reference_gene",
        "treatment",
        "type",
    ]:

        df[col] = (
            df[col]
            .astype(str)
            .str.strip()
        )


    # -------------------------------------------------------------
    # Numeric fields
    # -------------------------------------------------------------

    df["CT"] = pd.to_numeric(
        df["CT"],
        errors="coerce",
    )

    df["repeat"] = pd.to_numeric(
        df["repeat"],
        errors="coerce",
    )


    # -------------------------------------------------------------
    # Standardize genotype names
    # -------------------------------------------------------------

    df["genotype"] = df["genotype"].replace({
        "Akainriki": "Akashinriki",
        "Akaishinriki": "Akashinriki",
        "HOR1599": "HOR21599",
    })


    # -------------------------------------------------------------
    # Standardize reference gene names
    # -------------------------------------------------------------

    df["reference_gene"] = df["reference_gene"].replace({
        "actin": "Actin",
        "ACTIN": "Actin",
        "ubiquitin": "Ubiquitin",
        "UBIQUITIN": "Ubiquitin",
    })


    # -------------------------------------------------------------
    # Standardize treatment names
    # -------------------------------------------------------------

    df["treatment"] = df["treatment"].replace({
        "No treatment": UNTREATED_TREATMENT,
        "Empty cage": CLIP_CONTROL_TREATMENT,
        "Cage control": CLIP_CONTROL_TREATMENT,
        "R.padi": RP_TREATMENT,
        "M.persicae": MP_TREATMENT,
    })


    # -------------------------------------------------------------
    # Standardize data type
    # -------------------------------------------------------------

    df["type"] = df["type"].replace({
        "reference gene": "reference_gene",
        "reference_gene": "reference_gene",
        "sample": "sample",
    })


    # -------------------------------------------------------------
    # Remove invalid Ct values
    # -------------------------------------------------------------

    before_drop = len(df)

    df = df.dropna(
        subset=["CT"]
    ).copy()

    removed = before_drop - len(df)

    if removed > 0:

        print(
            f"Removed {removed} rows "
            "with missing/non-numeric Ct values."
        )


    # -------------------------------------------------------------
    # Keep only expected experimental groups
    # -------------------------------------------------------------

    df = df[
        df["genotype"].isin(
            GENOTYPE_ORDER
        )
    ].copy()

    df = df[
        df["reference_gene"].isin(
            REFERENCE_GENES
        )
    ].copy()

    df = df[
        df["treatment"].isin(
            ALL_TREATMENTS
        )
    ].copy()


    print(
        "Detected genotypes:",
        list(df["genotype"].unique()),
    )

    print(
        "Detected treatments:",
        list(df["treatment"].unique()),
    )

    print(
        "Detected reference genes:",
        list(df["reference_gene"].unique()),
    )

    print(
        "Detected data types:",
        list(df["type"].unique()),
    )


    return df


# =====================================================================
# Technical replicate processing
# =====================================================================

def calculate_per_repeat_delta_ct(df):
    """
    Average technical replicates first and calculate one Delta Ct value
    for each biological replicate.

    Delta Ct = mean target Ct - mean reference Ct
    """

    print(
        "Averaging technical replicates..."
    )


    # -------------------------------------------------------------
    # Technical replicate QC
    # -------------------------------------------------------------

    technical_qc = (
        df
        .groupby(
            [
                "genotype",
                "reference_gene",
                "repeat",
                "treatment",
                "type",
            ]
        )
        .agg(
            mean_CT=("CT", "mean"),
            n_valid_technical=("CT", "count"),
            SD_technical=("CT", "std"),
        )
        .reset_index()
    )


    # -------------------------------------------------------------
    # Target gene Ct
    # -------------------------------------------------------------

    sample_mean = (
        df[
            df["type"] == "sample"
        ]
        .groupby(
            [
                "genotype",
                "reference_gene",
                "repeat",
                "treatment",
            ]
        )
        .agg(
            sample_CT=("CT", "mean"),
            sample_n=("CT", "count"),
            sample_SD=("CT", "std"),
        )
        .reset_index()
    )


    # -------------------------------------------------------------
    # Reference gene Ct
    # -------------------------------------------------------------

    reference_mean = (
        df[
            df["type"] == "reference_gene"
        ]
        .groupby(
            [
                "genotype",
                "reference_gene",
                "repeat",
                "treatment",
            ]
        )
        .agg(
            reference_CT=("CT", "mean"),
            reference_n=("CT", "count"),
            reference_SD=("CT", "std"),
        )
        .reset_index()
    )


    # -------------------------------------------------------------
    # Merge target/reference values
    # -------------------------------------------------------------

    per_repeat = pd.merge(
        sample_mean,
        reference_mean,
        on=[
            "genotype",
            "reference_gene",
            "repeat",
            "treatment",
        ],
        how="outer",
    )


    per_repeat["DeltaCt"] = (
        per_repeat["sample_CT"]
        - per_repeat["reference_CT"]
    )


    missing_delta = per_repeat[
        per_repeat["DeltaCt"].isna()
    ].copy()


    if not missing_delta.empty:

        print(
            "\nWARNING: "
            "Some biological replicates "
            "have missing DeltaCt."
        )

        print(
            missing_delta[
                [
                    "genotype",
                    "reference_gene",
                    "repeat",
                    "treatment",
                    "sample_CT",
                    "reference_CT",
                ]
            ]
        )


    per_repeat_valid = per_repeat[
        per_repeat["DeltaCt"].notna()
    ].copy()


    print(
        "Valid biological-replicate records:",
        len(per_repeat_valid),
    )


    return (
        technical_qc,
        sample_mean,
        reference_mean,
        per_repeat,
        per_repeat_valid,
        missing_delta,
    )


# =====================================================================
# Delta-Delta Ct analysis
# =====================================================================

def calculate_untreated_ddct(per_repeat_valid):
    """
    Basal-expression analysis.

    Only untreated samples are analysed.

    Akashinriki untreated plants are used as the calibrator.
    """

    data = per_repeat_valid[
        per_repeat_valid["treatment"]
        == UNTREATED_TREATMENT
    ].copy()


    calibrator = (
        data[
            data["genotype"]
            == "Akashinriki"
        ]
        .groupby(
            "reference_gene"
        )
        .agg(
            calibrator_DeltaCt=("DeltaCt", "mean"),
            calibrator_SD=("DeltaCt", "std"),
            calibrator_n=("DeltaCt", "count"),
        )
        .reset_index()
    )


    missing_refs = (
        set(REFERENCE_GENES)
        - set(calibrator["reference_gene"])
    )


    if missing_refs:

        raise ValueError(
            "Missing untreated Akashinriki "
            f"calibrator for: {missing_refs}"
        )


    data = pd.merge(
        data,
        calibrator[
            [
                "reference_gene",
                "calibrator_DeltaCt",
            ]
        ],
        on="reference_gene",
        how="left",
    )


    data["DeltaDeltaCt"] = (
        data["DeltaCt"]
        - data["calibrator_DeltaCt"]
    )


    data["Relative_expression"] = (
        2 ** (
            -data["DeltaDeltaCt"]
        )
    )


    data["log2FC"] = (
        -data["DeltaDeltaCt"]
    )


    return (
        data,
        calibrator,
    )


def calculate_treatment_ddct(per_repeat_valid):
    """
    Aphid-response analysis.

    Each genotype is analysed separately using its corresponding
    clip-cage control as the calibrator.
    """

    data = per_repeat_valid[
        per_repeat_valid["treatment"].isin(
            TREATMENTS_FOR_APHID_CALCULATION
        )
    ].copy()


    calibrator = (
        data[
            data["treatment"]
            == CLIP_CONTROL_TREATMENT
        ]
        .groupby(
            [
                "genotype",
                "reference_gene",
            ]
        )
        .agg(
            calibrator_DeltaCt=("DeltaCt", "mean"),
            calibrator_SD=("DeltaCt", "std"),
            calibrator_n=("DeltaCt", "count"),
        )
        .reset_index()
    )


    expected = pd.MultiIndex.from_product(
        [
            GENOTYPE_ORDER,
            REFERENCE_GENES,
        ],
        names=[
            "genotype",
            "reference_gene",
        ],
    )


    observed = pd.MultiIndex.from_frame(
        calibrator[
            [
                "genotype",
                "reference_gene",
            ]
        ]
    )


    missing_calibrators = (
        expected.difference(observed)
    )


    if len(missing_calibrators) > 0:

        raise ValueError(
            "Missing genotype-specific "
            "clip-cage calibrators:\n"
            + str(
                list(missing_calibrators)
            )
        )


    data = pd.merge(
        data,
        calibrator[
            [
                "genotype",
                "reference_gene",
                "calibrator_DeltaCt",
            ]
        ],
        on=[
            "genotype",
            "reference_gene",
        ],
        how="left",
    )


    data["DeltaDeltaCt"] = (
        data["DeltaCt"]
        - data["calibrator_DeltaCt"]
    )


    data["Relative_expression"] = (
        2 ** (
            -data["DeltaDeltaCt"]
        )
    )


    data["log2FC"] = (
        -data["DeltaDeltaCt"]
    )


    return (
        data,
        calibrator,
    )


# =====================================================================
# Summary statistics
# =====================================================================

def make_summary(input_df):
    """
    Calculate biological-replicate summary statistics.
    """

    summary = (
        input_df
        .groupby(
            [
                "genotype",
                "reference_gene",
                "treatment",
            ]
        )
        .agg(
            mean_log2FC=("log2FC", "mean"),
            median_log2FC=("log2FC", "median"),
            SD_log2FC=("log2FC", "std"),
            SEM_log2FC=("log2FC", "sem"),
            mean_relative_expression=(
                "Relative_expression",
                "mean",
            ),
            SEM_relative_expression=(
                "Relative_expression",
                "sem",
            ),
            n=("log2FC", "count"),
        )
        .reset_index()
    )

    return summary


# =====================================================================
# ANOVA + Tukey HSD helper
# =====================================================================

def calculate_anova_tukey(
    groups,
    group_order,
    group_column_name,
    reference_gene,
    extra_columns=None,
):
    """
    Perform one-way ANOVA followed by Tukey HSD.

    The implementation supports unequal group sizes using the
    Tukey-Kramer standard error.

    Returns:
        ANOVA table
        Tukey HSD table
    """

    if extra_columns is None:
        extra_columns = {}


    group_names = [
        group_name
        for group_name in group_order
        if (
            group_name in groups
            and len(groups[group_name]) > 0
        )
    ]


    if len(group_names) < 2:

        return (
            pd.DataFrame(),
            pd.DataFrame(),
        )


    group_arrays = [
        np.asarray(
            groups[group_name],
            dtype=float,
        )
        for group_name
        in group_names
    ]


    # -------------------------------------------------------------
    # One-way ANOVA
    # -------------------------------------------------------------

    F_statistic, anova_p = (
        f_oneway(
            *group_arrays
        )
    )


    all_values = np.concatenate(
        group_arrays
    )

    grand_mean = np.mean(
        all_values
    )

    k = len(
        group_arrays
    )

    N = len(
        all_values
    )


    ss_between = sum(
        len(values)
        * (
            np.mean(values)
            - grand_mean
        ) ** 2
        for values
        in group_arrays
    )


    ss_within = sum(
        np.sum(
            (
                values
                - np.mean(values)
            ) ** 2
        )
        for values
        in group_arrays
    )


    df_between = k - 1
    df_within = N - k


    ms_between = (
        ss_between / df_between
        if df_between > 0
        else np.nan
    )


    ms_within = (
        ss_within / df_within
        if df_within > 0
        else np.nan
    )


    anova_record = {
        "reference_gene": reference_gene,
        **extra_columns,
        "source": group_column_name,
        "sum_sq": ss_between,
        "df": df_between,
        "mean_sq": ms_between,
        "F_statistic": F_statistic,
        "p_value": anova_p,
    }


    residual_record = {
        "reference_gene": reference_gene,
        **extra_columns,
        "source": "residual",
        "sum_sq": ss_within,
        "df": df_within,
        "mean_sq": ms_within,
        "F_statistic": np.nan,
        "p_value": np.nan,
    }


    anova_table = pd.DataFrame(
        [
            anova_record,
            residual_record,
        ]
    )


    if (
        df_within <= 0
        or pd.isna(ms_within)
    ):

        return (
            anova_table,
            pd.DataFrame(),
        )


    # -------------------------------------------------------------
    # Tukey HSD / Tukey-Kramer
    # -------------------------------------------------------------

    q_critical = (
        studentized_range.ppf(
            0.95,
            k,
            df_within,
        )
    )


    tukey_results = []


    for i in range(
        len(group_names)
    ):

        for j in range(
            i + 1,
            len(group_names)
        ):

            group1 = group_names[i]
            group2 = group_names[j]

            values1 = groups[group1]
            values2 = groups[group2]

            n1 = len(values1)
            n2 = len(values2)

            mean1 = np.mean(values1)
            mean2 = np.mean(values2)

            mean_diff = (
                mean2 - mean1
            )


            tukey_se = np.sqrt(
                ms_within
                / 2.0
                * (
                    1.0 / n1
                    + 1.0 / n2
                )
            )


            if (
                tukey_se == 0
                or pd.isna(tukey_se)
            ):

                q_statistic = np.nan
                p_adjusted = np.nan
                ci_lower = np.nan
                ci_upper = np.nan
                reject_h0 = False

            else:

                q_statistic = (
                    abs(mean_diff)
                    / tukey_se
                )

                p_adjusted = (
                    studentized_range.sf(
                        q_statistic,
                        k,
                        df_within,
                    )
                )

                half_width = (
                    q_critical
                    * tukey_se
                )

                ci_lower = (
                    mean_diff
                    - half_width
                )

                ci_upper = (
                    mean_diff
                    + half_width
                )

                reject_h0 = (
                    p_adjusted < 0.05
                )


            result = {
                "reference_gene": reference_gene,
                **extra_columns,
                f"{group_column_name}1": group1,
                f"{group_column_name}2": group2,
                "comparison":
                    f"{group1} vs {group2}",
                f"n_{group_column_name}1": n1,
                f"n_{group_column_name}2": n2,
                "mean_DeltaCt_difference":
                    mean_diff,
                "q_statistic":
                    q_statistic,
                "CI_lower":
                    ci_lower,
                "CI_upper":
                    ci_upper,
                "p_value_adjusted":
                    p_adjusted,
                "reject_H0":
                    reject_h0,
                "significance":
                    p_to_star(
                        p_adjusted
                    ),
            }


            tukey_results.append(
                result
            )


    tukey_table = pd.DataFrame(
        tukey_results
    )


    return (
        anova_table,
        tukey_table,
    )


# =====================================================================
# Figure 1 statistics
# =====================================================================

def significance_figure1(
    per_repeat_valid,
    ref_gene,
):
    """
    Compare untreated Delta Ct values among genotypes.
    """

    data = per_repeat_valid[
        (
            per_repeat_valid[
                "reference_gene"
            ]
            == ref_gene
        )
        &
        (
            per_repeat_valid[
                "treatment"
            ]
            == UNTREATED_TREATMENT
        )
    ].copy()


    groups = {}


    for genotype in GENOTYPE_ORDER:

        values = (
            data[
                data["genotype"]
                == genotype
            ]["DeltaCt"]
            .dropna()
            .astype(float)
            .values
        )

        if len(values) > 0:
            groups[genotype] = values


    return calculate_anova_tukey(
        groups=groups,
        group_order=GENOTYPE_ORDER,
        group_column_name="genotype",
        reference_gene=ref_gene,
    )


# =====================================================================
# Figure 2 statistics
# =====================================================================

def significance_figure2(
    per_repeat_valid,
    ref_gene,
):
    """
    Within each genotype, compare Delta Ct values among:

        Clip-cage control
        M. persicae infestation
        R. padi infestation
    """

    data = per_repeat_valid[
        per_repeat_valid[
            "reference_gene"
        ]
        == ref_gene
    ].copy()


    treatment_order = [
        CLIP_CONTROL_TREATMENT,
        MP_TREATMENT,
        RP_TREATMENT,
    ]


    anova_tables = []
    tukey_tables = []


    for genotype in GENOTYPE_ORDER:

        genotype_data = data[
            (
                data["genotype"]
                == genotype
            )
            &
            (
                data["treatment"]
                .isin(
                    treatment_order
                )
            )
        ].copy()


        groups = {}


        for treatment in treatment_order:

            values = (
                genotype_data[
                    genotype_data[
                        "treatment"
                    ]
                    == treatment
                ]["DeltaCt"]
                .dropna()
                .astype(float)
                .values
            )

            if len(values) > 0:
                groups[treatment] = values


        anova_df, tukey_df = (
            calculate_anova_tukey(
                groups=groups,
                group_order=treatment_order,
                group_column_name="treatment",
                reference_gene=ref_gene,
                extra_columns={
                    "genotype":
                        genotype
                },
            )
        )


        if not anova_df.empty:
            anova_tables.append(
                anova_df
            )


        if not tukey_df.empty:
            tukey_tables.append(
                tukey_df
            )


    anova_all = (
        pd.concat(
            anova_tables,
            ignore_index=True,
        )
        if anova_tables
        else pd.DataFrame()
    )


    tukey_all = (
        pd.concat(
            tukey_tables,
            ignore_index=True,
        )
        if tukey_tables
        else pd.DataFrame()
    )


    return (
        anova_all,
        tukey_all,
    )


# =====================================================================
# Figure 1: basal expression
# =====================================================================

def plot_untreated(
    raw_data,
    summary_data,
    ref_gene,
    significance_df,
    figure_dir,
):
    """
    Plot basal expression among untreated genotypes.
    """

    raw = raw_data[
        raw_data["reference_gene"]
        == ref_gene
    ].copy()


    summary = summary_data[
        summary_data["reference_gene"]
        == ref_gene
    ].copy()


    # -------------------------------------------------------------
    # Plot genotypes in ascending mean expression order
    # -------------------------------------------------------------

    plot_genotype_order = (
        summary[
            summary["genotype"].isin(
                GENOTYPE_ORDER
            )
        ]
        .set_index("genotype")
        .reindex(
            GENOTYPE_ORDER
        )
        .sort_values(
            "mean_log2FC",
            ascending=True,
            na_position="last",
        )
        .index
        .tolist()
    )


    print(
        f"Figure 1 plotting order "
        f"({ref_gene}, low to high):",
        plot_genotype_order,
    )


    fig, ax = plt.subplots(
        figsize=(8, 6)
    )


    x = np.arange(
        len(plot_genotype_order)
    )


    means = []
    sems = []


    for genotype in plot_genotype_order:

        temp = summary[
            summary["genotype"]
            == genotype
        ]

        if temp.empty:

            means.append(
                np.nan
            )

            sems.append(
                np.nan
            )

        else:

            means.append(
                temp[
                    "mean_log2FC"
                ].iloc[0]
            )

            sems.append(
                temp[
                    "SEM_log2FC"
                ].iloc[0]
            )


    bars = ax.bar(
        x,
        means,
        yerr=sems,
        width=0.65,
        capsize=4,
        edgecolor="black",
        linewidth=1.2,
    )


    for bar in bars:

        bar.set_facecolor(
            FIG1_BAR_COLOR
        )


    # -------------------------------------------------------------
    # Biological replicate points
    # -------------------------------------------------------------

    repeat_x_offsets = {
        1: -0.07,
        2: 0.00,
        3: 0.07,
    }


    for i, genotype in enumerate(
        plot_genotype_order
    ):

        genotype_data = raw[
            raw["genotype"]
            == genotype
        ]


        for repeat_id in REPEAT_ORDER:

            temp = genotype_data[
                genotype_data["repeat"]
                == repeat_id
            ]


            if temp.empty:
                continue


            value = temp[
                "log2FC"
            ].iloc[0]


            if pd.isna(value):
                continue


            ax.scatter(
                x[i]
                + repeat_x_offsets[
                    repeat_id
                ],
                value,
                marker=REPEAT_MARKERS[
                    repeat_id
                ],
                s=48,
                facecolor="black",
                edgecolor="black",
                linewidth=0.7,
                zorder=10,
            )


    # -------------------------------------------------------------
    # Determine annotation range
    # -------------------------------------------------------------

    upper_values = []
    lower_values = []


    for mean_value, sem_value in zip(
        means,
        sems,
    ):

        if pd.isna(mean_value):
            continue


        sem_value = (
            0
            if pd.isna(sem_value)
            else sem_value
        )


        upper_values.append(
            mean_value + sem_value
        )

        lower_values.append(
            mean_value - sem_value
        )


    raw_values = (
        raw["log2FC"]
        .dropna()
        .values
    )


    if len(raw_values) > 0:

        upper_values.extend(
            raw_values.tolist()
        )

        lower_values.extend(
            raw_values.tolist()
        )


    if upper_values:

        base_top = max(
            upper_values
        )

        base_bottom = min(
            lower_values
        )

    else:

        base_top = 1
        base_bottom = -1


    data_range = (
        base_top
        - base_bottom
    )


    if data_range == 0:
        data_range = 1


    # -------------------------------------------------------------
    # Six Tukey pairwise comparisons
    # -------------------------------------------------------------

    genotype_to_x = {
        genotype: i
        for i, genotype
        in enumerate(
            plot_genotype_order
        )
    }


    pair_lookup = {}


    for _, row in (
        significance_df.iterrows()
    ):

        key = frozenset(
            [
                row["genotype1"],
                row["genotype2"],
            ]
        )

        pair_lookup[key] = row


    if len(
        plot_genotype_order
    ) != 4:

        raise ValueError(
            "Figure 1 bracket layout "
            "requires exactly four genotypes."
        )


    g0, g1, g2, g3 = (
        plot_genotype_order
    )


    comparison_levels = [
        [
            (g0, g1),
            (g1, g2),
            (g2, g3),
        ],
        [
            (g0, g2),
        ],
        [
            (g1, g3),
        ],
        [
            (g0, g3),
        ],
    ]


    first_y = (
        base_top
        + data_range * 0.10
    )

    level_spacing = (
        data_range * 0.115
    )

    bracket_height = (
        data_range * 0.025
    )

    text_gap = (
        data_range * 0.025
    )

    highest_annotation = (
        base_top
    )


    for level, pairs in enumerate(
        comparison_levels
    ):

        y = (
            first_y
            + level
            * level_spacing
        )


        for g_left, g_right in pairs:

            row = pair_lookup.get(
                frozenset(
                    [
                        g_left,
                        g_right,
                    ]
                )
            )


            if row is None:
                continue


            label = (
                row["significance"]
                if (
                    pd.notna(
                        row["significance"]
                    )
                    and row[
                        "significance"
                    ] != ""
                )
                else "ns"
            )


            x1 = genotype_to_x[
                g_left
            ]

            x2 = genotype_to_x[
                g_right
            ]


            if x1 > x2:
                x1, x2 = x2, x1


            ax.plot(
                [
                    x1,
                    x1,
                    x2,
                    x2,
                ],
                [
                    y,
                    y + bracket_height,
                    y + bracket_height,
                    y,
                ],
                color="black",
                linewidth=1.05,
                clip_on=False,
            )


            text_y = (
                y
                + bracket_height
                + text_gap
            )


            ax.text(
                (x1 + x2) / 2,
                text_y,
                label,
                ha="center",
                va="bottom",
                fontsize=SIGNIFICANCE_FONT_SIZE,
                fontweight="bold",
            )


            highest_annotation = max(
                highest_annotation,
                text_y,
            )


    # -------------------------------------------------------------
    # Axis formatting
    # -------------------------------------------------------------

    ax.axhline(
        y=0,
        color="black",
        linewidth=0.8,
    )


    ax.set_xticks(
        x
    )


    ax.set_xticklabels(
        plot_genotype_order,
        fontsize=TICK_FONT_SIZE,
    )


    ax.set_xlabel(
        "Genotype",
        fontsize=AXIS_FONT_SIZE,
    )


    ax.set_ylabel(
        r"$\log_2(2^{-\Delta\Delta C_t})$",
        fontsize=AXIS_FONT_SIZE,
    )


    add_black_border(
        ax
    )


    ax.grid(
        False
    )


    ymin, ymax = (
        ax.get_ylim()
    )


    y_range = (
        ymax - ymin
    )


    if y_range == 0:
        y_range = 1


    final_top = max(
        ymax,
        highest_annotation
        + data_range * 0.08,
    )


    ax.set_ylim(
        ymin - y_range * 0.05,
        final_top,
    )


    fig.subplots_adjust(
        right=0.97,
        top=0.95,
        bottom=0.18,
        left=0.16,
    )


    output_name = (
        f"{ref_gene}_"
        "Figure1_Untreated_ANOVA_Tukey"
    )


    fig.savefig(
        figure_dir
        / f"{output_name}.tif",
        dpi=300,
    )


    plt.close(
        fig
    )


    print(
        f"Saved Figure 1: "
        f"{ref_gene}"
    )


# =====================================================================
# Figure 2 panel
# =====================================================================

def draw_treatments_panel(
    ax,
    raw_data,
    summary_data,
    ref_gene,
    significance_df,
):
    """
    Draw one aphid-response panel.
    """

    raw = raw_data[
        raw_data["reference_gene"]
        == ref_gene
    ].copy()


    summary = summary_data[
        summary_data["reference_gene"]
        == ref_gene
    ].copy()


    raw_plot = raw[
        raw["treatment"].isin(
            TREATMENTS_FOR_APHID_PLOT
        )
    ].copy()


    summary_plot = summary[
        summary["treatment"].isin(
            TREATMENTS_FOR_APHID_PLOT
        )
    ].copy()


    # -------------------------------------------------------------
    # Plot genotypes from low to high mean aphid response
    # -------------------------------------------------------------

    genotype_plot_mean = (
        summary_plot
        .groupby(
            "genotype"
        )["mean_log2FC"]
        .mean()
        .reindex(
            GENOTYPE_ORDER
        )
    )


    plot_genotype_order = (
        genotype_plot_mean
        .sort_values(
            ascending=True,
            na_position="last",
        )
        .index
        .tolist()
    )


    print(
        f"Figure 2 plotting order "
        f"({ref_gene}, low to high):",
        plot_genotype_order,
    )


    x = np.arange(
        len(plot_genotype_order)
    )


    n_treatments = len(
        TREATMENTS_FOR_APHID_PLOT
    )


    bar_width = 0.30


    offsets = (
        np.arange(
            n_treatments
        )
        - (
            n_treatments - 1
        ) / 2
    ) * bar_width


    bar_information = {}


    repeat_x_offsets = {
        1: -0.035,
        2: 0.000,
        3: 0.035,
    }


    # -------------------------------------------------------------
    # Bars
    # -------------------------------------------------------------

    for j, treatment in enumerate(
        TREATMENTS_FOR_APHID_PLOT
    ):

        means = []
        sems = []


        for genotype in (
            plot_genotype_order
        ):

            temp = summary_plot[
                (
                    summary_plot["genotype"]
                    == genotype
                )
                &
                (
                    summary_plot["treatment"]
                    == treatment
                )
            ]


            if temp.empty:

                means.append(
                    np.nan
                )

                sems.append(
                    np.nan
                )

            else:

                means.append(
                    temp[
                        "mean_log2FC"
                    ].iloc[0]
                )

                sems.append(
                    temp[
                        "SEM_log2FC"
                    ].iloc[0]
                )


        xpos = (
            x + offsets[j]
        )


        ax.bar(
            xpos,
            means,
            yerr=sems,
            width=bar_width,
            capsize=3,
            color=TREATMENT_COLORS[
                treatment
            ],
            edgecolor="black",
            linewidth=1.1,
            label=TREATMENT_DISPLAY[
                treatment
            ],
        )


        # ---------------------------------------------------------
        # Biological replicate points
        # ---------------------------------------------------------

        for i, genotype in enumerate(
            plot_genotype_order
        ):

            genotype_treatment_data = (
                raw_plot[
                    (
                        raw_plot["genotype"]
                        == genotype
                    )
                    &
                    (
                        raw_plot["treatment"]
                        == treatment
                    )
                ]
                .copy()
            )


            raw_group_values = (
                genotype_treatment_data[
                    "log2FC"
                ]
                .dropna()
                .astype(float)
                .values
            )


            mean_value = means[i]

            sem_value = (
                0
                if pd.isna(
                    sems[i]
                )
                else sems[i]
            )


            candidates = []


            if not pd.isna(
                mean_value
            ):

                candidates.append(
                    mean_value
                    + sem_value
                )


            if len(
                raw_group_values
            ) > 0:

                candidates.extend(
                    raw_group_values
                    .tolist()
                )


            local_top = (
                max(candidates)
                if candidates
                else np.nan
            )


            bar_information[
                (
                    genotype,
                    treatment,
                )
            ] = {
                "x": xpos[i],
                "mean": mean_value,
                "sem": sems[i],
                "raw_values":
                    raw_group_values,
                "local_top":
                    local_top,
            }


            for repeat_id in (
                REPEAT_ORDER
            ):

                temp = (
                    genotype_treatment_data[
                        genotype_treatment_data[
                            "repeat"
                        ]
                        == repeat_id
                    ]
                )


                if temp.empty:
                    continue


                value = temp[
                    "log2FC"
                ].iloc[0]


                if pd.isna(value):
                    continue


                ax.scatter(
                    xpos[i]
                    + repeat_x_offsets[
                        repeat_id
                    ],
                    value,
                    marker=REPEAT_MARKERS[
                        repeat_id
                    ],
                    s=45,
                    facecolor="black",
                    edgecolor="black",
                    linewidth=0.7,
                    zorder=10,
                )


    # -------------------------------------------------------------
    # Determine plotting range
    # -------------------------------------------------------------

    all_values = []


    for info in (
        bar_information.values()
    ):

        mean_value = info[
            "mean"
        ]

        sem_value = (
            0
            if pd.isna(
                info["sem"]
            )
            else info["sem"]
        )


        if not pd.isna(
            mean_value
        ):

            all_values.extend(
                [
                    mean_value
                    - sem_value,
                    mean_value
                    + sem_value,
                ]
            )


        all_values.extend(
            info[
                "raw_values"
            ].tolist()
        )


    if all_values:

        base_top = max(
            all_values
        )

        base_bottom = min(
            all_values
        )

    else:

        base_top = 1
        base_bottom = -1


    data_range = (
        base_top
        - base_bottom
    )


    if data_range == 0:
        data_range = 1


    star_gap = (
        data_range * 0.075
    )

    bracket_gap = (
        data_range * 0.11
    )

    bracket_height = (
        data_range * 0.03
    )

    bracket_text_gap = (
        data_range * 0.035
    )


    treatment_annotation_top = {}

    highest_annotation = (
        base_top
    )


    # -------------------------------------------------------------
    # 1. Aphid treatments vs clip-cage control
    # -------------------------------------------------------------

    control_rows = significance_df[
        (
            significance_df[
                "treatment1"
            ]
            == CLIP_CONTROL_TREATMENT
        )
        |
        (
            significance_df[
                "treatment2"
            ]
            == CLIP_CONTROL_TREATMENT
        )
    ].copy()


    for _, row in (
        control_rows.iterrows()
    ):

        genotype = row[
            "genotype"
        ]


        visible_treatment = (
            row["treatment2"]
            if (
                row["treatment1"]
                == CLIP_CONTROL_TREATMENT
            )
            else row["treatment1"]
        )


        key = (
            genotype,
            visible_treatment,
        )


        if (
            visible_treatment
            not in
            TREATMENTS_FOR_APHID_PLOT
            or key
            not in
            bar_information
        ):
            continue


        info = (
            bar_information[
                key
            ]
        )


        if pd.isna(
            info["local_top"]
        ):
            continue


        label = (
            row["significance"]
            if (
                pd.notna(
                    row["significance"]
                )
                and row[
                    "significance"
                ] != ""
            )
            else "ns"
        )


        y_text = (
            info["local_top"]
            + star_gap
        )


        ax.text(
            info["x"],
            y_text,
            label,
            ha="center",
            va="bottom",
            fontsize=SIGNIFICANCE_FONT_SIZE,
            fontweight="bold",
            zorder=20,
        )


        treatment_annotation_top[
            key
        ] = y_text


        highest_annotation = max(
            highest_annotation,
            y_text,
        )


    # -------------------------------------------------------------
    # 2. M. persicae vs R. padi
    # -------------------------------------------------------------

    aphid_pair_rows = significance_df[
        significance_df.apply(
            lambda row:
            {
                row["treatment1"],
                row["treatment2"],
            }
            ==
            {
                MP_TREATMENT,
                RP_TREATMENT,
            },
            axis=1,
        )
    ].copy()


    for _, row in (
        aphid_pair_rows.iterrows()
    ):

        genotype = row[
            "genotype"
        ]


        key_mp = (
            genotype,
            MP_TREATMENT,
        )

        key_rp = (
            genotype,
            RP_TREATMENT,
        )


        if (
            key_mp
            not in bar_information
            or key_rp
            not in bar_information
        ):
            continue


        info_mp = (
            bar_information[
                key_mp
            ]
        )

        info_rp = (
            bar_information[
                key_rp
            ]
        )


        local_candidates = [
            info_mp[
                "local_top"
            ],
            info_rp[
                "local_top"
            ],
            treatment_annotation_top.get(
                key_mp,
                -np.inf,
            ),
            treatment_annotation_top.get(
                key_rp,
                -np.inf,
            ),
        ]


        local_candidates = [
            value
            for value
            in local_candidates
            if (
                pd.notna(value)
                and np.isfinite(
                    value
                )
            )
        ]


        if not local_candidates:
            continue


        y = (
            max(
                local_candidates
            )
            + bracket_gap
        )


        x1 = info_mp["x"]
        x2 = info_rp["x"]


        ax.plot(
            [
                x1,
                x1,
                x2,
                x2,
            ],
            [
                y,
                y + bracket_height,
                y + bracket_height,
                y,
            ],
            color="black",
            linewidth=1.05,
            clip_on=False,
            zorder=15,
        )


        label = (
            row["significance"]
            if (
                pd.notna(
                    row["significance"]
                )
                and row[
                    "significance"
                ] != ""
            )
            else "ns"
        )


        text_y = (
            y
            + bracket_height
            + bracket_text_gap
        )


        ax.text(
            (x1 + x2) / 2,
            text_y,
            label,
            ha="center",
            va="bottom",
            fontsize=SIGNIFICANCE_FONT_SIZE,
            fontweight="bold",
            zorder=20,
        )


        highest_annotation = max(
            highest_annotation,
            text_y,
        )


    # -------------------------------------------------------------
    # Axes
    # -------------------------------------------------------------

    ax.axhline(
        y=0,
        color="black",
        linewidth=0.8,
    )


    ax.set_xticks(
        x
    )


    ax.set_xticklabels(
        plot_genotype_order,
        fontsize=TICK_FONT_SIZE,
    )


    ax.set_xlabel(
        "Genotype",
        fontsize=AXIS_FONT_SIZE,
    )


    ax.set_ylabel(
        r"$\log_2(2^{-\Delta\Delta C_t})$",
        fontsize=AXIS_FONT_SIZE,
    )


    add_black_border(
        ax
    )


    ax.grid(
        False
    )


    ymin, ymax = (
        ax.get_ylim()
    )


    y_range = (
        ymax - ymin
    )


    if y_range == 0:
        y_range = 1


    final_top = max(
        ymax,
        highest_annotation
        + data_range * 0.10,
    )


    ax.set_ylim(
        ymin
        - y_range * 0.08,
        final_top,
    )


    return (
        ax.get_legend_handles_labels()
    )


# =====================================================================
# Figure 2 individual panels
# =====================================================================

def plot_treatments(
    raw_data,
    summary_data,
    ref_gene,
    significance_df,
    figure_dir,
):
    """
    Save an individual aphid-response panel.
    """

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )


    handles, labels = (
        draw_treatments_panel(
            ax=ax,
            raw_data=raw_data,
            summary_data=summary_data,
            ref_gene=ref_gene,
            significance_df=significance_df,
        )
    )


    fig.legend(
        handles,
        labels,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        loc="upper left",
        bbox_to_anchor=(
            0.80,
            0.90,
        ),
        borderaxespad=0,
    )


    fig.subplots_adjust(
        right=0.76,
        top=0.95,
        bottom=0.18,
        left=0.16,
    )


    output_name = (
        f"{ref_gene}_"
        "Figure2_Aphid_ANOVA_Tukey"
    )


    fig.savefig(
        figure_dir
        / f"{output_name}.tif",
        dpi=300,
    )


    plt.close(
        fig
    )


    print(
        f"Saved Figure 2: "
        f"{ref_gene}"
    )


# =====================================================================
# Combined Figure 2
# =====================================================================

def combine_figure2_panels(
    reference_genes,
    treatment_data,
    treatment_summary,
    fig2_significance,
    figure_dir,
):
    """
    Combine Actin and Ubiquitin aphid-response panels into one figure.
    """

    if len(
        reference_genes
    ) != 2:

        print(
            "Skipping combined Figure 2: "
            "exactly two reference genes are required."
        )

        return


    fig, axes = plt.subplots(
        1,
        2,
        figsize=(16, 6.5),
    )


    for (
        panel_label,
        ax,
        ref_gene,
    ) in zip(
        [
            "A",
            "B",
        ],
        axes,
        reference_genes,
    ):

        draw_treatments_panel(
            ax=ax,
            raw_data=treatment_data,
            summary_data=treatment_summary,
            ref_gene=ref_gene,
            significance_df=(
                fig2_significance[
                    ref_gene
                ]
            ),
        )


        ax.text(
            -0.08,
            1.04,
            panel_label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=PANEL_LABEL_FONT_SIZE,
            fontweight="bold",
            color="black",
            clip_on=False,
            zorder=30,
        )


    # -------------------------------------------------------------
    # Shared legend
    # Order: R. padi, M. persicae
    # -------------------------------------------------------------

    shared_handles = [

        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=TREATMENT_COLORS[
                RP_TREATMENT
            ],
            edgecolor="black",
            linewidth=1.0,
        ),

        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=TREATMENT_COLORS[
                MP_TREATMENT
            ],
            edgecolor="black",
            linewidth=1.0,
        ),

    ]


    shared_labels = [
        TREATMENT_DISPLAY[
            RP_TREATMENT
        ],
        TREATMENT_DISPLAY[
            MP_TREATMENT
        ],
    ]


    fig.legend(
        shared_handles,
        shared_labels,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.025,
        ),
        ncol=2,
        columnspacing=2.0,
        handletextpad=0.6,
    )


    fig.subplots_adjust(
        left=0.07,
        right=0.985,
        top=0.87,
        bottom=0.23,
        wspace=0.22,
    )


    output_name = (
        "Figure2_Combined_AB"
    )


    fig.savefig(
        figure_dir
        / f"{output_name}.tif",
        dpi=300,
    )


    plt.close(
        fig
    )


    print(
        "Saved combined Figure 2: A/B"
    )


# =====================================================================
# Excel output
# =====================================================================

def save_excel_results(
    output_excel,
    cleaned_df,
    technical_qc,
    sample_mean,
    reference_mean,
    per_repeat,
    missing_delta,
    untreated_calibrator,
    untreated_data,
    untreated_summary,
    fig1_anova_all,
    fig1_tukey_all,
    treatment_calibrator,
    treatment_data,
    treatment_summary,
    fig2_anova_all,
    fig2_tukey_all,
):
    """
    Save all intermediate and final analysis tables.
    """

    print(
        "Saving Excel results..."
    )


    with pd.ExcelWriter(
        output_excel,
        engine="openpyxl",
    ) as writer:


        cleaned_df.to_excel(
            writer,
            sheet_name="Filtered_Raw_Data",
            index=False,
        )


        technical_qc.to_excel(
            writer,
            sheet_name="Technical_QC",
            index=False,
        )


        sample_mean.to_excel(
            writer,
            sheet_name="Sample_CT_Mean",
            index=False,
        )


        reference_mean.to_excel(
            writer,
            sheet_name="Reference_CT_Mean",
            index=False,
        )


        per_repeat.to_excel(
            writer,
            sheet_name="PerRepeat_DeltaCt",
            index=False,
        )


        untreated_calibrator.to_excel(
            writer,
            sheet_name="Fig1_Calibrator",
            index=False,
        )


        untreated_data.to_excel(
            writer,
            sheet_name="Fig1_PerRepeat",
            index=False,
        )


        untreated_summary.to_excel(
            writer,
            sheet_name="Fig1_Summary",
            index=False,
        )


        fig1_anova_all.to_excel(
            writer,
            sheet_name="Fig1_ANOVA",
            index=False,
        )


        fig1_tukey_all.to_excel(
            writer,
            sheet_name="Fig1_Tukey_HSD",
            index=False,
        )


        treatment_calibrator.to_excel(
            writer,
            sheet_name="Fig2_Clip_Calibrator",
            index=False,
        )


        treatment_data.to_excel(
            writer,
            sheet_name="Fig2_All_PerRepeat",
            index=False,
        )


        treatment_summary.to_excel(
            writer,
            sheet_name="Fig2_All_Summary",
            index=False,
        )


        fig2_plot_data = (
            treatment_data[
                treatment_data[
                    "treatment"
                ].isin(
                    TREATMENTS_FOR_APHID_PLOT
                )
            ]
            .copy()
        )


        fig2_plot_data.to_excel(
            writer,
            sheet_name="Fig2_Plot_PerRepeat",
            index=False,
        )


        fig2_plot_summary = (
            treatment_summary[
                treatment_summary[
                    "treatment"
                ].isin(
                    TREATMENTS_FOR_APHID_PLOT
                )
            ]
            .copy()
        )


        fig2_plot_summary.to_excel(
            writer,
            sheet_name="Fig2_Plot_Summary",
            index=False,
        )


        fig2_anova_all.to_excel(
            writer,
            sheet_name="Fig2_ANOVA",
            index=False,
        )


        fig2_tukey_all.to_excel(
            writer,
            sheet_name="Fig2_Tukey_HSD",
            index=False,
        )


        if not missing_delta.empty:

            missing_delta.to_excel(
                writer,
                sheet_name="Missing_DeltaCt",
                index=False,
            )


    print(
        f"Excel results saved to: "
        f"{output_excel}"
    )


# =====================================================================
# Main workflow
# =====================================================================

def main():

    args = parse_args()


    input_file = Path(
        args.input_file
    )


    output_dir = Path(
        args.output_dir
    )


    figure_dir = (
        output_dir
        / "figures"
    )


    output_excel = (
        output_dir
        / "rt-qPCR_analysis_results.xlsx"
    )


    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    figure_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    # -------------------------------------------------------------
    # Input data
    # -------------------------------------------------------------

    raw_df = load_input_data(
        input_file
    )


    cleaned_df = clean_input_data(
        raw_df
    )


    # -------------------------------------------------------------
    # Technical replicate averaging and Delta Ct
    # -------------------------------------------------------------

    (
        technical_qc,
        sample_mean,
        reference_mean,
        per_repeat,
        per_repeat_valid,
        missing_delta,
    ) = calculate_per_repeat_delta_ct(
        cleaned_df
    )


    # -------------------------------------------------------------
    # Delta-Delta Ct analyses
    # -------------------------------------------------------------

    (
        untreated_data,
        untreated_calibrator,
    ) = calculate_untreated_ddct(
        per_repeat_valid
    )


    (
        treatment_data,
        treatment_calibrator,
    ) = calculate_treatment_ddct(
        per_repeat_valid
    )


    untreated_summary = (
        make_summary(
            untreated_data
        )
    )


    treatment_summary = (
        make_summary(
            treatment_data
        )
    )


    # -------------------------------------------------------------
    # Statistics
    # -------------------------------------------------------------

    fig1_anova = {}
    fig1_significance = {}

    fig2_anova = {}
    fig2_significance = {}


    for ref_gene in (
        REFERENCE_GENES
    ):

        (
            fig1_anova[
                ref_gene
            ],
            fig1_significance[
                ref_gene
            ],
        ) = significance_figure1(
            per_repeat_valid,
            ref_gene,
        )


        (
            fig2_anova[
                ref_gene
            ],
            fig2_significance[
                ref_gene
            ],
        ) = significance_figure2(
            per_repeat_valid,
            ref_gene,
        )


        print(
            f"\nFigure 1 ANOVA "
            f"- {ref_gene}"
        )

        print(
            fig1_anova[
                ref_gene
            ]
        )


        print(
            f"\nFigure 1 Tukey HSD "
            f"- {ref_gene}"
        )

        print(
            fig1_significance[
                ref_gene
            ]
        )


        print(
            f"\nFigure 2 ANOVA "
            f"- {ref_gene}"
        )

        print(
            fig2_anova[
                ref_gene
            ]
        )


        print(
            f"\nFigure 2 Tukey HSD "
            f"- {ref_gene}"
        )

        print(
            fig2_significance[
                ref_gene
            ]
        )


    # -------------------------------------------------------------
    # Combine statistical output
    # -------------------------------------------------------------

    fig1_anova_all = pd.concat(
        [
            fig1_anova[
                ref_gene
            ]
            for ref_gene
            in REFERENCE_GENES
        ],
        ignore_index=True,
    )


    fig1_tukey_all = pd.concat(
        [
            fig1_significance[
                ref_gene
            ]
            for ref_gene
            in REFERENCE_GENES
        ],
        ignore_index=True,
    )


    fig2_anova_all = pd.concat(
        [
            fig2_anova[
                ref_gene
            ]
            for ref_gene
            in REFERENCE_GENES
        ],
        ignore_index=True,
    )


    fig2_tukey_all = pd.concat(
        [
            fig2_significance[
                ref_gene
            ]
            for ref_gene
            in REFERENCE_GENES
        ],
        ignore_index=True,
    )


    # -------------------------------------------------------------
    # Save tables
    # -------------------------------------------------------------

    save_excel_results(
        output_excel=output_excel,
        cleaned_df=cleaned_df,
        technical_qc=technical_qc,
        sample_mean=sample_mean,
        reference_mean=reference_mean,
        per_repeat=per_repeat,
        missing_delta=missing_delta,
        untreated_calibrator=untreated_calibrator,
        untreated_data=untreated_data,
        untreated_summary=untreated_summary,
        fig1_anova_all=fig1_anova_all,
        fig1_tukey_all=fig1_tukey_all,
        treatment_calibrator=treatment_calibrator,
        treatment_data=treatment_data,
        treatment_summary=treatment_summary,
        fig2_anova_all=fig2_anova_all,
        fig2_tukey_all=fig2_tukey_all,
    )


    # -------------------------------------------------------------
    # Generate figures
    # -------------------------------------------------------------

    print(
        "\nGenerating figures..."
    )


    for ref_gene in (
        REFERENCE_GENES
    ):

        plot_untreated(
            raw_data=untreated_data,
            summary_data=untreated_summary,
            ref_gene=ref_gene,
            significance_df=(
                fig1_significance[
                    ref_gene
                ]
            ),
            figure_dir=figure_dir,
        )


        plot_treatments(
            raw_data=treatment_data,
            summary_data=treatment_summary,
            ref_gene=ref_gene,
            significance_df=(
                fig2_significance[
                    ref_gene
                ]
            ),
            figure_dir=figure_dir,
        )


    combine_figure2_panels(
        reference_genes=REFERENCE_GENES,
        treatment_data=treatment_data,
        treatment_summary=treatment_summary,
        fig2_significance=fig2_significance,
        figure_dir=figure_dir,
    )


    # -------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------

    print(
        "\n"
        + "=" * 70
    )

    print(
        "Analysis completed successfully."
    )

    print(
        "=" * 70
    )


    print(
        "\nExcel output:"
    )

    print(
        output_excel
    )


    print(
        "\nFigure directory:"
    )

    print(
        figure_dir
    )


    print(
        "\nStatistics:"
    )

    print(
        "Figure 1: one-way ANOVA on biological-replicate "
        "Delta Ct values followed by Tukey HSD across untreated genotypes."
    )

    print(
        "Figure 2: within each genotype, one-way ANOVA on "
        "biological-replicate Delta Ct values followed by Tukey HSD "
        "among clip-cage control, M. persicae infestation, "
        "and R. padi infestation."
    )


    print(
        "\nSignificance thresholds:"
    )

    print(
        "** P < 0.01"
    )

    print(
        "*  P < 0.05"
    )

    print(
        "ns P >= 0.05"
    )


if __name__ == "__main__":
    main()
