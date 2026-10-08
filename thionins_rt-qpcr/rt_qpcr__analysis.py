#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RT-qPCR ABCD analysis.

A-C compare genotypes separately in untreated, M. persicae, and R. padi
conditions, for each reference gene. The genotype with the largest mean
DeltaCt is the common mean calibrator for that panel/reference gene.
D compares treatments within each genotype using its mean clip-cage DeltaCt.
Technical replicates are averaged first; ANOVA/Tukey use biological DeltaCt.
Bars show mean log2 expression +/- SEM; points show biological replicates.
Shared CLD letters indicate adjusted P >= 0.05, not equivalence.

Usage: python rt_qpcr_analysis.py --input_file input.xlsx --output_dir output
Dependencies: numpy pandas scipy matplotlib openpyxl.
"""
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import f_oneway, studentized_range


genotype_order = ['Akashinriki', 'HOR10350', 'HOR21599', 'Morex']

reference_gene_order = ['Actin', 'Ubiquitin']

untreated_treatment = 'Untreated control'

clip_control_treatment = 'Clip-cage control'

mp_treatment = 'M.persicae infestation'

rp_treatment = 'R.padi infestation'

all_treatments = [untreated_treatment, clip_control_treatment, mp_treatment, rp_treatment]

treatments_for_calculation = [clip_control_treatment, mp_treatment, rp_treatment]

treatments_for_plot = [mp_treatment, rp_treatment]

treatment_display = {clip_control_treatment: 'Clip-cage control', mp_treatment: '$M.\\ persicae$', rp_treatment: '$R.\\ padi$'}

FIG1_BAR_COLOR = '#D9D9D9'

treatment_colors = {mp_treatment: '#EED1CC', rp_treatment: '#C9DCC4'}

repeat_order = [1, 2, 3]

repeat_markers = {1: 'o', 2: 's', 3: '^'}

repeat_labels = {1: 'Repeat 1', 2: 'Repeat 2', 3: 'Repeat 3'}

P_STAR = 0.05

P_DOUBLE_STAR = 0.01

GENOTYPE_ORDER = genotype_order
REFERENCE_GENES = reference_gene_order
UNTREATED_TREATMENT = untreated_treatment
CLIP_CONTROL_TREATMENT = clip_control_treatment
MP_TREATMENT = mp_treatment
RP_TREATMENT = rp_treatment
ALL_TREATMENTS = all_treatments


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
        "Akahsinriki": "Akashinriki",
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
        "UBI": "Ubiquitin",
        "ubi": "Ubiquitin",
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

def p_to_star(p):
    """
    Convert P value to significance symbols.
    """
    if pd.isna(p):
        return ''
    if p < P_DOUBLE_STAR:
        return '**'
    if p < P_STAR:
        return '*'
    return 'ns'

def calculate_treatment_ddct(input_df):
    data = input_df[input_df['treatment'].isin(treatments_for_calculation)].copy()
    calibrator = data[data['treatment'] == clip_control_treatment].groupby(['genotype', 'reference_gene']).agg(calibrator_DeltaCt=('DeltaCt', 'mean'), calibrator_SD=('DeltaCt', 'std'), calibrator_n=('DeltaCt', 'count')).reset_index()
    expected = pd.MultiIndex.from_product([genotype_order, reference_gene_order], names=['genotype', 'reference_gene'])
    observed = pd.MultiIndex.from_frame(calibrator[['genotype', 'reference_gene']])
    missing_calibrators = expected.difference(observed)
    if len(missing_calibrators) > 0:
        raise ValueError('Missing genotype-specific Clip-cage calibrators:\n' + str(list(missing_calibrators)))
    data = pd.merge(data, calibrator[['genotype', 'reference_gene', 'calibrator_DeltaCt']], on=['genotype', 'reference_gene'], how='left')
    data['DeltaDeltaCt'] = data['DeltaCt'] - data['calibrator_DeltaCt']
    data['Relative_expression'] = 2 ** (-data['DeltaDeltaCt'])
    data['log2FC'] = -data['DeltaDeltaCt']
    return (data, calibrator)

def make_summary(input_df):
    return input_df.groupby(['genotype', 'reference_gene', 'treatment']).agg(mean_log2FC=('log2FC', 'mean'), median_log2FC=('log2FC', 'median'), SD_log2FC=('log2FC', 'std'), SEM_log2FC=('log2FC', 'sem'), mean_relative_expression=('Relative_expression', 'mean'), SEM_relative_expression=('Relative_expression', 'sem'), n=('log2FC', 'count')).reset_index()

def significance_figure1(input_df, ref_gene, condition=untreated_treatment):
    data = input_df[(input_df['reference_gene'] == ref_gene) & (input_df['treatment'] == condition)].copy()
    data = data[data['genotype'].isin(genotype_order)].dropna(subset=['DeltaCt', 'genotype']).copy()
    groups = {}
    for genotype in genotype_order:
        values = data[data['genotype'] == genotype]['DeltaCt'].dropna().astype(float).values
        if len(values) > 0:
            groups[genotype] = values
    if len(groups) < 2:
        raise ValueError(f'Not enough genotype groups for Figure 1 statistics: {ref_gene}')
    group_names = [genotype for genotype in genotype_order if genotype in groups]
    group_arrays = [groups[genotype] for genotype in group_names]
    F_statistic, anova_p = f_oneway(*group_arrays)
    all_values = np.concatenate(group_arrays)
    grand_mean = np.mean(all_values)
    k = len(group_arrays)
    N = len(all_values)
    ss_between = sum((len(values) * (np.mean(values) - grand_mean) ** 2 for values in group_arrays))
    ss_within = sum((np.sum((values - np.mean(values)) ** 2) for values in group_arrays))
    df_between = k - 1
    df_within = N - k
    ms_between = ss_between / df_between if df_between > 0 else np.nan
    ms_within = ss_within / df_within if df_within > 0 else np.nan
    anova_table = pd.DataFrame([{'reference_gene': ref_gene, 'source': 'genotype', 'sum_sq': ss_between, 'df': df_between, 'mean_sq': ms_between, 'F_statistic': F_statistic, 'p_value': anova_p}, {'reference_gene': ref_gene, 'source': 'residual', 'sum_sq': ss_within, 'df': df_within, 'mean_sq': ms_within, 'F_statistic': np.nan, 'p_value': np.nan}])
    if df_within <= 0 or pd.isna(ms_within):
        raise ValueError(f'Insufficient residual degrees of freedom for Tukey HSD: {ref_gene}')
    q_critical = studentized_range.ppf(0.95, k, df_within)
    tukey_results = []
    for i in range(len(group_names)):
        for j in range(i + 1, len(group_names)):
            genotype1 = group_names[i]
            genotype2 = group_names[j]
            values1 = groups[genotype1]
            values2 = groups[genotype2]
            n1 = len(values1)
            n2 = len(values2)
            mean1 = np.mean(values1)
            mean2 = np.mean(values2)
            mean_diff = mean2 - mean1
            tukey_se = np.sqrt(ms_within / 2.0 * (1.0 / n1 + 1.0 / n2))
            if tukey_se == 0 or pd.isna(tukey_se):
                q_statistic = np.nan
                p_adjusted = np.nan
                ci_lower = np.nan
                ci_upper = np.nan
                reject_h0 = False
            else:
                q_statistic = abs(mean_diff) / tukey_se
                p_adjusted = studentized_range.sf(q_statistic, k, df_within)
                half_width = q_critical * tukey_se
                ci_lower = mean_diff - half_width
                ci_upper = mean_diff + half_width
                reject_h0 = p_adjusted < 0.05
            tukey_results.append({'genotype1': genotype1, 'genotype2': genotype2, 'comparison': f'{genotype1} vs {genotype2}', 'n_genotype1': n1, 'n_genotype2': n2, 'mean_DeltaCt_difference': mean_diff, 'q_statistic': q_statistic, 'CI_lower': ci_lower, 'CI_upper': ci_upper, 'p_value_adjusted': p_adjusted, 'reject_H0': reject_h0, 'significance': p_to_star(p_adjusted)})
    tukey_df = pd.DataFrame(tukey_results)
    return (anova_table, tukey_df)

def significance_figure2(input_df, ref_gene):
    data = input_df[input_df['reference_gene'] == ref_gene].copy()
    anova_tables = []
    tukey_tables = []
    treatment_order_for_stats = [clip_control_treatment, mp_treatment, rp_treatment]
    for genotype in genotype_order:
        genotype_data = data[(data['genotype'] == genotype) & data['treatment'].isin(treatment_order_for_stats)].dropna(subset=['DeltaCt', 'treatment']).copy()
        groups = {}
        for treatment in treatment_order_for_stats:
            values = genotype_data[genotype_data['treatment'] == treatment]['DeltaCt'].dropna().astype(float).values
            if len(values) > 0:
                groups[treatment] = values
        group_names = [treatment for treatment in treatment_order_for_stats if treatment in groups]
        if len(group_names) < 2:
            anova_tables.append(pd.DataFrame([{'reference_gene': ref_gene, 'genotype': genotype, 'source': 'treatment', 'sum_sq': np.nan, 'df': np.nan, 'mean_sq': np.nan, 'F_statistic': np.nan, 'p_value': np.nan}]))
            continue
        group_arrays = [groups[treatment] for treatment in group_names]
        F_statistic, anova_p = f_oneway(*group_arrays)
        all_values = np.concatenate(group_arrays)
        grand_mean = np.mean(all_values)
        k = len(group_arrays)
        N = len(all_values)
        ss_between = sum((len(values) * (np.mean(values) - grand_mean) ** 2 for values in group_arrays))
        ss_within = sum((np.sum((values - np.mean(values)) ** 2) for values in group_arrays))
        df_between = k - 1
        df_within = N - k
        ms_between = ss_between / df_between if df_between > 0 else np.nan
        ms_within = ss_within / df_within if df_within > 0 else np.nan
        anova_tables.append(pd.DataFrame([{'reference_gene': ref_gene, 'genotype': genotype, 'source': 'treatment', 'sum_sq': ss_between, 'df': df_between, 'mean_sq': ms_between, 'F_statistic': F_statistic, 'p_value': anova_p}, {'reference_gene': ref_gene, 'genotype': genotype, 'source': 'residual', 'sum_sq': ss_within, 'df': df_within, 'mean_sq': ms_within, 'F_statistic': np.nan, 'p_value': np.nan}]))
        if df_within <= 0 or pd.isna(ms_within):
            continue
        q_critical = studentized_range.ppf(0.95, k, df_within)
        tukey_results = []
        for i in range(len(group_names)):
            for j in range(i + 1, len(group_names)):
                treatment1 = group_names[i]
                treatment2 = group_names[j]
                values1 = groups[treatment1]
                values2 = groups[treatment2]
                n1 = len(values1)
                n2 = len(values2)
                mean1 = np.mean(values1)
                mean2 = np.mean(values2)
                mean_diff = mean2 - mean1
                tukey_se = np.sqrt(ms_within / 2.0 * (1.0 / n1 + 1.0 / n2))
                if tukey_se == 0 or pd.isna(tukey_se):
                    q_statistic = np.nan
                    p_adjusted = np.nan
                    ci_lower = np.nan
                    ci_upper = np.nan
                    reject_h0 = False
                else:
                    q_statistic = abs(mean_diff) / tukey_se
                    p_adjusted = studentized_range.sf(q_statistic, k, df_within)
                    half_width = q_critical * tukey_se
                    ci_lower = mean_diff - half_width
                    ci_upper = mean_diff + half_width
                    reject_h0 = p_adjusted < 0.05
                tukey_results.append({'reference_gene': ref_gene, 'genotype': genotype, 'treatment1': treatment1, 'treatment2': treatment2, 'comparison': f'{treatment1} vs {treatment2}', 'n_treatment1': n1, 'n_treatment2': n2, 'mean_DeltaCt_difference': mean_diff, 'q_statistic': q_statistic, 'CI_lower': ci_lower, 'CI_upper': ci_upper, 'p_value_adjusted': p_adjusted, 'reject_H0': reject_h0, 'significance': p_to_star(p_adjusted)})
        if len(tukey_results) > 0:
            tukey_tables.append(pd.DataFrame(tukey_results))
    if len(anova_tables) > 0:
        anova_df = pd.concat(anova_tables, ignore_index=True)
    else:
        anova_df = pd.DataFrame()
    if len(tukey_tables) > 0:
        tukey_df = pd.concat(tukey_tables, ignore_index=True)
    else:
        tukey_df = pd.DataFrame()
    return (anova_df, tukey_df)

FRAME_LINE_WIDTH = 2.5
plt.rcParams['font.family'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 18
TICK_LINE_WIDTH = 2.0

plt.rcParams['axes.linewidth'] = FRAME_LINE_WIDTH

TICK_FONT_SIZE = 18

AXIS_FONT_SIZE = 22

SIGNIFICANCE_FONT_SIZE = 21

LEGEND_FONT_SIZE = 18

PANEL_LABEL_FONT_SIZE = 30

def add_black_border(ax):
    for spine_name in ['left', 'right', 'top', 'bottom']:
        ax.spines[spine_name].set_visible(True)
        ax.spines[spine_name].set_color('black')
        ax.spines[spine_name].set_linewidth(FRAME_LINE_WIDTH)
    ax.tick_params(axis='both', which='both', direction='out', width=TICK_LINE_WIDTH, length=6, color='black', labelcolor='black', labelsize=TICK_FONT_SIZE)

def get_repeat_legend_handles():
    handles = []
    for repeat_id in repeat_order:
        handle = Line2D([0], [0], marker=repeat_markers[repeat_id], linestyle='None', markerfacecolor='black', markeredgecolor='black', markeredgewidth=0.7, markersize=7, label=repeat_labels[repeat_id])
        handles.append(handle)
    return handles

def compact_letters(group_names, tukey_df, first_col, second_col):
    """Cover all non-significant pairs with maximal cliques; never merge significant pairs.

    A shared letter means Tukey adjusted P >= 0.05, not proof of equivalence.
    Missing/non-finite pairwise P values yield 'NA' for this comparison family.
    """
    from itertools import combinations
    names = list(group_names)
    if len(names) < 2:
        return {name: 'NA' for name in names}
    lookup = {
        frozenset((row[first_col], row[second_col])): row['p_value_adjusted']
        for _, row in tukey_df.iterrows()
    }
    for pair in combinations(names, 2):
        p = lookup.get(frozenset(pair), np.nan)
        if not np.isfinite(p):
            return {name: 'NA' for name in names}
    cliques = []
    for size in range(1, len(names) + 1):
        for subset in combinations(names, size):
            if all(lookup[frozenset(pair)] >= P_STAR
                   for pair in combinations(subset, 2)):
                cliques.append(frozenset(subset))
    maximal = [c for c in cliques if not any(c < other for other in cliques)]
    maximal.sort(key=lambda c: tuple(i for i, name in enumerate(names) if name in c))
    labels = {name: '' for name in names}
    for index, clique in enumerate(maximal):
        letter = chr(ord('a') + index)
        for name in clique:
            labels[name] += letter
    return labels

def _draw_untreated_panel(ax, raw_data, summary_data, ref_gene, significance_df):
    raw = raw_data[raw_data['reference_gene'] == ref_gene].copy()
    summary = summary_data[summary_data['reference_gene'] == ref_gene].copy()

    plot_genotype_order = (
        summary[summary['genotype'].isin(genotype_order)]
        .set_index('genotype')
        .reindex(genotype_order)
        .sort_values('mean_log2FC', ascending=True, na_position='last')
        .index.tolist()
    )
    print(f'Figure 1 plotting order ({ref_gene}, low to high):', plot_genotype_order)

    x = np.arange(len(plot_genotype_order))

    means, sems = [], []
    for genotype in plot_genotype_order:
        temp = summary[summary['genotype'] == genotype]
        if temp.empty:
            means.append(np.nan)
            sems.append(np.nan)
        else:
            means.append(temp['mean_log2FC'].iloc[0])
            sems.append(temp['SEM_log2FC'].iloc[0])

    bars = ax.bar(
        x, means, yerr=sems, width=0.65, capsize=4,
        edgecolor='black', linewidth=1.2
    )
    for bar in bars:
        bar.set_facecolor(FIG1_BAR_COLOR)

    repeat_x_offsets = {1: -0.07, 2: 0.0, 3: 0.07}
    for i, genotype in enumerate(plot_genotype_order):
        genotype_data = raw[raw['genotype'] == genotype]
        for repeat_id in repeat_order:
            temp = genotype_data[genotype_data['repeat'] == repeat_id]
            if temp.empty:
                continue
            value = temp['log2FC'].iloc[0]
            if pd.isna(value):
                continue
            ax.scatter(
                x[i] + repeat_x_offsets[repeat_id], value,
                marker=repeat_markers[repeat_id], s=48,
                facecolor='black', edgecolor='black', linewidth=0.7,
                zorder=10
            )

    # Determine top and bottom of plotted data, including biological replicates.
    upper_values, lower_values = [], []
    for mean_value, sem_value in zip(means, sems):
        if pd.isna(mean_value):
            continue
        sem_value = 0 if pd.isna(sem_value) else sem_value
        upper_values.append(mean_value + sem_value)
        lower_values.append(mean_value - sem_value)

    raw_values = raw['log2FC'].dropna().values
    if len(raw_values) > 0:
        upper_values.extend(raw_values.tolist())
        lower_values.extend(raw_values.tolist())

    if upper_values:
        base_top = max(upper_values)
        base_bottom = min(lower_values)
    else:
        base_top, base_bottom = 1, -1

    data_range = base_top - base_bottom
    if data_range == 0:
        data_range = 1

    labels = compact_letters(plot_genotype_order, significance_df, 'genotype1', 'genotype2')
    highest_annotation = base_top
    for i, genotype in enumerate(plot_genotype_order):
        values = raw.loc[raw['genotype'] == genotype, 'log2FC'].dropna().tolist()
        if pd.notna(means[i]):
            values.append(means[i] + (0 if pd.isna(sems[i]) else sems[i]))
        if not values:
            continue
        text_y = max(0, max(values)) + data_range * 0.07
        ax.text(x[i], text_y, labels[genotype], ha='center', va='bottom',
                fontsize=SIGNIFICANCE_FONT_SIZE, fontweight='bold')
        highest_annotation = max(highest_annotation, text_y)

    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_genotype_order, fontsize=TICK_FONT_SIZE, rotation=20, ha='right')
    ax.set_xlabel('Genotype', fontsize=AXIS_FONT_SIZE)
    ax.set_ylabel(r'$\log_2(2^{-\Delta\Delta C_t})$', fontsize=AXIS_FONT_SIZE)

    add_black_border(ax)
    ax.grid(False)

    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    if y_range == 0:
        y_range = 1
    final_top = max(ymax, highest_annotation + data_range * 0.08)
    ax.set_ylim(ymin - y_range * 0.05, final_top)

    return labels

def _draw_treatments_panel(ax, raw_data, summary_data, ref_gene, significance_df):
    """Draw one Figure 2 panel on an existing axes and return legend handles/labels."""
    raw = raw_data[raw_data['reference_gene'] == ref_gene].copy()
    summary = summary_data[summary_data['reference_gene'] == ref_gene].copy()

    raw_plot = raw[raw['treatment'].isin(treatments_for_plot)].copy()
    summary_plot = summary[summary['treatment'].isin(treatments_for_plot)].copy()

    genotype_plot_mean = (
        summary_plot.groupby('genotype')['mean_log2FC']
        .mean()
        .reindex(genotype_order)
    )
    plot_genotype_order = genotype_plot_mean.sort_values(
        ascending=True, na_position='last'
    ).index.tolist()
    print(f'Figure 2 plotting order ({ref_gene}, low to high):', plot_genotype_order)

    x = np.arange(len(plot_genotype_order))
    n_treatments = len(treatments_for_plot)
    bar_width = 0.30
    offsets = (np.arange(n_treatments) - (n_treatments - 1) / 2) * bar_width

    bar_information = {}
    repeat_x_offsets = {1: -0.035, 2: 0.0, 3: 0.035}

    for j, treatment in enumerate(treatments_for_plot):
        means, sems = [], []
        for genotype in plot_genotype_order:
            temp = summary_plot[
                (summary_plot['genotype'] == genotype)
                & (summary_plot['treatment'] == treatment)
            ]
            if temp.empty:
                means.append(np.nan)
                sems.append(np.nan)
            else:
                means.append(temp['mean_log2FC'].iloc[0])
                sems.append(temp['SEM_log2FC'].iloc[0])

        xpos = x + offsets[j]
        ax.bar(
            xpos, means, yerr=sems, width=bar_width, capsize=3,
            color=treatment_colors[treatment], edgecolor='black',
            linewidth=1.1, label=treatment_display[treatment]
        )

        for i, genotype in enumerate(plot_genotype_order):
            genotype_treatment_data = raw_plot[
                (raw_plot['genotype'] == genotype)
                & (raw_plot['treatment'] == treatment)
            ].copy()
            raw_group_values = genotype_treatment_data['log2FC'].dropna().astype(float).values

            mean_value = means[i]
            sem_value = 0 if pd.isna(sems[i]) else sems[i]
            candidates = []
            if not pd.isna(mean_value):
                candidates.append(mean_value + sem_value)
            if len(raw_group_values) > 0:
                candidates.extend(raw_group_values.tolist())
            local_top = max(candidates) if candidates else np.nan

            bar_information[(genotype, treatment)] = {
                'x': xpos[i],
                'mean': mean_value,
                'sem': sems[i],
                'raw_values': raw_group_values,
                'local_top': local_top,
            }

            for repeat_id in repeat_order:
                temp = genotype_treatment_data[
                    genotype_treatment_data['repeat'] == repeat_id
                ]
                if temp.empty:
                    continue
                value = temp['log2FC'].iloc[0]
                if pd.isna(value):
                    continue
                ax.scatter(
                    xpos[i] + repeat_x_offsets[repeat_id], value,
                    marker=repeat_markers[repeat_id], s=45,
                    facecolor='black', edgecolor='black', linewidth=0.7,
                    zorder=10
                )

    all_values = []
    for info in bar_information.values():
        mean_value = info['mean']
        sem_value = 0 if pd.isna(info['sem']) else info['sem']
        if not pd.isna(mean_value):
            all_values.extend([mean_value - sem_value, mean_value + sem_value])
        all_values.extend(info['raw_values'].tolist())

    if all_values:
        base_top = max(all_values)
        base_bottom = min(all_values)
    else:
        base_top, base_bottom = 1, -1

    data_range = base_top - base_bottom
    if data_range == 0:
        data_range = 1

    control_letters = {}
    highest_annotation = base_top
    for genotype in plot_genotype_order:
        group_names = [clip_control_treatment, mp_treatment, rp_treatment]
        pairs = significance_df[significance_df['genotype'] == genotype] if not significance_df.empty else significance_df
        letters = compact_letters(group_names, pairs, 'treatment1', 'treatment2')
        control_letters[genotype] = letters[clip_control_treatment]
        for treatment in treatments_for_plot:
            info = bar_information[(genotype, treatment)]
            if pd.isna(info['local_top']):
                continue
            text_y = max(0, info['local_top']) + data_range * 0.075
            ax.text(info['x'], text_y, letters[treatment],
                    ha='center', va='bottom', fontsize=SIGNIFICANCE_FONT_SIZE,
                    fontweight='bold', zorder=20)
            highest_annotation = max(highest_annotation, text_y)

    ax.axhline(y=0, color='black', linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_genotype_order, fontsize=TICK_FONT_SIZE, rotation=20, ha='right')
    ax.set_xlabel('Genotype', fontsize=AXIS_FONT_SIZE)
    ax.set_ylabel(r'$\log_2(2^{-\Delta\Delta C_t})$', fontsize=AXIS_FONT_SIZE)

    add_black_border(ax)
    ax.grid(False)

    ymin, ymax = ax.get_ylim()
    y_range = ymax - ymin
    if y_range == 0:
        y_range = 1
    final_top = max(ymax, highest_annotation + data_range * 0.10)
    ax.set_ylim(ymin - y_range * 0.08, final_top)

    # Reserve a clear top band for the control annotation, above all data letters.
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + data_range * 0.18)
    if all(letter == 'a' for letter in control_letters.values()):
        control_note = 'Clip-cage controls: a'
    else:
        # Never label controls 'a' when the actual CLD is different or unavailable.
        control_note = 'Clip-cage controls: ' + '; '.join(
            f'{genotype}: {control_letters[genotype]}' for genotype in plot_genotype_order)
        print(f'WARNING ({ref_gene}): control letters are not all a; displaying actual letters.')
    ax.text(0.02, 0.98, control_note, transform=ax.transAxes,
            ha='left', va='top', fontsize=LEGEND_FONT_SIZE,
            bbox=dict(facecolor='white', edgecolor='none', pad=2), zorder=30,
            wrap=True)

    return ax.get_legend_handles_labels()

ABC_CONDITIONS = {'A': untreated_treatment, 'B': mp_treatment, 'C': rp_treatment}
ABC_TITLES = {'A': 'Untreated control', 'B': r'$M.\ persicae$', 'C': r'$R.\ padi$'}

def calculate_between_genotypes(input_df):
    """Lowest expression = largest mean DeltaCt; tie resolved by genotype_order.

    A common mean calibrator shifts all groups equally and does not change their
    DeltaCt ANOVA/Tukey results. It is not a replicate-matched calibrator.
    """
    rows, calibrators = [], []
    for panel, condition in ABC_CONDITIONS.items():
        for ref in reference_gene_order:
            data = input_df[(input_df.treatment == condition) &
                            (input_df.reference_gene == ref)].copy()
            means = data.groupby('genotype')['DeltaCt'].mean().reindex(genotype_order)
            if means.isna().any():
                raise ValueError(f'{panel}/{ref}: missing genotypes: {means[means.isna()].index.tolist()}')
            calibrator_genotype = means.idxmax()
            values = data.loc[data.genotype == calibrator_genotype, 'DeltaCt']
            baseline = values.mean()
            data['panel'] = panel
            data['calibrator_genotype'] = calibrator_genotype
            data['calibrator_DeltaCt'] = baseline
            data['DeltaDeltaCt'] = data.DeltaCt - baseline
            data['log2FC'] = -data.DeltaDeltaCt
            data['Relative_expression'] = np.exp2(data.log2FC)
            rows.append(data)
            calibrators.append(dict(panel=panel, treatment=condition, reference_gene=ref,
                calibrator_genotype=calibrator_genotype, calibrator_DeltaCt=baseline,
                calibrator_SD=values.std(), calibrator_n=len(values)))
            print(f'{panel}/{ref}: lowest-expression calibrator = {calibrator_genotype}')
    return pd.concat(rows, ignore_index=True), pd.DataFrame(calibrators)

def save_figure(fig, name, figure_dir):
    for extension in ('png', 'tif', 'pdf'):
        options = {'pil_kwargs': {'compression': 'tiff_lzw'}} if extension == 'tif' else {}
        fig.savefig(os.path.join(figure_dir, name + '.' + extension), dpi=300,
                    bbox_inches='tight', pad_inches=0.15, **options)

def draw_abc(ax, panel, ref, abc_data, abc_summary, abc_tests):
    condition = ABC_CONDITIONS[panel]
    raw = abc_data[abc_data.treatment == condition]
    summary = abc_summary[abc_summary.treatment == condition]
    _draw_untreated_panel(ax, raw, summary, ref, abc_tests[(panel, ref)])
    color = FIG1_BAR_COLOR if panel == 'A' else treatment_colors[condition]
    for patch in ax.patches:
        patch.set_facecolor(color)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    figure_dir = output_dir / 'figures'
    figure_dir.mkdir(parents=True, exist_ok=True)
    excel_output_path = output_dir / 'rt-qPCR_analysis_results.xlsx'
    df = clean_input_data(load_input_data(Path(args.input_file)))
    if df.empty:
        raise ValueError('No usable experimental records remain after filtering.')
    if df['repeat'].isna().any() or not df['repeat'].isin(repeat_order).all():
        raise ValueError('Biological repeat identifiers must be 1, 2, or 3; edit REPEAT_ORDER and plotting configuration for another design.')
    if not df['type'].isin(['sample', 'reference_gene']).all():
        raise ValueError('type must contain sample or reference_gene.')
    technical_qc, sample_mean, reference_mean, per_repeat, per_repeat_valid, missing_delta = calculate_per_repeat_delta_ct(df)
    if not np.isfinite(per_repeat_valid['DeltaCt']).all():
        raise ValueError('DeltaCt values must be finite.')
    treatment_data, treatment_calibrator = calculate_treatment_ddct(per_repeat_valid)
    treatment_summary = make_summary(treatment_data)
    fig2_anova, fig2_significance = {}, {}
    for ref in reference_gene_order:
        fig2_anova[ref], fig2_significance[ref] = significance_figure2(per_repeat_valid, ref)

    abc_data, abc_calibrators = calculate_between_genotypes(per_repeat_valid)
    abc_summary = make_summary(abc_data)

    abc_anova, abc_tukey, abc_letters = [], [], []

    abc_tests = {}

    for panel, condition in ABC_CONDITIONS.items():
        for ref in reference_gene_order:
            anova, tukey = significance_figure1(per_repeat_valid, ref, condition)
            abc_tests[(panel, ref)] = tukey
            abc_anova.append(anova.assign(panel=panel, treatment=condition))
            abc_tukey.append(tukey.assign(panel=panel, treatment=condition, reference_gene=ref))
            order = abc_summary[(abc_summary.treatment == condition) &
                                (abc_summary.reference_gene == ref)].sort_values('mean_log2FC').genotype.tolist()
            letters = compact_letters(order, tukey, 'genotype1', 'genotype2')
            abc_letters.extend(dict(panel=panel, treatment=condition, reference_gene=ref,
                                    genotype=g, CLD=letter) for g, letter in letters.items())

    d_letters = []

    for ref in reference_gene_order:
        for genotype in genotype_order:
            pairs = fig2_significance[ref]
            pairs = pairs[pairs.genotype == genotype]
            letters = compact_letters(treatments_for_calculation, pairs, 'treatment1', 'treatment2')
            d_letters.extend(dict(reference_gene=ref, genotype=genotype, treatment=t, CLD=l)
                             for t, l in letters.items())

    for ref in reference_gene_order:
        tag = 'ACTIN' if ref == 'Actin' else 'UBI'
        for panel in ABC_CONDITIONS:
            fig, ax = plt.subplots(figsize=(8, 7))
            draw_abc(ax, panel, ref, abc_data, abc_summary, abc_tests)
            ax.text(-0.14, 1.06, panel, transform=ax.transAxes, fontsize=PANEL_LABEL_FONT_SIZE, fontweight='bold')
            fig.tight_layout()
            save_figure(fig, f'{tag}_{panel}', figure_dir)
            plt.close(fig)
        fig, ax = plt.subplots(figsize=(14, 7))
        handles, labels = _draw_treatments_panel(ax, treatment_data, treatment_summary, ref, fig2_significance[ref])
        ax.text(-0.08, 1.04, 'D', transform=ax.transAxes, fontsize=PANEL_LABEL_FONT_SIZE, fontweight='bold')
        fig.legend(handles, labels, loc='lower center', ncol=2, frameon=False, fontsize=LEGEND_FONT_SIZE)
        fig.subplots_adjust(bottom=0.25, top=0.88, left=0.1, right=0.98)
        save_figure(fig, f'{tag}_D', figure_dir)
        plt.close(fig)
    
        fig = plt.figure(figsize=(21, 14))
        grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.1])
        axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
        for panel, ax in zip(ABC_CONDITIONS, axes):
            draw_abc(ax, panel, ref, abc_data, abc_summary, abc_tests)
            ax.text(-0.14, 1.08, panel, transform=ax.transAxes, fontsize=PANEL_LABEL_FONT_SIZE, fontweight='bold')
        ax = fig.add_subplot(grid[1, :])
        handles, labels = _draw_treatments_panel(ax, treatment_data, treatment_summary, ref, fig2_significance[ref])
        ax.text(-0.045, 1.04, 'D', transform=ax.transAxes, fontsize=PANEL_LABEL_FONT_SIZE, fontweight='bold')
        natural_handle = plt.Rectangle((0, 0), 1, 1, facecolor=FIG1_BAR_COLOR,
                                       edgecolor='black', linewidth=1.1)
        fig.legend([natural_handle] + handles, ['Untreated'] + labels,
                   loc='lower center', bbox_to_anchor=(0.5, 0.01),
                   ncol=3, frameon=False, fontsize=LEGEND_FONT_SIZE)
        fig.subplots_adjust(left=0.065, right=0.98, top=0.93, bottom=0.14, wspace=0.40, hspace=0.60)
        save_figure(fig, f'{tag}_Combined_ABCD', figure_dir)
        plt.close(fig)

    with pd.ExcelWriter(excel_output_path, engine='openpyxl') as writer:
        tables = {
            'Filtered_Raw_Data': df, 'Technical_QC': technical_qc,
            'PerRepeat_DeltaCt': per_repeat, 'ABC_Calibrators': abc_calibrators,
            'ABC_PerRepeat': abc_data, 'ABC_Summary': abc_summary,
            'ABC_ANOVA': pd.concat(abc_anova, ignore_index=True),
            'ABC_Tukey_HSD': pd.concat(abc_tukey, ignore_index=True),
            'ABC_CLD': pd.DataFrame(abc_letters), 'D_Clip_Calibrators': treatment_calibrator,
            'D_PerRepeat': treatment_data, 'D_Summary': treatment_summary,
            'D_ANOVA': pd.concat(list(fig2_anova.values()), ignore_index=True),
            'D_Tukey_HSD': pd.concat(list(fig2_significance.values()), ignore_index=True),
            'D_CLD': pd.DataFrame(d_letters), 'Missing_DeltaCt': missing_delta,
        }
        for sheet, table in tables.items():
            table.to_excel(writer, sheet_name=sheet, index=False)

    print('Completed: ACTIN/UBI individual A-D panels and two combined ABCD figures.')

    print('Figures:', figure_dir)

    print('Statistics:', excel_output_path)

    print('A-C letters compare genotypes within the same panel and reference gene.')

    print('D letters compare treatments ONLY within each genotype; controls remain hidden.')

    print('All bars: mean log2 expression +/- SEM; points: biological replicates.')


if __name__ == '__main__':
    main()
