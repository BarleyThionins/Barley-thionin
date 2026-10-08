# RT-qPCR Delta-Delta Ct Analysis

## Overview

This script analyzes RT-qPCR Ct data using the Delta-Delta Ct method and generates figures that separate expression differences among genotypes from aphid-induced expression changes within each genotype.

The analysis contains four panels for each reference gene:

| Panel | Analysis | Calibrator |
|---|---|---|
| A | Expression among genotypes under untreated conditions | Lowest-expression genotype within this treatment and reference gene |
| B | Expression among genotypes under *M. persicae* infestation | Lowest-expression genotype within this treatment and reference gene |
| C | Expression among genotypes under *R. padi* infestation | Lowest-expression genotype within this treatment and reference gene |
| D | Aphid-induced expression changes within each genotype | Corresponding genotype-specific clip-cage control |

The script was developed for barley thionin expression analysis with Actin and Ubiquitin as reference genes. The current configuration uses four genotypes and three biological replicates. Other experimental designs require changes to the configuration and associated validation and plotting code.

## Method Description

### Technical replicates and Delta Ct

Raw Ct values are separated into target-gene and reference-gene measurements using the `type` column. Technical replicates are averaged separately for each genotype, treatment, biological replicate, reference gene, and sample type.

One Delta Ct value is calculated for each biological replicate:

```text
DeltaCt = mean target Ct - mean reference Ct
```

Technical replicates are not treated as independent observations. Statistical tests use biological-replicate Delta Ct values.

### Panels A-C: Expression differences among genotypes

Genotypes are compared separately under untreated conditions, *M. persicae* infestation, and *R. padi* infestation, for each reference gene.

Within each treatment and reference gene, the genotype with the largest mean Delta Ct is selected as the lowest-expression calibrator. Its mean Delta Ct across valid biological replicates is used as a common baseline for all genotypes in that comparison:

```text
calibrator genotype = genotype with the largest mean DeltaCt
DeltaDeltaCt = sample DeltaCt - mean DeltaCt of the calibrator genotype
```

If genotype means are tied, the first genotype in `genotype_order` is selected.

The calibrator genotype has a mean plotted log2 expression of zero. Its individual biological replicates retain their variation and are not individually set to zero. This is a common mean calibrator, not a replicate-matched calibrator.

Because A-C use separately selected calibrators, bar heights across these panels cannot directly quantify induction between treatments. Panel D quantifies induction relative to clip-cage controls.

### Panel D: Aphid-induced expression changes

Each genotype and reference gene is analyzed separately. The mean Delta Ct of its corresponding clip-cage control is used as the calibrator:

```text
DeltaDeltaCt = sample DeltaCt - mean DeltaCt of the genotype-specific clip-cage control
```

Clip-cage control, *M. persicae* infestation, and *R. padi* infestation are included in the statistical analysis. Only the two aphid treatments are displayed as bars. Clip-cage control letters are shown inside panel D.

Positive log2 values indicate increased expression relative to the clip-cage control; negative values indicate decreased expression.

### Relative expression

For all panels:

```text
Relative expression = 2^(-DeltaDeltaCt)
log2 relative expression = -DeltaDeltaCt
```

Bars show the mean log2 relative expression, with SEM calculated across biological replicates. This is the mean of log2-transformed replicate values, not the log2 of the arithmetic mean relative expression.

## Input Requirements

The input Excel file must contain these columns:

| Column | Description |
|---|---|
| `genotype` | Plant genotype or accession name |
| `CT` | RT-qPCR Ct value |
| `reference_gene` | Reference gene used for normalization |
| `repeat` | Biological replicate ID: 1, 2, or 3 in the current configuration |
| `treatment` | Experimental treatment |
| `type` | `sample` for target-gene measurements or `reference_gene` for reference-gene measurements |

Expected treatment labels:

```text
Untreated control
Clip-cage control
M.persicae infestation
R.padi infestation
```

Configured genotypes:

```text
Akashinriki
HOR10350
HOR21599
Morex
```

Configured reference genes:

```text
Actin
Ubiquitin
```

Example input:

| genotype | CT | reference_gene | repeat | treatment | type |
|---|---|---|---|---|---|
| Akashinriki | 23.927 | Actin | 1 | Clip-cage control | sample |
| Akashinriki | 23.907 | Actin | 1 | Clip-cage control | reference_gene |
| Akashinriki | 23.297 | Actin | 1 | R.padi infestation | sample |
| Akashinriki | 21.737 | Actin | 1 | R.padi infestation | reference_gene |

Multiple technical-replicate rows can share the same genotype, reference gene, biological replicate, treatment, and sample type.

The script standardizes selected spelling variants and aliases, removes missing or nonnumeric Ct values, and filters records to configured genotypes, reference genes, and treatments. Biological replicates without a valid target/reference pairing are excluded from calculations and recorded in `Missing_DeltaCt`.

All four genotypes must have valid Delta Ct data for each A-C treatment/reference-gene comparison. Panel D requires a valid clip-cage calibrator for every configured genotype/reference-gene combination. Provide all three D treatment groups to obtain complete treatment comparisons. Missing or nonfinite pairwise P values produce `NA` rather than a significance letter for that comparison family.

## Installation

```bash
pip install pandas numpy matplotlib scipy openpyxl
```

## Usage

```bash
python rt_qpcr_analysis.py --input_file input.xlsx --output_dir rt_qpcr_results
```

The command-line arguments remain `--input_file` and `--output_dir`. Personal file paths do not need to be edited in the script.

## Output

The output directory contains `rt-qPCR_analysis_results.xlsx` and a `figures/` directory.

### Excel workbook

| Sheet | Description |
|---|---|
| `Filtered_Raw_Data` | Cleaned and filtered Ct records |
| `Technical_QC` | Mean Ct, valid technical-replicate count, and technical-replicate SD |
| `PerRepeat_DeltaCt` | Mean target/reference Ct, counts, SDs, and Delta Ct for each biological replicate |
| `ABC_Calibrators` | Selected calibrator genotype, mean Delta Ct, SD, and replicate count for each A-C panel/reference gene |
| `ABC_PerRepeat` | Biological-replicate A-C calculations, including calibrator identity, DeltaDeltaCt, relative expression, and log2 expression |
| `ABC_Summary` | A-C summary statistics for each genotype, reference gene, and treatment |
| `ABC_ANOVA` | One-way ANOVA results among genotypes for each A-C panel/reference gene |
| `ABC_Tukey_HSD` | Pairwise genotype comparisons for each A-C panel/reference gene |
| `ABC_CLD` | Compact letter displays for A-C genotype comparisons |
| `D_Clip_Calibrators` | Genotype-specific clip-cage mean Delta Ct, SD, and replicate count |
| `D_PerRepeat` | Biological-replicate D calculations, including clip-cage controls |
| `D_Summary` | D summary statistics, including clip-cage controls |
| `D_ANOVA` | One-way ANOVA results among treatments within each genotype/reference gene |
| `D_Tukey_HSD` | Pairwise treatment comparisons within each genotype/reference gene |
| `D_CLD` | Compact letter displays for all three D treatments |
| `Missing_DeltaCt` | Biological-replicate records without a valid target/reference pairing |

### Figures

Each reference gene has four individual panels and one combined ABCD figure. Each figure is saved in PNG, TIFF, and PDF formats.

| Reference gene | Individual figure names | Combined figure name |
|---|---|---|
| Actin | `ACTIN_A`, `ACTIN_B`, `ACTIN_C`, `ACTIN_D` | `ACTIN_Combined_ABCD` |
| Ubiquitin | `UBI_A`, `UBI_B`, `UBI_C`, `UBI_D` | `UBI_Combined_ABCD` |

For example, `ACTIN_A.png`, `ACTIN_A.tif`, and `ACTIN_A.pdf` contain the Actin-normalized untreated panel. The default configuration generates 30 figure files.

Bars represent mean log2 relative expression +/- SEM. Black symbols identify biological replicates:

| Biological replicate | Symbol |
|---|---|
| 1 | Circle |
| 2 | Square |
| 3 | Triangle |

A-C genotypes are ordered from low to high mean log2 expression within each panel and reference gene. D genotypes are ordered from low to high mean response across the two aphid treatments.

## Statistical Analysis

Statistical tests use biological-replicate Delta Ct values. Delta-Delta Ct and relative-expression values are used for visualization.

For each A-C treatment and reference gene:

```text
One-way ANOVA: DeltaCt ~ genotype
```

For each D genotype and reference gene:

```text
One-way ANOVA: DeltaCt ~ treatment
```

D includes clip-cage control and both aphid treatments. All pairwise comparisons are calculated using Tukey HSD, with the Tukey-Kramer standard error when group sizes differ. Adjusted P < 0.05 defines a significant pairwise difference.

### Compact letter displays

Figures use letters instead of significance brackets and stars:

- Groups sharing at least one letter have no statistically significant difference at adjusted P < 0.05.
- Groups with no shared letter differ significantly at adjusted P < 0.05.
- For example, `ab` shares a letter with both `a` and `b`; groups labeled `a` and `b` differ significantly.
- Shared letters do not establish equivalence.
- `NA` indicates that the comparison family lacks a complete set of finite pairwise P values.

A-C letters compare genotypes only within the same panel and reference gene. D letters compare treatments only within the same genotype and reference gene. Letters must not be compared across these separate analysis families.

The Tukey result tables also retain `significance` fields (`**`, `*`, and `ns`); figures display compact letters.

## Notes

- A-C and D address different questions: expression differences among genotypes and treatment-induced changes within genotypes, respectively.
- Selecting a common mean calibrator shifts all A-C log2 expression values equally within a comparison and does not change the Delta Ct ANOVA or Tukey results.
- Each reference gene is analyzed independently; normalization results are not pooled.
- The 2^(-DeltaDeltaCt) calculation assumes approximately equal target and reference amplification efficiencies near 100%.
- One-way ANOVA treats biological replicates as independent observations. Repeat IDs identify plotted symbols; the model does not include a repeat or blocking term.
- To adapt the script, update `genotype_order`, `reference_gene_order`, treatment labels, and associated display settings. Changing biological replicate IDs also requires updating validation, markers, labels, and plotting offsets.
