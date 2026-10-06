# RT-qPCR Delta-Delta Ct Analysis

## Overview

This script analyses RT-qPCR Ct data using the Delta-Delta Ct method and generates publication-style figures for basal gene expression and aphid-responsive expression.

The analysis is separated into two parts:

1. **Basal expression under normal growth conditions**, using untreated plants.
2. **Expression responses following aphid infestation**, using genotype-specific clip-cage controls.

The script was developed for barley thionin expression analysis using two reference genes, Actin and Ubiquitin, but the workflow can be adapted to other genes, plant species, treatments, or reference genes.

---

## Method Description

Raw Ct values are separated into target-gene and reference-gene measurements according to the `type` column.

For each genotype, treatment, biological repeat, and reference gene, technical replicates are first averaged separately for the target gene and reference gene. One Delta Ct value is then calculated for each biological replicate as:

```text
Delta Ct = mean Ct_target - mean Ct_reference
```

Statistical analyses are performed using these biological-replicate Delta Ct values.

Two separate Delta-Delta Ct analyses are then performed.

### 1. Basal expression analysis

Only untreated plants are included in the basal-expression analysis.

For each reference gene, untreated `Akashinriki` is used as the calibrator. The mean Delta Ct across biological replicates of untreated Akashinriki is used as the baseline:

```text
DeltaDelta Ct =
Delta Ct_sample - mean Delta Ct_Akashinriki untreated
```

Relative expression is calculated as:

```text
2^-DeltaDeltaCt
```

and the plotted log2-transformed relative expression is:

```text
log2(2^-DeltaDeltaCt) = -DeltaDeltaCt
```

Differences among the four barley genotypes are assessed using one-way ANOVA on biological-replicate Delta Ct values, followed by Tukey's HSD test for all pairwise comparisons.

### 2. Aphid-response analysis

For aphid-response analysis, each genotype is analysed separately.

The corresponding genotype-specific `Clip-cage control` is used as the calibrator:

```text
DeltaDelta Ct =
Delta Ct_sample - mean Delta Ct_genotype-specific clip-cage control
```

Relative expression and log2-transformed relative expression are calculated as described above.

For each genotype, differences among:

- `Clip-cage control`
- `M.persicae infestation`
- `R.padi infestation`

are assessed using one-way ANOVA on biological-replicate Delta Ct values, followed by Tukey's HSD test for all pairwise comparisons.

Only the two aphid treatments are displayed as bars in the final aphid-response figures. Statistical annotations indicate comparisons of each aphid treatment with the corresponding clip-cage control, together with the comparison between *M. persicae* and *R. padi*.

---

## Features

- Technical replicates averaged before biological-replicate analysis
- Delta Ct calculation at the biological-replicate level
- Separate basal-expression and aphid-response analyses
- Untreated Akashinriki used as the basal-expression calibrator
- Genotype-specific clip-cage controls used for aphid-response analysis
- Delta-Delta Ct relative expression analysis
- Support for multiple reference genes
- One-way ANOVA using biological-replicate Delta Ct values
- Tukey HSD post-hoc testing
- Biological-replicate mean ± SEM visualization
- Individual biological replicates shown with distinct marker shapes
- Automatic significance annotation
- Excel output containing intermediate calculations and statistical results
- Publication-style figure generation
- Combined Actin and Ubiquitin aphid-response figure
- Italicized aphid species names in figure legends

---

## Input Requirements

### Input Excel file (`--input_file`)

The input file must contain the following columns:

| Column | Description |
|---|---|
| `genotype` | Plant genotype or accession name |
| `CT` | RT-qPCR Ct value |
| `reference_gene` | Reference gene used for normalization |
| `repeat` | Biological replicate ID |
| `treatment` | Experimental treatment |
| `type` | Either `sample` or `reference_gene` |

The expected treatments are:

```text
Untreated control
Clip-cage control
M.persicae infestation
R.padi infestation
```

The current analysis uses four barley genotypes:

```text
Akashinriki
HOR10350
HOR21599
Morex
```

and two reference genes:

```text
Actin
Ubiquitin
```

Example input:

```text
genotype        CT      reference_gene  repeat  treatment                  type
Akashinriki     23.927  Actin           1       Clip-cage control          sample
Akashinriki     23.907  Actin           1       Clip-cage control          reference_gene
Akashinriki     23.297  Actin           1       R.padi infestation         sample
Akashinriki     21.737  Actin           1       R.padi infestation         reference_gene
```

Multiple technical replicate rows can be provided for the same genotype, reference gene, biological repeat, treatment, and sample type.

---

## Installation

Install the required Python packages using:

```bash
pip install pandas numpy matplotlib scipy openpyxl
```

Required packages:

- pandas
- numpy
- matplotlib
- scipy
- openpyxl

---

## Usage

Run the script from the command line:

```bash
python rt_qpcr_analysis.py \
    --input_file input.xlsx \
    --output_dir output_directory
```

Example:

```bash
python rt_qpcr_analysis.py \
    --input_file input.xlsx \
    --output_dir rt_qpcr_results
```

---

## Output

The script generates:

```text
rt-qPCR_analysis_results.xlsx
```

and a `figures/` directory.

### Excel output

The Excel workbook contains the following sheets:

| Sheet | Description |
|---|---|
| `Filtered_Raw_Data` | Cleaned and filtered input Ct data |
| `Technical_QC` | Mean, number, and SD of technical replicate Ct values |
| `Sample_CT_Mean` | Mean target-gene Ct for each biological replicate |
| `Reference_CT_Mean` | Mean reference-gene Ct for each biological replicate |
| `PerRepeat_DeltaCt` | Biological-replicate Delta Ct values |
| `Fig1_Calibrator` | Untreated Akashinriki calibrator values |
| `Fig1_PerRepeat` | Biological-replicate basal-expression results |
| `Fig1_Summary` | Summary statistics for basal expression |
| `Fig1_ANOVA` | One-way ANOVA results for basal expression |
| `Fig1_Tukey_HSD` | Tukey HSD pairwise comparisons among genotypes |
| `Fig2_Clip_Calibrator` | Genotype-specific clip-cage calibrator values |
| `Fig2_All_PerRepeat` | Biological-replicate aphid-response calculations |
| `Fig2_All_Summary` | Summary statistics for aphid-response analysis |
| `Fig2_Plot_PerRepeat` | Biological-replicate values displayed in aphid-response figures |
| `Fig2_Plot_Summary` | Summary values displayed in aphid-response figures |
| `Fig2_ANOVA` | One-way ANOVA results within each genotype |
| `Fig2_Tukey_HSD` | Tukey HSD treatment comparisons within each genotype |
| `Missing_DeltaCt` | Records with incomplete target/reference Ct data, if present |

---

## Figures

### Figure 1: Basal expression

Separate figures are generated for Actin and Ubiquitin normalization:

```text
Actin_Figure1_Untreated_ANOVA_Tukey.tif
Ubiquitin_Figure1_Untreated_ANOVA_Tukey.tif
```

Each bar represents the mean log2-transformed relative expression:

```text
log2(2^-DeltaDeltaCt)
```

and error bars represent the standard error of the mean (SEM) across biological replicates.

Individual biological replicates are shown using different symbols:

```text
Repeat 1 = circle
Repeat 2 = square
Repeat 3 = triangle
```

Genotypes are ordered according to their mean expression values.

All pairwise genotype comparisons are displayed using Tukey HSD significance brackets.

---

### Figure 2: Aphid-response expression

Separate aphid-response figures are generated for each reference gene:

```text
Actin_Figure2_Aphid_ANOVA_Tukey.tif
Ubiquitin_Figure2_Aphid_ANOVA_Tukey.tif
```

Only *M. persicae* and *R. padi* infestation treatments are displayed as bars.

Each bar represents mean log2-transformed relative expression relative to the corresponding genotype-specific clip-cage control, with SEM shown as error bars.

Individual biological replicates are shown as black symbols.

For each genotype:

- significance labels above individual aphid-treatment bars indicate comparison with the corresponding clip-cage control;
- brackets between the two aphid treatments indicate the comparison between *M. persicae* and *R. padi*.

Genotypes are ordered according to their mean response across the two aphid treatments.

---

### Combined aphid-response figure

A combined two-panel figure is also generated:

```text
Figure2_Combined_AB.tif
```

where:

```text
A = Actin
B = Ubiquitin
```

A shared legend is used for the two aphid treatments.

---

## Statistical Analysis

Statistical analyses are performed using biological-replicate Delta Ct values rather than Delta-Delta Ct or transformed relative-expression values.

### Basal expression

For each reference gene:

```text
One-way ANOVA:
Delta Ct ~ genotype
```

followed by Tukey's HSD test for all pairwise genotype comparisons.

### Aphid-response analysis

For each genotype and reference gene:

```text
One-way ANOVA:
Delta Ct ~ treatment
```

using:

```text
Clip-cage control
M.persicae infestation
R.padi infestation
```

followed by Tukey's HSD test for all pairwise treatment comparisons.

Significance is displayed as:

```text
**  P < 0.01
*   P < 0.05
ns  P >= 0.05
```

---

## Calculation Summary

### Delta Ct

```text
Delta Ct = mean Ct_target - mean Ct_reference
```

### Basal-expression DeltaDelta Ct

```text
DeltaDelta Ct =
Delta Ct_sample -
mean Delta Ct_Akashinriki untreated
```

### Aphid-response DeltaDelta Ct

```text
DeltaDelta Ct =
Delta Ct_sample -
mean Delta Ct_genotype-specific clip-cage control
```

### Relative expression

```text
Relative expression = 2^-DeltaDeltaCt
```

### Log2-transformed relative expression

```text
log2(2^-DeltaDeltaCt) = -DeltaDeltaCt
```

---

## Notes

Technical replicates are not treated as independent observations. They are averaged first to produce a single target-gene Ct and reference-gene Ct value for each biological replicate.

All statistical testing is therefore performed using biological replicates as the experimental units.

Delta Ct values are used for ANOVA and Tukey HSD testing, whereas Delta-Delta Ct-derived values are used for visualization and biological interpretation.

The basal-expression and aphid-response analyses use different calibrators because they address different biological questions. Untreated Akashinriki is used to compare basal expression among genotypes, whereas genotype-specific clip-cage controls are used to quantify aphid-induced expression changes within each genotype.
