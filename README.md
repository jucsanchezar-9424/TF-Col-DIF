# TF-Col-DIF

Reproducible code and supplementary materials for the TF-Col differential item functioning study across clinical and socioeconomic subgroups in Colombian oncology.

## Associated manuscript

> Sánchez-Arévalo JC, López SA, Sierra-Matamoros FA. Differential item functioning
> of the TF-Col financial toxicity measure across clinical and socioeconomic
> subgroups in Colombian oncology. *Manuscript in preparation.*

## Related publication

Analysis code accompanying:

> Sánchez-Arévalo JC, López SA, Sierra-Matamoros FA. Differential item functioning
> of the TF-Col financial toxicity measure across clinical and socioeconomic
> subgroups in Colombian oncology. *Submitted.*

Companion to the instrument development and validation study:

> Sánchez-Arévalo JC, Sierra-Matamoros FA, Pérez Sosa NM, Cuéllar Rivera DI,
> Sánchez Pedraza R. Development and psychometric validation of the TF-Col:
> a patient- and caregiver-reported measure of financial toxicity in oncology.
> *Value in Health Regional Issues* 2026;101695.
> https://doi.org/10.1016/j.vhri.2026.101695

---

## Data

The participant-level dataset is **not distributed here**: it contains sensitive
patient information and is governed by the confidentiality conditions of the
participating institutions. De-identified data may be requested from the
corresponding author with appropriate institutional approval.

**One sheet.** The analysis reads a single sheet, `ANALISIS`, of the workbook
`BDTFCol_ANALISIS.xlsx`. That sheet is sheet `V4` of the source workbook — 348
records, the 20 items, the grouping variables and the summed scores — with two
columns appended from sheet `theta_TRI`, `theta_EAP` and `SE_theta`, merged on
the participant key (`Record ID` ↔ `ID`). Nothing was recoded, renamed,
reordered or imputed; the merge was verified to be one-to-one and the three
summed scores were checked to agree row by row across the two sheets.

Two sheets of the source workbook are easy to confuse: `V3` also has 348 rows
and the same 20 items, and `V1`/`V2` contain only the 200 adult respondents.
`ANALISIS` is the only sheet the analysis uses, and the script stops if the data
it receives lack `theta_EAP` and `SE_theta`.

### Two ways to point the script at the data

Through environment variables, never an absolute path in the code:

```bash
export TFCOL_DATA="/path/to/BDTFCol_ANALISIS.xlsx"
export TFCOL_OUT="./output"
```

Or, if the path does not resolve conveniently on your machine, import the sheet
by hand into an object named `datos` and source the script. Section 2 detects
the object and skips the file read altogether:

```r
datos <- readxl::read_excel("BDTFCol_ANALISIS.xlsx", sheet = "ANALISIS")
Sys.setenv(TFCOL_OUT = "./output")
source("DIF_TFCol_Analysis.R")
```

## Run order

A single script does everything:

```bash
export TFCOL_DATA="/path/to/BDTFCol_ANALISIS.xlsx"
export TFCOL_OUT="./output"
Rscript DIF_TFCol_Analysis.R
```

It produces, in order: item-level DIF for the six comparisons (Table 2), the
domain-level sensitivity analysis, differential test functioning (Table S1 /
Online Resource 6), the four figures (Online Resources 2–5), Table 1, the
trait distribution by group (Table S2 / Online Resource 7), and
`SessionInfo.txt`.

## Requirements

R ≥ 4.3, with `lordif`, `mirt`, `readxl`, `dplyr` and `MASS`. `mirt` ships as a
dependency of `lordif`. The script checks for all five at startup and stops with
the exact `install.packages()` call if any is missing. Figures use base graphics
only — no `ggplot2` dependency.

## Analytical decisions worth knowing before reading the code

**Two thresholds, one criterion.** `lordif` applies an internal purification
threshold of ΔR² ≥ 0.020 when deciding which items to drop from the anchor set.
The prespecified criterion for *meaningful* DIF in this study is the published
Jodoin and Gierl value of ΔR² ≥ 0.035. In the full-scale analysis four items
reached 0.020 — IF07 and IPS08 under respondent type/site, IF10 and IPS10 under
insurance regime — and **none** reached 0.035. Those four are potential signals,
not DIF.

**The one value above the criterion.** In the domain-level sensitivity analysis
(Section 5), IPS10 under insurance regime reaches ΔR² = 0.036 within the
Psychosocial Impact domain — 0.001 above the criterion, and the only value above
it anywhere in the study. It is the same item and the same comparison that
produced the largest full-scale value (0.032), estimated on a ten-item
calibration whose trait estimate is less precise than the twenty-item one. It is
reported as found; it is not treated as overturning the full-scale result, and
it is not smoothed away.

**Matching variable.** `lordif` does not accept an external matching variable: it
calibrates a graded response model on the item set it is given and conditions on
the EAP trait estimate, re-derived at each purification iteration. Any variable
computed as a "total score" and passed alongside is not used.

**One metric for DTF.** Section 6 uses a single calibration of
the full analytic sample for every comparison, with no per-comparison rescaling,
so that all comparisons and the SE-derived reference value share one scale.
SE(θ) from the validation calibration (column `SE_theta`) is linearly mapped
onto that metric before any comparison against it, and the script prints the
correlation between the two θ estimates as a check on the link.

**Impact is not bias.** Section 8.1 reports the trait and summed-score
distribution of every comparison group. Groups differ in how much financial
toxicity they actually experience, which is impact and is exactly what the DIF
models condition on; it carries no implication that any item is biased. The
table exists because the two are routinely conflated, and because a null DIF
result means something quite different depending on whether the groups being
compared differ on the construct at all.

**What the all-items-free model is not.** It is a conservative sensitivity
analysis — a stress test of possible cumulative score distortion. It is not an
upper bound: item-level differences can cancel and the procedure also absorbs
sampling variability.

**Locale.** Column names and category labels in the source workbook contain
accented characters. Matching them against literal UTF-8 strings silently
produces empty subgroups when the session locale is not UTF-8, so every script
matches on ASCII-only substrings and asserts that no group is empty (Section 2.1).

## Reproducibility notes

- One seed (`set.seed(2025)`) is set at the top and re-set before every model fit.
- A built-in check verifies that ΔR²(total) ≥ ΔR²(uniform) in all 120 cells, which
  is guaranteed by model nesting and therefore catches extraction errors.
- Figures are written with `cairo_pdf()` so that Greek symbols embed correctly
  regardless of viewer or locale; a `pdf()` + plotmath fallback is included.
- No absolute paths anywhere.

## Still to be completed before submission

- [x] Run the script end to end under R 4.5.2 / lordif 0.4.2 and reconcile every
      number in the manuscript against the exported tables.
- [ ] Confirm that θ renders correctly in Online Resources 4 and 5; on the run
      of record `cairo_pdf` warned that it could not load the cairo library.
- [ ] Archive the power/precision simulation (script, seed, replications,
      raw output) or leave the corresponding argument resting on the
      published sample-size guidance, as it currently does.
- [ ] Clean the duplicated diagnosis label in the source data
      (`Cáncer de piel` vs `Cáncer piel`).
- [ ] Completed COSMIN reporting checklist (version 2.0) as Online Resource 1.
- [ ] Supply the repository URL, which appears twice in the manuscript.

## License

Code released under the MIT License. See `LICENSE`.
