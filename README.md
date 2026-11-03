# EMT-ERphagy-Paper-DiMonaco
Scripts used in Di Monaco et al., 2026
# Analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23882537.svg)](https://doi.org/10.5281/zenodo.23882537)

Analysis code for:

> **FAM134A-dependent ER-phagy restrains procollagen secretion during Epithelial-Mesenchymal Transition**
> Di Monaco M, Vascon A, Sen O, Saksono G, Teo JY, Brunner V, Schellhaas C,
> Ulrich S, Squires I, Filipe Soares R, Peteri M, Corn J,
> Liang JR.
> *Developmental Cell* (TBD_YEAR). DOI: TBD_PAPER_DOI

Lead contact: Jin Rui Liang (jliang001@dundee.ac.uk), MRC Protein
Phosphorylation and Ubiquitylation Unit, University of Dundee.

## Contents

- [What is here](#what-is-here)
- [Quick start](#quick-start)
- [Points a reader should know before reusing this code](#points-a-reader-should-know-before-reusing-this-code)
- [Software](#software)
- [Deposited data](#deposited-data)
- [What is not here](#what-is-not-here)
- [1. Image quantification](#1-image-quantification)
- [2. Genome-wide CRISPRa screen](#2-genome-wide-crispra-screen)
- [3. Whole-cell-lysate proteomics](#3-whole-cell-lysate-proteomics)
- [Citation](#citation)
- [Licence](#licence)

## What is here

Three independent analyses. None reads another's output, and they can be run in
any order.

| | Scripts | Detail |
|---|---|---|
| **Image quantification** | `quantify_puncta.py`, `puncta_core.py`, `split_gfp/`, `transwell/` | [§1](#1-image-quantification) |
| **Genome-wide CRISPRa screen** | `Screen_GO_simplified_bubbleplots.R` | [§2](#2-genome-wide-crispra-screen) |
| **Whole-cell-lysate proteomics** | `proteomics_clustering_and_GO.R`, `proteomics_differential_expression_and_GO.R` | [§3](#3-whole-cell-lysate-proteomics) |

Each R script writes a `sessionInfo_*.txt` recording the package versions it
ran under.

## Quick start

```bash
# image quantification
python quantify_puncta.py --images /path/to/raw --out results --qc
```

```r
# screen GO enrichment
source("Screen_GO_simplified_bubbleplots.R")

# proteomics, from the directory holding data_from_mass_spec_facility/
source("proteomics_clustering_and_GO.R")
source("proteomics_differential_expression_and_GO.R")
```

Fiji macros are run from `Plugins → Macros → Run…`, and prompt for their input
and output folders.

## Points a reader should know before reusing this code

Each is documented in full in the relevant section below; they are collected
here because each one changes how a result should be read.

- **The facility's `log2FC` sign is inverted for the two ConA contrasts.** The
  proteomics differential-expression script corrects it on import. Any
  reanalysis starting from the deposited CSV files must apply the same
  correction, or every ConA effect will be inverted. See
  [§3](#important-log2fc-sign-correction).
- **The number of proteomics clusters was fixed, not selected.** `Mclust` is
  called with `G = 15` rather than choosing the number by BIC, and
  `set.seed(123)` fixes cluster *numbering* — the manuscript refers to
  particular clusters by number.
- **The screen's GO term size cap was raised** from clusterProfiler's default
  of 500 to 1000, after observing that MAPK cascade (777 genes) fell outside
  the default. This is a post-hoc adjustment, not a default setting.
- **Pearson's r in the split-GFP analysis is computed over whole cell
  regions**, not restricted to thresholded pixels, so it is inflated by
  co-absent background. Manders' M1/M2, reported alongside, are the
  better-behaved statistic.
- **Per-cell imaging tables were filtered by hand** after the scripts ran, to
  remove a small number of cell-segmentation merge artifacts. Re-running the
  code on the raw images will not reproduce the published per-cell tables
  exactly.
- **Imaging statistics pool cells across replicates** without a field- or
  replicate-level random effect, so the per-cell n is not fully independent.
  Effect directions and magnitudes are unaffected; p-values are
  anticonservative.

## Software

| | Version |
|---|---|
| Python | 3.12.2 (Windows, 64-bit); package requirements in `requirements.txt` |
| R | 4.4.2 (screen), 4.5.1 (proteomics) |
| clusterProfiler | 4.14.6 (screen), 4.16.0 (proteomics) |
| org.Hs.eg.db | 3.20.0 |
| GOSemSim | 2.32.0 |
| mclust | 6.1.2 |
| ggplot2 | 4.0.2 |
| ImageJ | 1.54p (Fiji distribution), with the Bio-Formats importer |
| DIA-NN | 2.0 (Academia) |
| limma | 3.99.0 |

Per-analysis package lists are given in each section. `requirements.txt`
records the versions the Python code is known to run under. Two version
bounds are load-bearing rather than incidental:

- **ggplot2 ≥ 4.0.** The proteomics scripts call
  `scale_fill_continuous(palette = …)`, which only accepts a palette argument
  from ggplot2 4.0 onwards; on 3.x those calls error.
- **scikit-image < 2.0.** `puncta_core.py` passes `min_size` to
  `remove_small_objects` in `detect_puncta_split` and `segment_nuclei`. That
  argument is removed in 2.0, and its replacement drops objects smaller than
  *or equal to* the threshold where the present one drops only smaller ones —
  so a punctum of exactly 6 px², or a nucleus of exactly 150 px², would be
  dropped rather than kept.

## Deposited data

| Data | Accession |
|---|---|
| CRISPRa screen, raw sequencing | GEO **GSE315931** |
| Mass spectrometry proteomics | ProteomeXchange/PRIDE **PXD073014** |

Raw microscopy fields are not deposited and are available from the lead contact
on request.

## What is not here

The differential expression itself — protein identification and quantification
(DIA-NN v2.0) and the moderated *t*-tests (limma v3.99.0) — was performed by
the proteomics facility; the proteomics scripts consume its per-comparison
output. The genome-wide CRISPRa screen was analysed with
[ScreenProcessing](https://github.com/mhorlbeck/ScreenProcessing) before the
script here takes over. LAMP2/mCherry-KDEL colocalisation was measured
interactively in Imaris and has no script; its settings are given in the
manuscript STAR Methods.

---

# 1. Image quantification

| Script | Produces |
|---|---|
| `quantify_puncta.py` + `puncta_core.py` | GR-COL4A1 reporter puncta counts and mCherry-only fractions, per field and per cell |
| `split_gfp/Split_GFP_per_cell_count_lsm.ijm`, `..._czi.ijm` | per-cell split-GFP colocalisation (Manders' M1/M2, Pearson's r) |
| `transwell/count_nuclei_transwell_manual.ijm` | blinded nuclei counts for transwell invasion indices |

The two split-GFP macros are the same analysis for different input formats. The
`.czi` version accepts several extensions at once through its `fileExts` array
if a folder holds a mixture.

## Running the puncta quantification

```bash
# process a folder, taking the condition from the filename
python quantify_puncta.py --images /path/to/raw \
    --condition-regex '^(?P<condition>[^_]+)_(?P<replicate>REP\d+)' \
    --out results --qc
```

Reads `.lsm`, `.czi` and `.tif`. Writes `per_image.csv` and `per_cell.csv`, and
with `--qc` a detection overlay per field.

`--condition-regex` is applied to each filename; its named group `condition` is
required and `replicate` is optional. With no regex, every file in the folder
is treated as one condition.

Alternatively `--manifest` takes a CSV naming the files to process and how to
label them, which is useful when filenames do not carry the condition or when a
subset is wanted:

```
filename,condition,replicate,batch,include,note
field_001.czi,sgNT,1,,yes,
field_002.czi,sgFAM134A,1,,yes,
```

Only rows marked `include` are processed, and conditions and replicate labels
come from the manifest rather than the filename.

### Detection parameters

Defaults are the published values. All are exposed as flags so that a
sensitivity analysis needs no edit to the code.

| Flag / constant | Value | |
|---|---|---|
| `--k-rfp` | 20 × noise | mCherry detection threshold |
| `--k-gfp` | 8 × noise | eGFP colocalisation threshold |
| `--min-area` | 6 px² | minimum punctum area |
| `--max-area` | 250 px² | above this, a component is split rather than kept whole |
| `--min-cell-area` | 300 px² | minimum cell footprint for per-cell output |
| denoise σ | 1.5 px | in `puncta_core.py` |
| background σ | 15 px | in `puncta_core.py` |
| blob-split peak separation | 4 px | in `puncta_core.py` |
| DAPI σ | 3 px | in `puncta_core.py` |
| minimum nucleus area | 150 px² | in `puncta_core.py` |
| nucleus seed separation | 30 px | in `puncta_core.py` |
| cytoplasm σ | 10 px | in `puncta_core.py` |
| foreground threshold | 0.25 × Otsu | in `puncta_core.py` |

## Acquisition

Single confocal slices, three channels (mCherry/RFP, eGFP, DAPI), pixel size
70.58 nm throughout. The value is hard-coded as `PIXEL_SIZE_NM` in
`puncta_core.py`, and as `DEFAULT_PIXEL_SIZE_UM` in the Fiji macros, where it
is used only if Bio-Formats cannot read calibration from the file.

Channel order is positional in both pipelines and should be checked against any
new exports:

| | mCherry/RFP | eGFP | DAPI |
|---|---|---|---|
| `puncta_core.py` (`CH_IDX`, 0-based) | 0 | 1 | 2 |
| Fiji macros (1-based) | 1 | 2 | 3 |

The Fiji macros print a per-file channel-means line and warn if DAPI is not the
brightest channel, which is the quickest check that the indices are right for a
given file.

No membrane marker was acquired, so both pipelines approximate whole-cell
territories from nuclei — by an intensity-driven marker-controlled watershed in
Python, and by fixed-radius expansion of each nuclear outline in Fiji.

## How the puncta pipeline works

The reporter carries tandem mCherry and eGFP on endogenous COL4A1. eGFP is
quenched in the acidic lysosomal lumen while mCherry is not, so mCherry-only
puncta represent COL4A1 that has reached the lysosome by ER-phagy, and
double-positive puncta a pre-degradative, non-lysosomal pool.

1. Each channel is lightly denoised (Gaussian, σ = 1.5 px). A wider Gaussian
   blur of the **denoised** image (σ = 15 px) is then used as a local
   background estimate and subtracted from it. Note the order — background is
   estimated from the denoised image, not the raw one.
2. A per-image, per-channel noise estimate is taken from the negative residual
   values only: `median(|negative residuals|) / 0.6745`. Using only the
   negative tail avoids the downward bias real puncta would introduce, since
   puncta push the residual positive.
3. mCherry puncta are connected components exceeding 20× that estimate, of area
   6–250 px². Components above 250 px² are fused puncta rather than background,
   and are split by a local-maxima watershed within that component's footprint
   (intensity smoothed σ = 1.0 px, minimum peak separation 4 px, sub-regions
   below 2 px² discarded). A component that cannot be split into at least two
   maxima is kept whole, so no above-threshold signal is discarded.
4. A punctum is scored as colocalised if its mean eGFP residual exceeds 8× the
   eGFP noise estimate.
5. Nuclei: DAPI smoothed (σ = 3 px), three-class multi-Otsu keeping the lower
   cut, small-object removal at 150 px², binary closing, then a
   distance-transform watershed (distance map smoothed σ = 4 px, minimum seed
   separation 30 px). Multi-Otsu rather than plain Otsu because a single
   bimodal threshold calibrates to "bright versus everything else" when a field
   contains one or two much brighter nuclei, dropping most dim but real nuclei
   from the mask.
6. Cells: marker-controlled watershed seeded from the nuclear masks, on the
   combined smoothed (σ = 10 px) eGFP+mCherry signal, foreground at 25% of that
   channel's Otsu value, with nuclear pixels always included. Regions below
   300 px² are excluded from per-cell output as segmentation fragments; they
   are retained at field level.
7. Each punctum is assigned to the cell whose territory contains its centroid.

`detect_puncta_split` is the detector used throughout. The older
`detect_puncta` is retained in `puncta_core.py` for reference only: it discards
any connected component larger than `max_area` instead of splitting it, which
removes every punctum inside a merged blob — around 31% and 40% of all
above-threshold signal pixels in two representative dense fields.
`per_image.csv` reports `n_blobs_split` and `n_subpuncta_recovered` so the
contribution of the splitting step is visible per field.

### Threshold sensitivity

The two detection multipliers were set during development by visual comparison
of detection and classification overlays against raw fields. Their influence on
the conclusions was tested by re-running detection and classification across a
±25% range of each (mCherry 15–25×, eGFP 6–10×) on the complete sg-knockout
panel: 54 fields, four independent biological replicates, n = 676–714 cells
depending on threshold.

Both significant effects held at every threshold tested (Dunnett's test,
p < 0.01 throughout): reduced mCherry-only fraction on FAM134A knockdown, and
reduced total puncta on P4HA2 knockdown. The P4HA2 mCherry-only fraction, not
significant at the chosen thresholds (p = 0.150), does cross into significance
at the upper end of the mCherry range and the lower end of the eGFP range; no
claim depends on it being null.

The analysis is reproduced by varying `--k-rfp` and `--k-gfp`.

## Split-GFP colocalisation

Nuclei are segmented from DAPI (Gaussian σ = 2 px, Otsu, watershed), filtered
to 8–30 µm equivalent diameter, with edge-touching nuclei excluded. Each
nucleus is then expanded by a fixed 10 µm radius to approximate a cell region.
eGFP and mCherry are median-filtered (radius 1 px) and background-subtracted
(rolling ball, radius 25 px), and thresholded **once per field** by Otsu on the
whole image rather than per cell, because a single cell's pixel sample is too
small for a stable threshold.

Reported metrics are Manders' M1 and M2 and Pearson's r, per cell. Two points
on interpretation:

- **eGFP puncta counts are deliberately not reported.** In the first batch,
  punctum count correlated with the per-field Otsu threshold at Spearman
  r = −0.85 (p < 0.0001), so counts were dominated by threshold instability
  rather than punctum density. Manders' coefficients and Pearson's r showed no
  such dependence (|r| < 0.15, p > 0.5 against both threshold and nucleus
  count), which is why those are reported.
- **Pearson's r is computed over the whole expanded cell region** of the
  background-subtracted channels, not restricted to thresholded pixels. A
  whole-region coefficient is inflated by co-absent background, so Manders'
  M1/M2 — computed against the shared per-field masks — are the better-behaved
  statistic for co-distribution.

Two nuclei closer together than roughly twice the expansion radius receive
overlapping regions, and a punctum in the overlap is assigned to whichever
region is tested first. QC overlays are written to `QC_overlays/` for every
field; crowded fields should be inspected rather than trusted blindly.

Requires the Bio-Formats importer and `Auto Local Threshold`, which ships with
Fiji but not with plain ImageJ 1.x.

## Transwell invasion

`count_nuclei_transwell_manual.ijm` counts stained nuclei across a folder of
transwell images, with the threshold set by hand for each image and everything
else automatic. It reads `.tif`, `.tiff`, `.png`, `.jpg` and `.jpeg`; for a
multi-slice file only the first slice is used, with a warning. Both
dark-background stains (DAPI or Hoechst, the default) and light-background
stains (crystal violet) are supported through the `darkBackground` setting,
which also determines whether object signal is taken as the mean or its
inverse.

Each image is converted to 8-bit and background-subtracted (rolling ball,
radius 50 px). A duplicate is Gaussian-blurred (σ = 1.0 px) and presented for
manual thresholding, after which the mask is watershed-split.

**Blinding.** Images are presented in a random order under coded titles of the
form "Image *k* of *n*", so thresholds are set without knowledge of the
condition. Opening and pre-processing happen with batch mode on so that
filenames never appear on screen, and the results tables are assembled only
after every image has been thresholded. Filenames are revealed solely in the
saved output. The threshold used for each image is recorded.

**Object filtering, in two stages.** `Analyze Particles` applies the area
filter (80–2,000 px², or µm² for a calibrated image), the circularity filter
(≥ 0.50) and edge exclusion; objects failing these are never enumerated.
Solidity (≥ 0.85) and signal (≥ 40 on the 8-bit scale, measured on the
background-subtracted image rather than the blurred copy) are then applied to
the survivors. The `Rejected` count in the output therefore refers only to
objects rejected at this second stage.

**Outputs**, written into the image folder: `nuclei_counts.csv` (per image:
nuclei kept, objects rejected, threshold used), `nuclei_objects.csv` (per
object: area, signal, solidity, kept or not), and a `QC/` folder containing
every image with kept objects outlined in green and rejected objects in red.

**Reproducing a run exactly.** The randomised presentation order is unseeded,
and the starting threshold offered for each image is the one accepted for the
previous image. A fresh run on the same images will therefore present them in a
different order and may settle on slightly different thresholds. The per-image
thresholds written to `nuclei_counts.csv` are the authoritative record of the
manual step; the file from the published run is not part of this deposit and is
available from the lead contact on request.

The invasion index is the count on the lower membrane surface as a percentage
of the count on the upper surface of the same chamber, quantified per field and
averaged within each biological replicate. Because invasion indices are ratio
data whose variance scales with the mean, replicate means were log-transformed
before testing and effects are reported as fold changes (geometric means).
Images from membranes damaged during scraping were excluded: 4 of 204
transwells across three independent experiments.

## Imaging reproducibility notes

**Per-cell tables were filtered by hand.** A small number of
cell-segmentation merge artifacts were removed manually from the combined
per-cell table after the scripts had run, and this filtering is not reapplied
by the code. Re-running on the raw images will therefore not reproduce the
published per-cell tables exactly.

**Statistics pool cells across replicates** without a field- or
replicate-level random effect, so the per-cell n is not fully independent.
Effect directions and magnitudes are unaffected; p-values are anticonservative.

**Fonts.** `puncta_core.py` resolves a TrueType font for QC overlay annotation
on first use, falling back to Pillow's built-in default if none is found.
Quantification does not depend on it and no numeric output changes if the
fallback is used.

---

# 2. Genome-wide CRISPRa screen

| Script | Produces |
|---|---|
| `Screen_GO_simplified_bubbleplots.R` | **Fig 2A, 2B** (GO Biological Process and Molecular Function enrichment of the screen's enhancing and inhibiting hits), **Supp Table 1** sheets `GO_BP_*` and `GO_MF_*`, and `symbol_map.csv` |

It writes `sessionInfo_screen_GO.txt` recording the package versions it ran
under.

## Input data

The script reads **one** file, and accepts either of two forms. It works out
which from the extension, so nothing else changes.

| `input_file` | What it is |
|---|---|
| `Supp_Table1-CRISPRa_starved_screen.xlsx` *(default)* | Supplementary Table 1 as published. Sheet `Gene_level_results` holds the 18,915 library genes, `Negative_controls` their 18,915 matched control pseudogenes. |
| `CRISPRa_starved_genetable_collapsed_namesfixed.txt` | ScreenProcessing's collapsed gene table for the amino-acid starvation arm — the original input. |

Both carry the same numbers, and the two routes have been checked to produce
**byte-identical** GO tables. (The only difference anywhere in the output is
the row order of `symbol_map.csv`: 82 of 18,915 rows sit in a different
position, all adjacent pairs such as `RBM4` / `RBM14-RBM4` that the two sources
sort differently. Every mapping is the same, and nothing downstream depends on
the order.) Table S1 is the default because it is self-contained: the raw gene
table is not redistributed here, since everything the script needs from it is
already published.

Set `input_file` at the top of the script to wherever your copy lives.

Three entry points, depending on how far back you want to start:

| Starting from | Where | What you still need to run |
|---|---|---|
| Raw sequencing reads | **GEO: GSE315931** | [ScreenProcessing](https://github.com/mhorlbeck/ScreenProcessing), then this script |
| Gene-level phenotypes | **Supp Table 1**, sheets `Gene_level_results` + `Negative_controls` | this script |
| Finished enrichment | **Supp Table 1**, sheets `GO_BP_*` / `GO_MF_*` | nothing |

### Warning: do not open the .txt gene table in Excel

The hCRISPRa-v2 library contains genes whose symbols Excel silently converts to
dates — `SEPT7` becomes `Sep.07`, `MARCH1` becomes `Mar.01`, and so on. The
`_namesfixed` suffix means those entries have been restored to their library
symbols. Opening the file in Excel and saving it will reintroduce the
corruption, and the affected genes will then be dropped from the enrichment
without any error being raised.

The script also drops the library's 18,915 matched negative-control `pseudo*`
rows before enrichment: they are the screen's null distribution, not candidate
hits.

### Gene symbol updating

hCRISPRa-v2 was built in 2016, so some library symbols are now aliases, and
`enrichGO` discards symbols it cannot map **silently**. The script therefore
resolves every library symbol to a current HGNC symbol first, updating only
where an alias maps to exactly one current symbol and leaving ambiguous ones
alone. The decision for every gene is written to
`screen_GO_output/supp_tables/symbol_map.csv`, and hits that still cannot be
mapped are printed at the end of the run.

## Running it

The script takes its paths from two variables at the top (`input_file`,
`output_dir`) and calls no `setwd()`, so it can be run from anywhere:

```r
source("Screen_GO_simplified_bubbleplots.R")
```

### Parameters

Three parameters change which GO terms appear, and are marked `PARAMETER` where
they are set:

| Parameter | Value | Effect |
|---|---|---|
| `go_max_size` | 1000 | Raised from clusterProfiler's default of 500 so that "MAPK cascade" (777 genes) is tested at all. Higher values let very broad parent terms dominate. |
| `phenotype_cutoff` | 0.7 | Hit threshold on the average phenotype of the strongest 3 sgRNAs. Relaxing it adds weaker hits and broadens the enrichment. |
| `n_top` | 15 | Number of most significant non-redundant terms drawn per panel. |

MAPK/ERK terms are additionally force-included if significant, even when
`simplify()` merged them into a relative. EMT terms are **only** emphasised in
bold and are never force-included, so every EMT term in Fig 2A reached the top
list on its own adjusted *p*-value.

## What the enrichment does

Gene-level phenotypes come from ScreenProcessing (Horlbeck et al., 2016;
Gilbert et al., 2014), not from this repository. For each gene the phenotype is
the mean of its three strongest sgRNAs (five per gene), averaged across two
biological replicates, comparing the enhanced and inhibited sorted populations.
Negative-control pseudogenes — 18,915 groups of randomly sampled non-targeting
sgRNAs — are excluded. Genes with an averaged phenotype ≥ 0.7 are classified as
enhancers (99 genes) and ≤ −0.7 as inhibitors (90 genes); no p-value threshold
is applied at this step.

Over-representation analysis runs on each hit list separately via `enrichGO`
against the org.Hs.eg.db annotation (GO annotations dated September 2024), for
Biological Process and Molecular Function. One-sided hypergeometric test
against all genes annotated in the respective ontology, Benjamini–Hochberg
correction, adjusted p < 0.05 taken as significant. GO terms with 10–1,000
annotated genes are tested. Redundant terms are removed with `simplify` (Wang
semantic similarity via GOSemSim, cutoff 0.7), keeping the lowest adjusted
p-value from each redundant group. Plots show the 15 most significant
non-redundant terms per list, together with any significant MAPK/ERK-related
terms.

## Screen software

| Package | Used for |
|---|---|
| clusterProfiler | `enrichGO`, `simplify` |
| org.Hs.eg.db | gene symbol annotation and alias resolution |
| AnnotationDbi | annotation interface (called with its namespace, since `dplyr::select` masks `select()`) |
| GOSemSim | semantic similarity behind `simplify()` |
| ggplot2, patchwork, stringr | figure construction |
| readxl | *optional*; only to read Supplementary Table 1 |
| ragg | *optional*; PNG export |
| systemfonts | *optional*; checks that Arial Narrow is installed |

Verified to reproduce the published tables under **R 4.4.2** with
clusterProfiler 4.14.6, org.Hs.eg.db 3.20.0, GOSemSim 2.32.0 and ggplot2 4.0.2.

The figure uses the system font **Arial Narrow** and, on macOS, the Quartz PDF
device. Where either is unavailable the script falls back to the default font
and to `cairo_pdf`, with a message; the layout is unchanged and only glyph
rendering differs slightly from the published figure. The font check itself
cannot stop the run.

---

# 3. Whole-cell-lysate proteomics

Two independent R scripts. Neither reads the other's output, and they can be
run in either order.

| Script | Produces |
|---|---|
| `proteomics_clustering_and_GO.R` | **Fig 4C** (CC enrichment dot plots for clusters of interest), **Supp Fig 4A** (standardised expression profiles of all 15 clusters), **Supp Table 2** (`ClusterData`, `GOenrich_full`, `GOenrich_simplified`) |
| `proteomics_differential_expression_and_GO.R` | **Fig 4B** (volcano plots: TGFβ vs untreated, and TGFβ+ConA vs TGFβ), **Supp Fig 4D** (BP enrichment of the TGFβ-vs-basal upregulated set) |

Each script writes a `sessionInfo_*.txt` recording the package versions it ran
under.

## Input data

Both proteomics scripts expect a directory `data_from_mass_spec_facility/`
beside them, holding the facility's per-comparison CSV files:

```
data_from_mass_spec_facility/
  CTRL.ConA_vs_CTRL.CTRL.csv
  TGFb.ConA_vs_CTRL.ConA.csv
  TGFb.ConA_vs_TGFb.CTRL.csv
  TGFb.CTRL_vs_CTRL.CTRL.csv
```

**This directory is not in the repository.** The underlying mass spectrometry
data are deposited at ProteomeXchange via PRIDE under accession **PXD073014**.

Treatment codes: `CTRL`/`TGFb` is the first factor (TGFβ, 5 ng/mL, 48 h) and
`CTRL`/`ConA` the second (concanamycin A, final 16 h), six biological
replicates per condition.

## Important: log2FC sign correction

The facility's analysis code reports `log2FC` with **the wrong sign for the two
ConA contrasts**, `CTRL.ConA_vs_CTRL.CTRL` and `TGFb.ConA_vs_TGFb.CTRL`.
`proteomics_differential_expression_and_GO.R` negates both columns on import
(the lines marked `SIGN CORRECTION`), so that a positive `log2FC` means "higher
with ConA", which is how the figures and the manuscript read it.

**Any reanalysis starting from the deposited CSVs must apply the same
correction**, or every ConA effect will appear inverted.

The clustering script is unaffected: it uses only the per-sample intensity
columns and discards the facility's `log2FC` and *p*-value columns.

## Running them

From the directory containing `data_from_mass_spec_facility/`:

```r
source("proteomics_clustering_and_GO.R")
source("proteomics_differential_expression_and_GO.R")
```

Run each top to bottom. In the differential expression script the GO enrichment
blocks operate on whichever contrast the preceding block assigned to `df`, so
re-running a block out of order will silently enrich the wrong comparison.

`set.seed(123)` precedes the `Mclust` call. This fixes cluster *numbering*, not
just numerical output: Gaussian mixture cluster labels are arbitrary, and the
manuscript refers to particular clusters by number (Cluster 8, Cluster 11). A
different seed or a different `mclust` version may permute the labels even when
the partition itself is unchanged.

Gene sets for enrichment are defined by adjusted p < 0.05 and the sign of
`log2FC`, with no fold-change cut; the volcano plots additionally apply
|log2FC| > 1 for labelling only.

## Proteomics software

Run under **R v4.5.1**.

| Package | Version | Used for |
|---|---|---|
| tidyverse | — | data handling, plotting |
| mclust | 6.1.2 | Gaussian mixture clustering, `G = 15`, `modelNames = "EII"` |
| clusterProfiler | 4.16.0 | `enrichGO`, `simplify` |
| org.Hs.eg.db | 3.20.0 | gene symbol annotation |
| GOSemSim | 2.32.0 | semantic similarity for `simplify(measure = "Wang")` |
| AnnotationDbi | — | annotation interface |
| ggplot2 | 4.0.2 | figures |
| ggrastr, ggrepel, ggprism, gridExtra | — | point rasterisation, labels, theme, panel layout |

## Known limitation

In `proteomics_differential_expression_and_GO.R` the GO background
(`universe_genes`) splits composite protein groups on `";"` while the
foreground sets (`genes_to_run`) do not. The 59 of 7,889 protein groups with a
composite `Genes` entry — 0.7%, all paralogue pairs or pseudogenes such as
`GATD3;GATD3B` — therefore cannot match a symbol in the background and are
effectively absent from the test set.

This is left exactly as it was run, so that the code reproduces the published
enrichment. No affected entry bears on the conclusions. The clustering script
splits in both places and is not affected.

## Note: input files are selected by name, not position

`proteomics_differential_expression_and_GO.R` imports the four CSVs with
`list.files()`, which sorts by the active locale's collation. These filenames
order differently in different locales, because `ConA` and `CTRL` first differ
at `o` versus `T`:

| position | English locale | C/POSIX locale |
|---|---|---|
| `[[1]]` | CTRL.ConA_vs_CTRL.CTRL | CTRL.ConA_vs_CTRL.CTRL |
| `[[2]]` | TGFb.ConA_vs_CTRL.ConA | **TGFb.CTRL_vs_CTRL.CTRL** |
| `[[3]]` | TGFb.ConA_vs_TGFb.CTRL | **TGFb.ConA_vs_CTRL.ConA** |
| `[[4]]` | TGFb.CTRL_vs_CTRL.CTRL | **TGFb.ConA_vs_TGFb.CTRL** |

The script therefore indexes the list by name:

```r
CTRL.ConA_vs_CTRL.CTRL <- diff_exprs_sets$CTRL.ConA_vs_CTRL.CTRL
```

The list is named from the filenames on import, so this resolves to the same
comparison in any locale.

Earlier versions of this script selected by position (`diff_exprs_sets[[1]]`
and so on). In an English locale — the locale the published analysis was run in
— the positional and name-based mappings are identical, so the figures in the
paper are unaffected. Under `LC_COLLATE=C` three of the four variables would
have received the wrong comparison, with no error raised, since the variable
names were asserted rather than read from the data.

---

## Citation

If you use this code, please cite the paper above and this repository:

> https://doi.org/10.5281/zenodo.23882537

That is the concept DOI: it always resolves to the latest archived version.
Version-specific DOIs are listed on the Zenodo record if you need to cite the
exact snapshot you used.

## Licence

Released under the MIT Licence; see [LICENSE](LICENSE).

---

## Before release — delete this section once complete

- `TBD_YEAR`, `TBD_PAPER_DOI` — on acceptance. Cut release v1.1.0 at the same
  time so the archived copy carries the final citation. The concept DOI
  (10.5281/zenodo.23882537) does not change and needs no update.
- Confirm `LICENSE` and `requirements.txt` are committed.
- Confirm the Zenodo record's title and author list match the paper; Zenodo
  fills both from GitHub by default.
