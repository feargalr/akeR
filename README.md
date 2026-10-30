# akeR

Calculate microbiome distances to reference samples, or within and between
groups, in R.

One function, `akeR()`, returns a tidy data.frame of distances with the
samples' metadata attached, ready to plot and analyse:

- **`mode = "reference"`**: the distance from each sample to its reference
  sample, such as its FMT donor or the subject's own baseline. Your metadata
  needs a column (by default called `reference`) holding the sample ID of
  each sample's reference. Distances to the other reference samples are
  included too, as a baseline for comparison.
- **`mode = "group"`**: the distance between every pair of samples,
  labelled as within or between groups (family, subject, site, donor...).

It supports compositional (Aitchison, robust Aitchison), abundance-based
(Bray-Curtis, Hellinger, Jensen-Shannon, Manhattan, ...), presence/absence
(Jaccard, Sorensen) and phylogenetic (UniFrac) distances, and accepts count
matrices, `phyloseq` objects, `(Tree)SummarizedExperiment` objects (e.g. from
curatedMetagenomicData) or precomputed `dist` objects.

akeR deliberately does not run statistical tests, because the right test
depends on the study design. The vignettes explain the main pitfall:
pairwise distances are not independent.

## Installation

```r
# install.packages("remotes")
remotes::install_github("feargalr/akeR", build_vignettes = TRUE)
```

## Quick start

### Distance to FMT donor

```r
library(akeR)
library(ggpubr)
data(fmt_example)  # simulated FMT study

# fmt_example$metadata$reference holds each recipient's donor sample ID
fmt <- akeR(fmt_example$counts, fmt_example$metadata,
  mode = "reference", method = "jaccard")

ggboxplot(fmt, x = "timepoint", y = "distance", fill = "comparison")
```

### Within vs between groups

```r
data(family_example)  # simulated mothers and infants

fam <- akeR(family_example$counts, family_example$metadata,
  mode = "group", column = "family", method = "jaccard")
fam <- fam[fam$role_1 != fam$role_2, ]  # mother-infant pairs only

ggboxplot(fam, x = "comparison", y = "distance", add = "jitter")
```

## Functions

| Function | Purpose |
|---|---|
| `akeR()` | Distances to reference samples, or within vs between groups |
| `akeR_methods()` | List the available distance metrics |
| `calc_dist()` | Distance matrix on its own |

## Vignettes

- `vignette("akeR")`: getting started, choosing a metric, within vs between
  groups.
- `vignette("fmt")`: distance to donor in FMT studies.

## Example data

- `fmt_example`: a simulated FMT study with two sites, donors, responders
  and non-responders (`data-raw/make_fmt_example.R`).
- `family_example`: a simulated set of 12 mothers and their infants
  (`data-raw/make_family_example.R`).

Both are simulated. The vignettes also show how to run akeR on real
studies downloaded from
[curatedMetagenomicData](https://bioconductor.org/packages/curatedMetagenomicData/).
