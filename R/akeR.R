#' Distances to reference samples, or within and between groups
#'
#' The main akeR function. It calculates distances between samples and
#' returns them as a data.frame with the samples' metadata attached, ready
#' for plotting and for whatever statistical analysis suits your study
#' design.
#'
#' @section `mode = "reference"`:
#' The metadata column (default `"reference"`) gives, for each sample, the ID
#' of the sample it should be compared with. IDs must match sample names in
#' the count table. Samples with `NA` (for example the reference samples
#' themselves) are not compared. Typical uses:
#' \itemize{
#'   \item FMT: each recipient sample's reference is its donor's sample.
#'   \item Longitudinal: each sample's reference is the same subject's
#'     baseline sample.
#'   \item A single reference for everyone, such as one donor or a mock
#'     community: every row gets the same ID.
#' }
#'
#' Each sample is compared with its own reference (`comparison = "own"`) and
#' with every other reference sample in the column (`comparison = "other"`).
#' The `other` distances are a natural baseline: in FMT, a recipient who is
#' closer to their own donor than to other donors shows donor-specific
#' change rather than a general shift. Filter to `comparison == "own"` if you
#' only need the own-reference distances.
#'
#' @section `mode = "group"`:
#' The metadata column (default `"group"`) defines groups. Every pair of
#' samples is returned once and labelled `"within"` (same group) or
#' `"between"` (different groups). Every metadata column is included for
#' both samples (`<name>_1`, `<name>_2`), so pairs are easy to filter, e.g.
#' `res[res$role_1 != res$role_2, ]` for mother-infant pairs only.
#'
#' @section Statistical testing:
#' akeR does not test the distances, because the appropriate test depends on
#' the study design. Note that pairwise distances are not independent: each
#' sample appears in many pairs. A Wilcoxon or t-test on the rows of the
#' output (as done by `ggpubr::stat_compare_means()`) treats them as
#' independent and gives p-values that are too small. Consider permutation
#' tests that shuffle labels across samples, PERMANOVA
#' ([vegan::adonis2()]), or models with random effects for samples or
#' subjects.
#'
#' @param x Abundance data: a numeric matrix or data.frame of counts (taxa in
#'   rows, samples in columns; see `taxa_are_rows` in [calc_dist()]), a
#'   `phyloseq` object, a `SummarizedExperiment`/`TreeSummarizedExperiment`,
#'   or a precomputed `dist` object.
#' @param metadata Sample metadata with sample names as rownames (or see
#'   `sample_col`). Optional for phyloseq and SummarizedExperiment input,
#'   where the object's own sample data are used.
#' @param mode `"reference"` or `"group"`; see the sections below.
#' @param column Metadata column to use. Defaults to `"reference"` or
#'   `"group"` to match `mode`.
#' @param method Distance metric; see [akeR_methods()].
#' @param sample_col Metadata column holding sample names, if they are not
#'   the rownames.
#' @param ... Passed to [calc_dist()], e.g. `taxa_are_rows`, `relative`,
#'   `pseudocount`.
#' @return A data.frame. Its `"method"` attribute records the distance
#'   metric.
#'   \describe{
#'     \item{Reference mode}{One row per sample-reference pair: `sample`,
#'       `reference`, `distance`, `comparison` (`"own"` or `"other"`),
#'       followed by the sample's metadata.}
#'     \item{Group mode}{One row per pair of samples: `sample_1`,
#'       `sample_2`, `distance`, `comparison` (`"within"` or `"between"`),
#'       `group_pair` (the group for within pairs, e.g. `"A"`, or both groups
#'       for between pairs, e.g. `"A vs B"`), then `<name>_1` and `<name>_2`
#'       for every metadata column.}
#'   }
#' @export
#' @examples
#' # FMT: distance from each recipient to their own donor and other donors
#' data(fmt_example)
#' fmt <- akeR(fmt_example$counts, fmt_example$metadata, mode = "reference",
#'   method = "jaccard")
#' head(fmt[, c("sample", "reference", "distance", "comparison",
#'   "timepoint")])
#' aggregate(distance ~ timepoint + comparison, data = fmt, FUN = median)
#'
#' # Groups: mother-infant pairs within and between families
#' data(family_example)
#' fam <- akeR(family_example$counts, family_example$metadata,
#'   mode = "group", column = "family", method = "jaccard")
#' fam <- fam[fam$role_1 != fam$role_2, ]
#' aggregate(distance ~ comparison, data = fam, FUN = median)
#'
#' # Plot with ggplot2
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(fmt, aes(timepoint, distance, fill = comparison)) +
#'     geom_boxplot()
#' }
akeR <- function(x, metadata = NULL, mode = c("reference", "group"),
                 column = NULL, method = "aitchison", sample_col = NULL,
                 ...) {
  mode <- match.arg(mode)
  if (is.null(column)) column <- mode

  d <- if (inherits(x, "dist")) x else calc_dist(x, method = method, ...)
  if (is.null(attr(d, "Labels"))) {
    stop("A supplied dist object must have sample labels.", call. = FALSE)
  }
  metadata <- .sample_metadata(x, metadata, sample_col)
  .check_cols(metadata, column)
  metadata <- .align_metadata(metadata, attr(d, "Labels"))

  out <- if (mode == "reference") {
    .reference_mode(as.matrix(d), metadata, column)
  } else {
    .group_mode(d, metadata, column)
  }
  rownames(out) <- NULL
  attr(out, "method") <- if (is.null(attr(d, "method"))) "user supplied" else
    attr(d, "method")
  out
}
