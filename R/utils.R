# Internal helpers for turning supported inputs into a samples x taxa matrix
# and a metadata data.frame whose rownames are sample IDs.

.is_phyloseq <- function(x) inherits(x, "phyloseq")
.is_se <- function(x) methods::is(x, "SummarizedExperiment")

# Return a numeric matrix with samples as rows and taxa as columns
.sample_matrix <- function(x, taxa_are_rows = TRUE, assay_name = NULL) {
  if (.is_phyloseq(x)) {
    .require("phyloseq")
    otu <- phyloseq::otu_table(x)
    m <- methods::as(otu, "matrix")
    if (phyloseq::taxa_are_rows(otu)) m <- t(m)
    return(m)
  }
  if (.is_se(x)) {
    .require("SummarizedExperiment")
    if (is.null(assay_name)) assay_name <- 1L
    m <- as.matrix(SummarizedExperiment::assay(x, assay_name))
    return(t(m))
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix/data.frame, a phyloseq object, ",
      "a SummarizedExperiment, or a dist object.", call. = FALSE)
  }
  m <- if (taxa_are_rows) t(x) else x
  if (is.null(rownames(m))) {
    stop("Samples must be named (column names when `taxa_are_rows = TRUE`).",
      call. = FALSE)
  }
  if (any(m < 0, na.rm = TRUE)) {
    stop("Abundances must be non-negative.", call. = FALSE)
  }
  if (anyNA(m)) stop("Abundances contain missing values.", call. = FALSE)
  m
}

# Metadata for a supported input, falling back on what the object carries
.sample_metadata <- function(x, metadata = NULL, sample_col = NULL) {
  if (is.null(metadata)) {
    if (.is_phyloseq(x)) {
      .require("phyloseq")
      metadata <- methods::as(phyloseq::sample_data(x), "data.frame")
    } else if (.is_se(x)) {
      .require("SummarizedExperiment")
      metadata <- as.data.frame(SummarizedExperiment::colData(x))
    } else {
      stop("`metadata` is required for this type of input.", call. = FALSE)
    }
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!is.null(sample_col)) {
    .check_cols(metadata, sample_col)
    rownames(metadata) <- as.character(metadata[[sample_col]])
  }
  metadata
}

# Align metadata rows to the labels of a dist object
.align_metadata <- function(metadata, samples) {
  missing <- setdiff(samples, rownames(metadata))
  if (length(missing) > 0) {
    stop(length(missing), " sample(s) have no metadata row, e.g. '",
      missing[1], "'. Metadata rownames (or `sample_col`) must match sample ",
      "names.", call. = FALSE)
  }
  metadata[samples, , drop = FALSE]
}

.check_cols <- function(df, cols) {
  cols <- cols[!is.null(cols)]
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop("Column(s) not found in metadata: ",
      paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this input/method. ",
      "Install it with BiocManager::install('", pkg, "').", call. = FALSE)
  }
}
