.akeR_methods <- data.frame(
  method = c("aitchison", "robust.aitchison", "bray", "jaccard", "sorensen",
    "hellinger", "jsd", "manhattan", "euclidean", "canberra", "kulczynski",
    "horn", "unifrac", "wunifrac"),
  type = c("compositional", "compositional", "abundance", "presence/absence",
    "presence/absence", "abundance", "abundance", "abundance", "abundance",
    "abundance", "abundance", "abundance", "phylogenetic", "phylogenetic"),
  description = c(
    "Euclidean distance of CLR-transformed abundances (pseudocount added)",
    "Aitchison distance on CLR of each sample's non-zero taxa (no pseudocount)",
    "Bray-Curtis dissimilarity",
    "Jaccard distance on presence/absence",
    "Sorensen (binary Bray-Curtis) distance on presence/absence",
    "Euclidean distance of Hellinger-transformed abundances",
    "Jensen-Shannon distance (square root of JS divergence, log base 2)",
    "Manhattan (city block) distance",
    "Euclidean distance",
    "Canberra distance",
    "Kulczynski dissimilarity",
    "Morisita-Horn dissimilarity",
    "Unweighted UniFrac (phyloseq input with a tree)",
    "Weighted UniFrac (phyloseq input with a tree)"
  ),
  uses_relative = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE,
    TRUE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

#' Distance metrics available in akeR
#'
#' Lists the distance and dissimilarity metrics that [calc_dist()] supports,
#' with a short description of each.
#'
#' @return A data.frame with columns `method`, `type`, `description` and
#'   `uses_relative`. `uses_relative` is `TRUE` for metrics that are computed
#'   on relative abundances when `relative = TRUE` in [calc_dist()].
#' @export
#' @examples
#' akeR_methods()
akeR_methods <- function() {
  .akeR_methods
}

#' Calculate a distance matrix from microbiome abundance data
#'
#' Computes sample-to-sample distances using one of the metrics listed by
#' [akeR_methods()].
#'
#' Microbiome data are compositional, so the choice of metric matters:
#' \describe{
#'   \item{`"aitchison"`}{Centred log-ratio (CLR) transform followed by
#'     Euclidean distance. Zeros are handled with a pseudocount.}
#'   \item{`"robust.aitchison"`}{CLR calculated over the non-zero taxa of each
#'     sample, with zeros left at zero, then Euclidean distance. No
#'     pseudocount is needed. Each sample is transformed on its own; unlike
#'     `vegan::vegdist(method = "robust.aitchison")` in recent vegan versions,
#'     zeros are not imputed from the other samples, which would let samples
#'     influence each other's distances.}
#'   \item{`"bray"`, `"manhattan"`, `"euclidean"`, `"canberra"`,
#'     `"kulczynski"`}{Computed on relative abundances by default, so that
#'     differences in sequencing depth do not drive the distance.}
#'   \item{`"jaccard"`, `"sorensen"`}{Presence/absence. Sensitive to
#'     sequencing depth; consider rarefying first.}
#'   \item{`"hellinger"`}{Euclidean distance on square-root proportions.}
#'   \item{`"jsd"`}{Jensen-Shannon distance between sample profiles.}
#'   \item{`"unifrac"`, `"wunifrac"`}{Need a phyloseq object with a
#'     phylogenetic tree.}
#' }
#'
#' @param x Abundance data: a numeric matrix or data.frame (taxa in rows by
#'   default), a `phyloseq` object, or a `SummarizedExperiment` /
#'   `TreeSummarizedExperiment`.
#' @param method Distance metric; one of `akeR_methods()$method`.
#' @param taxa_are_rows For matrix input, are taxa in rows (`TRUE`, the
#'   default) or columns?
#' @param relative Convert to relative abundances before computing
#'   abundance-based metrics (see Details). Default `TRUE`.
#' @param pseudocount Pseudocount for `"aitchison"`. `"auto"` (default) uses 1
#'   for integer counts and half the smallest non-zero value otherwise.
#' @param assay_name For `SummarizedExperiment` input, the assay to use.
#'   Defaults to the first assay.
#' @return A `dist` object labelled with sample names.
#' @export
#' @examples
#' data(family_example)
#' d <- calc_dist(family_example$counts, method = "aitchison")
#' d_bray <- calc_dist(family_example$counts, method = "bray")
calc_dist <- function(x, method = "aitchison", taxa_are_rows = TRUE,
                      relative = TRUE, pseudocount = "auto",
                      assay_name = NULL) {
  method <- match.arg(method, .akeR_methods$method)

  if (method %in% c("unifrac", "wunifrac")) {
    if (!.is_phyloseq(x) || is.null(x@phy_tree)) {
      stop("UniFrac needs a phyloseq object with a phylogenetic tree.",
        call. = FALSE)
    }
    .require("phyloseq")
    return(phyloseq::UniFrac(x, weighted = method == "wunifrac"))
  }

  m <- .sample_matrix(x, taxa_are_rows = taxa_are_rows,
    assay_name = assay_name)
  m <- m[, colSums(m) > 0, drop = FALSE]
  empty <- rowSums(m) == 0
  if (any(empty)) {
    stop("Sample(s) with no counts: ",
      paste(rownames(m)[empty], collapse = ", "), call. = FALSE)
  }

  uses_relative <- .akeR_methods$uses_relative[.akeR_methods$method == method]
  if ((relative && uses_relative) || method == "jsd") {
    m <- m / rowSums(m)
  }

  d <- switch(method,
    aitchison = stats::dist(.clr(m, pseudocount)),
    robust.aitchison = stats::dist(.rclr(m)),
    jaccard = vegan::vegdist(m, method = "jaccard", binary = TRUE),
    sorensen = vegan::vegdist(m, method = "bray", binary = TRUE),
    hellinger = stats::dist(vegan::decostand(m, method = "hellinger")),
    jsd = .jsd(m),
    euclidean = stats::dist(m),
    vegan::vegdist(m, method = method)
  )
  attr(d, "method") <- method
  d
}

.clr <- function(m, pseudocount = "auto") {
  if (identical(pseudocount, "auto")) {
    pseudocount <- if (all(m == round(m))) 1 else min(m[m > 0]) / 2
  }
  if (!is.numeric(pseudocount) || pseudocount <= 0) {
    stop("`pseudocount` must be \"auto\" or a positive number.", call. = FALSE)
  }
  lm <- log(m + pseudocount)
  lm - rowMeans(lm)
}

# Robust CLR: CLR over each row's non-zero values, zeros stay zero
.rclr <- function(m) {
  lm <- log(m)
  lm[!is.finite(lm)] <- NA
  out <- lm - rowMeans(lm, na.rm = TRUE)
  out[is.na(out)] <- 0
  out
}

# Jensen-Shannon distance between rows of a proportion matrix
.jsd <- function(p) {
  n <- nrow(p)
  out <- matrix(0, n, n, dimnames = list(rownames(p), rownames(p)))
  plogp <- function(x) ifelse(x > 0, x * log2(x), 0)
  h <- rowSums(plogp(p))
  for (i in seq_len(n - 1)) {
    js <- i + seq_len(n - i)
    mid <- sweep(p[js, , drop = FALSE], 2, p[i, ], "+") / 2
    div <- rowSums(plogp(mid))
    jsd <- (h[i] + h[js]) / 2 - div
    out[js, i] <- sqrt(pmax(jsd, 0))
  }
  stats::as.dist(out)
}
