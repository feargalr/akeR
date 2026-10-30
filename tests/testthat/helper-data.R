# Small simulated dataset: 3 groups x 4 samples, 30 taxa, with group signal
sim_data <- function(seed = 1) {
  set.seed(seed)
  n_taxa <- 30
  groups <- rep(c("A", "B", "C"), each = 4)
  base <- lapply(c("A", "B", "C"), function(g) stats::rgamma(n_taxa, 0.5))
  names(base) <- c("A", "B", "C")
  counts <- vapply(groups, function(g) {
    stats::rpois(n_taxa, 2000 * base[[g]] / sum(base[[g]]) * 10)
  }, numeric(n_taxa))
  colnames(counts) <- paste0("s", seq_along(groups))
  rownames(counts) <- paste0("taxon", seq_len(n_taxa))
  metadata <- data.frame(
    group = groups,
    batch = rep(c("x", "y"), 6),
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}
