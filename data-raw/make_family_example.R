# Simulate the synthetic mother-infant dataset shipped with akeR
# (family_example).
#
# Run from the package root:  Rscript data-raw/make_family_example.R
#
# Design
#   - 12 families, each with one mother and one infant sample.
#   - Mothers have an adult-type microbiome; infants have a less rich,
#     infant-type microbiome dominated by a different set of taxa, so mothers
#     and infants differ strongly as groups.
#   - Vertical transmission: each infant carries about 25% of the taxa
#     present in their own mother, at an abundance close to the mother's.
#     This makes mother-infant pairs from the same family more alike than
#     pairs from different families.
#   - Reads are drawn from a negative binomial around each sample's
#     proportions, with sequencing depth varying between samples.

set.seed(2017)

n_taxa <- 200
taxa <- sprintf("Species_%03d", seq_len(n_taxa))
n_families <- 12

# Adult and infant community types: shared taxa pool, different prevalence
# and abundance. Taxa 1-25 are typical infant taxa.
adult_log <- rnorm(n_taxa, 0, 2)
adult_prev <- rbeta(n_taxa, 0.8, 0.8)
infant_log <- rnorm(n_taxa, -1, 2)
infant_prev <- rbeta(n_taxa, 0.4, 2)
infant_taxa <- 1:25
infant_log[infant_taxa] <- infant_log[infant_taxa] + 4
infant_prev[infant_taxa] <- rbeta(length(infant_taxa), 4, 1.5)

new_individual <- function(log_mu, prev) {
  list(
    log = log_mu + rnorm(n_taxa, 0, 1.2),
    present = runif(n_taxa) < prev
  )
}

to_prop <- function(ind) {
  p <- exp(ind$log) * ind$present
  p / sum(p)
}

samples <- list()
meta <- list()
for (f in seq_len(n_families)) {
  family <- sprintf("F%02d", f)

  mother <- new_individual(adult_log, adult_prev)
  infant <- new_individual(infant_log, infant_prev)

  # Vertical transmission from the mother
  transmitted <- mother$present & runif(n_taxa) < 0.25
  infant$present <- infant$present | transmitted
  infant$log[transmitted] <- mother$log[transmitted] - 1 +
    rnorm(sum(transmitted), 0, 0.5)

  for (role in c("mother", "infant")) {
    id <- paste0(family, "_", role)
    samples[[id]] <- to_prop(if (role == "mother") mother else infant)
    meta[[id]] <- data.frame(sample_id = id, family = family, role = role,
      stringsAsFactors = FALSE)
  }
}

props <- do.call(cbind, samples)
rownames(props) <- taxa
depth <- round(rlnorm(ncol(props), log(2e5), 0.4))
counts <- sapply(seq_len(ncol(props)), function(j) {
  rnbinom(n_taxa, mu = props[, j] * depth[j], size = 2)
})
dimnames(counts) <- dimnames(props)
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, ]

metadata <- do.call(rbind, meta)
rownames(metadata) <- metadata$sample_id

family_example <- list(counts = counts, metadata = metadata)
save(family_example, file = "data/family_example.rda", compress = "xz")
