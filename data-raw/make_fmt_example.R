# Simulate the synthetic FMT dataset shipped with akeR (fmt_example).
#
# Run from the package root:  Rscript data-raw/make_fmt_example.R
#
# Design
#   - 2 cohorts (sites), each with 3 donors and 12 recipients (4 per donor).
#   - Recipients are sampled before FMT (day 0) and on days 7, 30 and 90.
#   - Taxa have a shared baseline abundance plus a cohort-specific shift, so
#     samples from the same cohort are more alike. This creates a "cohort
#     effect" that can make recipients look closer to their own donor even
#     before FMT when donors from other cohorts are used for comparison.
#   - Recipients start with a depleted (less rich) microbiome. After FMT
#     they acquire donor taxa: about 60% of absent donor taxa on day 7 for
#     responders, fading to about 40% by day 90. Non-responders (1 in 4)
#     acquire about 10%.
#   - Reads are drawn from a negative binomial around each sample's
#     proportions, with sequencing depth varying between samples.

set.seed(2025)

n_taxa <- 200
taxa <- sprintf("Species_%03d", seq_len(n_taxa))
cohorts <- c("site_A", "site_B")
days <- c(0, 7, 30, 90)

# Shared log-abundance and prevalence of each taxon
base_log <- rnorm(n_taxa, 0, 2)
prevalence <- rbeta(n_taxa, 0.8, 0.8)

# Cohort-specific shifts in abundance and prevalence
cohort_log <- sapply(cohorts, function(c) rnorm(n_taxa, 0, 0.8))
cohort_prev <- sapply(cohorts, function(c) {
  plogis(qlogis(pmin(pmax(prevalence, 0.01), 0.99)) + rnorm(n_taxa, 0, 1))
})

# A new individual: log-abundances and presence/absence
new_individual <- function(cohort, richness = 1) {
  list(
    log = base_log + cohort_log[, cohort] + rnorm(n_taxa, 0, 1.2),
    present = runif(n_taxa) < cohort_prev[, cohort] * richness
  )
}

# Proportions from an individual
to_prop <- function(ind) {
  p <- exp(ind$log) * ind$present
  p / sum(p)
}

samples <- list()
meta <- list()
add_sample <- function(id, ind, row) {
  samples[[id]] <<- to_prop(ind)
  meta[[id]] <<- data.frame(sample_id = id, row, stringsAsFactors = FALSE)
}

donor_n <- 0
recipient_n <- 0
for (cohort in cohorts) {
  for (k in 1:3) {
    donor_n <- donor_n + 1
    donor_id <- sprintf("D%02d", donor_n)
    donor <- new_individual(cohort)
    donor_sample <- paste0(donor_id, "_donor")
    add_sample(donor_sample, donor, data.frame(subject_id = donor_id,
      cohort = cohort, role = "donor", donor_id = donor_id, day = NA,
      responder = NA, reference = NA_character_))

    for (r in 1:4) {
      recipient_n <- recipient_n + 1
      subject <- sprintf("R%02d", recipient_n)
      responder <- r != 4
      baseline <- new_individual(cohort, richness = 0.6)

      for (day in days) {
        ind <- baseline
        # Day-to-day variation within a person
        ind$log <- ind$log + rnorm(n_taxa, 0, 0.4)
        if (day > 0) {
          engraft <- if (responder) {
            c(`7` = 0.6, `30` = 0.5, `90` = 0.4)[[as.character(day)]]
          } else {
            0.1
          }
          gained <- donor$present & !baseline$present &
            runif(n_taxa) < engraft
          lost <- baseline$present & !donor$present &
            runif(n_taxa) < engraft / 2
          ind$present <- (baseline$present | gained) & !lost
          # Gained taxa take the donor's abundance; shared taxa move towards it
          ind$log[gained] <- donor$log[gained] + rnorm(sum(gained), 0, 0.4)
          shared <- baseline$present & donor$present
          ind$log[shared] <- (1 - engraft) * ind$log[shared] +
            engraft * donor$log[shared]
        }
        add_sample(sprintf("%s_d%02d", subject, day), ind,
          data.frame(subject_id = subject, cohort = cohort,
            role = "recipient", donor_id = donor_id, day = day,
            responder = responder, reference = donor_sample))
      }
    }
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
metadata$timepoint <- factor(
  ifelse(metadata$role == "donor", "donor", paste0("day_", metadata$day)),
  levels = c("donor", paste0("day_", days))
)
rownames(metadata) <- metadata$sample_id
metadata <- metadata[, c("sample_id", "subject_id", "cohort", "role",
  "donor_id", "day", "timepoint", "responder", "reference")]

fmt_example <- list(counts = counts, metadata = metadata)
save(fmt_example, file = "data/fmt_example.rda", compress = "xz")
