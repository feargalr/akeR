#' Simulated faecal microbiota transplantation (FMT) data
#'
#' A synthetic dataset of FMT donors and recipients, simulated so that the
#' true effects are known. It is not real data.
#'
#' Two sites each have 3 donors and 12 recipients (4 per donor). Recipients
#' are sampled before FMT (day 0) and on days 7, 30 and 90. Recipients start
#' with a less diverse microbiome and, after FMT, acquire taxa from their
#' donor: responders acquire about 60% of the donor taxa they lack by day 7,
#' fading to about 40% by day 90; non-responders (one recipient per donor)
#' acquire about 10%. Samples from the same site share a site-specific
#' profile, which creates a cohort effect when donors from both sites are
#' compared (see `vignette("fmt", package = "akeR")`).
#'
#' @format A list with two elements:
#' \describe{
#'   \item{`counts`}{Integer matrix of simulated read counts, taxa (rows)
#'     by samples (columns).}
#'   \item{`metadata`}{data.frame with one row per sample (rownames are
#'     sample IDs):
#'     \describe{
#'       \item{`sample_id`, `subject_id`}{Sample and subject identifiers.}
#'       \item{`cohort`}{Site: `"site_A"` or `"site_B"`.}
#'       \item{`role`}{`"donor"` or `"recipient"`.}
#'       \item{`donor_id`}{The donor (for donors, themselves).}
#'       \item{`day`}{Days since FMT (`NA` for donors).}
#'       \item{`timepoint`}{Factor: `"donor"`, `"day_0"`, `"day_7"`,
#'         `"day_30"`, `"day_90"`.}
#'       \item{`responder`}{Whether the recipient was simulated as a
#'         responder (`NA` for donors).}
#'       \item{`reference`}{For recipients, the sample ID of their donor
#'         (`NA` for donors). The column [akeR()] uses by default.}
#'     }}
#' }
#' @source Simulated with `data-raw/make_fmt_example.R` in the package
#'   source.
#' @examples
#' data(fmt_example)
#' table(fmt_example$metadata$cohort, fmt_example$metadata$timepoint)
#'
#' # Distance from each recipient to their own donor and to other donors
#' fmt <- akeR(fmt_example$counts, fmt_example$metadata, mode = "reference",
#'   method = "jaccard")
#' aggregate(distance ~ timepoint + comparison, data = fmt, FUN = median)
"fmt_example"

#' Simulated mother-infant data
#'
#' A synthetic dataset of 12 families, each with one mother and one infant
#' sample, simulated so that the true effects are known. It is not real
#' data.
#'
#' Infants have a less rich microbiome dominated by different taxa from
#' adults, so mothers and infants differ strongly as groups. Each infant also
#' carries about 25% of the taxa present in their own mother, at a lower
#' abundance, so mother-infant pairs from the same family are more alike
#' than pairs from different families.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{`counts`}{Integer matrix of simulated read counts, taxa (rows)
#'     by samples (columns).}
#'   \item{`metadata`}{data.frame with one row per sample (rownames are
#'     sample IDs) and columns `sample_id`, `family` and `role`
#'     (`"mother"` or `"infant"`).}
#' }
#' @source Simulated with `data-raw/make_family_example.R` in the package
#'   source.
#' @examples
#' data(family_example)
#' head(family_example$metadata)
#'
#' # Mother-infant pairs within and between families
#' fam <- akeR(family_example$counts, family_example$metadata,
#'   mode = "group", column = "family", method = "jaccard")
#' fam <- fam[fam$role_1 != fam$role_2, ]
#' aggregate(distance ~ comparison, data = fam, FUN = median)
"family_example"
