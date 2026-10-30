fmt_data <- function() {
  set.seed(2)
  counts <- matrix(rpois(30 * 8, 20), 30, 8,
    dimnames = list(paste0("t", 1:30),
      c("dA", "dB", "r1_0", "r1_7", "r2_0", "r2_7", "r3_0", "r3_7")))
  metadata <- data.frame(
    row.names = colnames(counts),
    reference = c(NA, NA, "dA", "dA", "dA", "dA", "dB", "dB"),
    day = c(NA, NA, 0, 7, 0, 7, 0, 7),
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}

test_that("reference mode returns own and other distances", {
  sim <- fmt_data()
  res <- akeR(sim$counts, sim$metadata, method = "bray")
  expect_s3_class(res, "data.frame")
  expect_equal(attr(res, "method"), "bray")
  dm <- as.matrix(calc_dist(sim$counts, "bray"))
  expect_equal(res$distance, dm[cbind(res$sample, res$reference)])
  expect_equal(sum(res$comparison == "own"), 6)
  expect_equal(sum(res$comparison == "other"), 6)
  own <- res[res$comparison == "own", ]
  expect_equal(own$reference, sim$metadata[own$sample, "reference"])
  expect_true("day" %in% names(res))
  expect_false(any(c("dA", "dB") %in% res$sample))
})

test_that("reference mode works with a single reference and other columns", {
  sim <- fmt_data()
  sim$metadata$baseline <- c(NA, NA, NA, "r1_0", NA, "r2_0", NA, "r3_0")
  res <- akeR(sim$counts, sim$metadata, column = "baseline")
  expect_equal(nrow(res[res$comparison == "own", ]), 3)

  sim$metadata$single <- c(NA, NA, rep("dA", 6))
  res <- akeR(sim$counts, sim$metadata, column = "single")
  expect_true(all(res$comparison == "own"))
  expect_equal(nrow(res), 6)
})

test_that("reference mode validates IDs and columns", {
  sim <- fmt_data()
  sim$metadata$reference[3] <- "not_a_sample"
  expect_error(akeR(sim$counts, sim$metadata), "not samples")
  expect_error(akeR(sim$counts, sim$metadata, column = "nope"), "not found")
  expect_error(akeR(sim$counts, sim$metadata[-1, ]), "no metadata")
})

test_that("group mode returns every pair once with labels", {
  sim <- sim_data()
  res <- akeR(sim$counts, sim$metadata, mode = "group", method = "bray")
  n <- ncol(sim$counts)
  expect_equal(nrow(res), n * (n - 1) / 2)
  dm <- as.matrix(calc_dist(sim$counts, "bray"))
  expect_equal(res$distance, dm[cbind(res$sample_1, res$sample_2)])
  expect_equal(res$comparison == "within", res$group_1 == res$group_2)
  expect_equal(sum(res$comparison == "within"), 3 * choose(4, 2))
  expect_true(all(c("batch_1", "batch_2") %in% names(res)))
  expect_true(all(res$group_pair[res$comparison == "between"] %in%
    c("A vs B", "A vs C", "B vs C")))
})

test_that("group mode drops NA groups and accepts dist and phyloseq input", {
  sim <- sim_data()
  sim$metadata$group[1] <- NA
  res <- akeR(calc_dist(sim$counts, "jaccard"), sim$metadata,
    mode = "group")
  expect_false("s1" %in% c(res$sample_1, res$sample_2))
  expect_equal(attr(res, "method"), "jaccard")

  skip_if_not_installed("phyloseq")
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(sim$counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(sim$metadata)
  )
  expect_s3_class(akeR(ps, mode = "group"), "data.frame")
})
