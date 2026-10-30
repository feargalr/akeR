test_that("every listed method returns a labelled dist", {
  sim <- sim_data()
  for (m in setdiff(akeR_methods()$method, c("unifrac", "wunifrac"))) {
    d <- calc_dist(sim$counts, method = m)
    expect_s3_class(d, "dist")
    expect_equal(attr(d, "Labels"), colnames(sim$counts))
    expect_equal(attr(d, "method"), m)
    expect_true(all(d >= 0))
  }
})

test_that("aitchison matches a manual CLR + euclidean", {
  sim <- sim_data()
  x <- t(sim$counts) + 1
  clr <- log(x) - rowMeans(log(x))
  expect_equal(as.vector(calc_dist(sim$counts, "aitchison")),
    as.vector(dist(clr)))
})

test_that("bray uses relative abundances unless relative = FALSE", {
  sim <- sim_data()
  rel <- t(sim$counts) / colSums(sim$counts)
  expect_equal(as.vector(calc_dist(sim$counts, "bray")),
    as.vector(vegan::vegdist(rel, "bray")))
  expect_equal(as.vector(calc_dist(sim$counts, "bray", relative = FALSE)),
    as.vector(vegan::vegdist(t(sim$counts), "bray")))
})

test_that("jsd is 0 for identical and 1 for disjoint profiles", {
  m <- cbind(a = c(1, 1, 0, 0), b = c(2, 2, 0, 0), c = c(0, 0, 3, 1))
  d <- as.matrix(calc_dist(m, "jsd"))
  expect_equal(d["a", "b"], 0)
  expect_equal(d["a", "c"], 1)
})

test_that("taxa_are_rows = FALSE and data.frame input are accepted", {
  sim <- sim_data()
  d1 <- calc_dist(sim$counts, "bray")
  d2 <- calc_dist(as.data.frame(t(sim$counts)), "bray", taxa_are_rows = FALSE)
  expect_equal(as.vector(d1), as.vector(d2))
})

test_that("bad input gives informative errors", {
  sim <- sim_data()
  expect_error(calc_dist(sim$counts, "nonsense"))
  bad <- sim$counts
  bad[1, 1] <- -1
  expect_error(calc_dist(bad), "non-negative")
  empty <- sim$counts
  empty[, 1] <- 0
  expect_error(calc_dist(empty), "no counts")
  expect_error(calc_dist(sim$counts, "unifrac"), "phylogenetic tree")
})

test_that("phyloseq input works", {
  skip_if_not_installed("phyloseq")
  sim <- sim_data()
  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(sim$counts, taxa_are_rows = TRUE),
    phyloseq::sample_data(sim$metadata)
  )
  expect_equal(as.vector(calc_dist(ps, "bray")),
    as.vector(calc_dist(sim$counts, "bray")))
})

test_that("SummarizedExperiment input works", {
  skip_if_not_installed("SummarizedExperiment")
  sim <- sim_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = sim$counts),
    colData = sim$metadata
  )
  expect_equal(as.vector(calc_dist(se, "aitchison")),
    as.vector(calc_dist(sim$counts, "aitchison")))
})

test_that("robust aitchison treats each sample on its own", {
  sim <- sim_data()
  sim$counts[sample(length(sim$counts), 60)] <- 0
  full <- as.matrix(calc_dist(sim$counts, "robust.aitchison"))
  pair <- calc_dist(sim$counts[, c("s1", "s2")], "robust.aitchison")
  expect_equal(full["s1", "s2"], as.vector(pair))
  x <- sim$counts[, "s1"]
  y <- sim$counts[, "s2"]
  rclr <- function(v) ifelse(v > 0, log(v) - mean(log(v[v > 0])), 0)
  expect_equal(as.vector(pair), sqrt(sum((rclr(x) - rclr(y))^2)))
})
