# Internal workers for the two akeR() modes. Both return a data.frame of
# distances with the samples' metadata attached.

.reference_mode <- function(dm, metadata, column) {
  samples <- rownames(dm)
  refs <- as.character(metadata[[column]])
  refs[refs == ""] <- NA

  missing <- setdiff(stats::na.omit(refs), samples)
  if (length(missing) > 0) {
    stop(length(missing), " reference ID(s) in column '", column,
      "' are not samples in the count table, e.g. '", missing[1], "'.",
      call. = FALSE)
  }
  has_ref <- !is.na(refs)
  if (!any(has_ref)) {
    stop("Column '", column, "' has no reference IDs.", call. = FALSE)
  }

  # Own reference, plus every other reference sample
  all_refs <- unique(refs[has_ref])
  query <- samples[has_ref]
  own <- refs[has_ref]
  other <- lapply(own, function(r) setdiff(all_refs, r))
  out <- data.frame(
    sample = c(query, rep(query, lengths(other))),
    reference = c(own, unlist(other)),
    comparison = factor(c(rep("own", length(query)),
      rep("other", sum(lengths(other)))), levels = c("own", "other")),
    stringsAsFactors = FALSE
  )
  out$distance <- dm[cbind(out$sample, out$reference)]
  out <- out[order(match(out$sample, samples), out$comparison),
    c("sample", "reference", "distance", "comparison")]
  extra <- setdiff(names(metadata), names(out))
  out <- cbind(out, metadata[out$sample, extra, drop = FALSE])
  droplevels(out)
}

.group_mode <- function(d, metadata, column) {
  # Drop samples without a group label
  keep <- !is.na(metadata[[column]])
  if (!all(keep)) {
    d <- stats::as.dist(as.matrix(d)[keep, keep])
    metadata <- metadata[keep, , drop = FALSE]
  }
  labels <- as.character(metadata[[column]])

  samples <- attr(d, "Labels")
  idx <- which(lower.tri(matrix(NA, length(samples), length(samples))),
    arr.ind = TRUE)
  i <- idx[, "col"]
  j <- idx[, "row"]
  within <- labels[i] == labels[j]

  out <- data.frame(
    sample_1 = samples[i],
    sample_2 = samples[j],
    distance = as.vector(d),
    comparison = factor(ifelse(within, "within", "between"),
      levels = c("within", "between")),
    group_pair = ifelse(within, labels[i],
      paste(pmin(labels[i], labels[j]), "vs", pmax(labels[i], labels[j]))),
    stringsAsFactors = FALSE
  )
  for (col in names(metadata)) {
    out[[paste0(col, "_1")]] <- metadata[[col]][i]
    out[[paste0(col, "_2")]] <- metadata[[col]][j]
  }
  out
}
