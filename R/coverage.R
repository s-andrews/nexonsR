# Coverage allows fractional counts and does not require DESeq2.
nexons_coverage <- function(object, file, level, threshold) {
  validate_string(file, "file")
  validate_string(level, "level")
  if (!file %in% c("gene", "partial", "unique")) {
    stop("`file` must be 'gene', 'partial', or 'unique'.", call. = FALSE)
  }
  if (!level %in% c("gene", "isoform")) {
    stop("`level` must be 'gene' or 'isoform'.", call. = FALSE)
  }
  if (file == "gene" && level == "isoform") {
    stop("Isoform coverage requires a 'partial' or 'unique' file.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0) {
    stop("`threshold` must be a single finite positive number.", call. = FALSE)
  }
  data <- object$get_data(file)
  metadata <- object$get_metadata()
  samples <- metadata$sample
  n_annotation <- if (file == "gene") 6L else 7L
  columns <- seq.int(n_annotation + 1L, ncol(data))
  positions <- match(samples, names(data)[columns])
  if (anyNA(positions)) stop("Count table samples do not match metadata.", call. = FALSE)
  counts <- as.matrix(data[, columns[positions], drop = FALSE])
  if (length(counts) &&
      (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0))) {
    stop("Coverage requires finite, non-negative counts.", call. = FALSE)
  }
  id_column <- if (level == "gene") "Gene_ID" else "Transcript_ID"
  ids <- data[[id_column]]
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(id_column, " must not contain missing or empty identifiers.", call. = FALSE)
  }
  if (level == "gene" && file != "gene") {
    if (nrow(counts)) counts <- rowsum(counts, ids, reorder = FALSE)
  } else if (anyDuplicated(ids)) {
    stop(id_column, " must contain unique identifiers for this analysis.", call. = FALSE)
  }
  result <- metadata
  count_name <- make.unique(c(names(result), level))[ncol(result) + 1L]
  result[[count_name]] <- unname(colSums(counts >= threshold))
  result
}
