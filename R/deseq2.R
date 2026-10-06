# Assemble inputs separately from fitting so alignment and aggregation can be
# checked without requiring the optional Bioconductor dependency.
prepare_deseq2 <- function(object, file, level, samples) {
  prepare_model_counts(object, file, level, samples, "DESeq2")
}

prepare_model_counts <- function(object, file, level, samples, package) {
  validate_string(file, "file")
  validate_string(level, "level")
  if (!file %in% c("gene", "partial", "unique")) {
    stop("`file` must be 'gene', 'partial', or 'unique'.", call. = FALSE)
  }
  if (!level %in% c("gene", "isoform")) {
    stop("`level` must be 'gene' or 'isoform'.", call. = FALSE)
  }
  if (file == "gene" && level == "isoform") {
    stop("Isoform analysis requires a 'partial' or 'unique' file.", call. = FALSE)
  }
  if (is.na(object$files[[file]])) {
    stop("No ", file, " output file was found for this object.", call. = FALSE)
  }
  metadata <- object$get_metadata()
  if (is.null(samples)) samples <- metadata$sample
  if (!is.character(samples) || !length(samples) || anyNA(samples) ||
      anyDuplicated(samples) || !all(samples %in% metadata$sample)) {
    stop("`samples` must be a non-empty character vector of unique known sample names.",
         call. = FALSE)
  }
  metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]
  metadata <- droplevels(metadata)
  rownames(metadata) <- samples
  data <- object$get_data(file)
  # Select only count columns by position, avoiding collisions with annotations.
  n_annotation <- if (file == "gene") 6L else 7L
  count_columns <- seq.int(n_annotation + 1L, ncol(data))
  positions <- match(samples, names(data)[count_columns])
  if (anyNA(positions)) stop("Count table samples do not match metadata.", call. = FALSE)
  counts <- as.matrix(data[, count_columns[positions], drop = FALSE])
  colnames(counts) <- samples
  if (!nrow(counts)) stop("The selected count table has no rows.", call. = FALSE)
  if (!is.numeric(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts))) {
    stop(package, " requires finite, non-negative whole-number counts; counts are not rounded.",
         call. = FALSE)
  }
  id_column <- if (level == "gene") "Gene_ID" else "Transcript_ID"
  ids <- data[[id_column]]
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(id_column, " must not contain missing or empty identifiers.", call. = FALSE)
  }
  if (level == "gene" && file != "gene") {
    counts <- rowsum(counts, group = ids, reorder = FALSE)
  } else {
    if (anyDuplicated(ids)) {
      stop(id_column, " must contain unique identifiers for this analysis.", call. = FALSE)
    }
    rownames(counts) <- ids
  }
  if (any(counts > .Machine$integer.max)) {
    stop("Counts exceed the integer range supported by ", package, ".", call. = FALSE)
  }
  storage.mode(counts) <- "integer"
  list(counts = counts, metadata = metadata)
}
