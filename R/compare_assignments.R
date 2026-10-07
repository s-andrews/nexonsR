#' Compare quantitations across assignment files
#'
#' @name compare_assignments
#' @param level Either `"gene"` (default) or `"isoform"`.
#' @param samples Optional non-empty character vector of unique sample names.
#'   `NULL` selects all samples in metadata order.
#' @param units Either `"counts"` (default) or `"log2RPM"`.
#' @return A data frame with one row per feature and sample. Gene output has
#'   `Gene_ID`, `Number_of_Isoforms`, `sample`, and available `gene`, `partial`,
#'   `unique` columns. Isoform output has `Transcript_ID`, `Gene_ID`, `sample`,
#'   `partial`, and `unique` columns.
#' @details Gene quantitations sum raw isoform counts by gene before log2RPM
#'   conversion. Each file uses its own sample totals and the formula
#'   `log2(count / sample_total * 1e6 + 1)`. Normalization requires finite,
#'   non-negative counts and positive sample totals, as in [get_data()].
#'   Isoform counts include zero-count transcripts and use the larger distinct
#'   transcript count from partial and unique, or the sole available count.
#'   They are `NA` if neither isoform file contains the gene.
#'
#'   Gene output omits columns for unavailable files. Isoform output requires
#'   both partial and unique. Features are aligned by gene ID or the
#'   transcript/gene pair. All features are retained; missing quantitations
#'   are `NA`. Differing feature sets produce one warning per call.
#'   Missing/empty identifiers and duplicate feature keys within a file raise
#'   errors. Rows form sample blocks in selected sample order, with features
#'   in first-occurrence order across gene, partial, then unique files.
#'   Empty tables return typed zero-row output for counts.
#'   Raw tables are cached through [get_data()] and remain unchanged.
#' @rawRd \usage{\special{x$compare_assignments(level = "gene", samples = NULL,
#'   units = "counts")}}
#' @seealso [read_nexons()], [get_data()]
#' @examples
#' \dontrun{
#' x <- read_nexons("results")
#' x$compare_assignments()
#' x$compare_assignments("isoform", samples = "sample-1", units = "log2RPM")
#' }
NULL

nexons_compare_assignments <- function(object, level, samples, units) {
  validate_string(level, "level")
  validate_string(units, "units")
  if (!level %in% c("gene", "isoform")) {
    stop("`level` must be 'gene' or 'isoform'.", call. = FALSE)
  }
  if (!units %in% c("counts", "log2RPM")) {
    stop("`units` must be 'counts' or 'log2RPM'.", call. = FALSE)
  }
  if (is.null(samples)) samples <- object$get_metadata()$sample
  types <- c("gene", "partial", "unique")
  if (level == "isoform") {
    types <- c("partial", "unique")
    if (anyNA(object$files[types])) {
      stop("Isoform comparison requires both partial and unique output files.",
           call. = FALSE)
    }
  } else {
    types <- types[!is.na(object$files[types])]
  }
  tables <- lapply(types, function(type) {
    data <- object$get_data(type, samples = samples, units = "counts")
    id_columns <- if (type == "gene") "Gene_ID" else c("Transcript_ID", "Gene_ID")
    for (column in id_columns) {
      ids <- data[[column]]
      if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
        stop(column, " must not contain missing or empty identifiers.", call. = FALSE)
      }
    }
    # Length-prefix the transcript ID so arbitrary ID characters cannot
    # create collisions between distinct transcript/gene pairs.
    keys <- if (type == "gene") data$Gene_ID else if (!nrow(data)) character() else
      paste0(nchar(data$Transcript_ID), ":", data$Transcript_ID, data$Gene_ID)
    if (anyDuplicated(keys)) {
      stop("Duplicate feature keys in ", type, " output file.", call. = FALSE)
    }
    n_annotation <- if (type == "gene") 6L else 7L
    counts <- as.matrix(data[, n_annotation + seq_along(samples), drop = FALSE])
    if (units == "log2RPM" &&
        (any(!is.finite(counts)) || any(counts < 0))) {
      stop("log2RPM requires finite, non-negative counts.", call. = FALSE)
    }
    isoforms <- NULL
    if (level == "gene") {
      genes <- unique(data$Gene_ID)
      if (type != "gene") {
        isoforms <- tabulate(match(data$Gene_ID, genes), nbins = length(genes))
        if (nrow(counts)) counts <- rowsum(counts, data$Gene_ID, reorder = FALSE)
      }
      keys <- genes
      features <- data.frame(Gene_ID = genes, stringsAsFactors = FALSE)
    } else {
      features <- data[, c("Transcript_ID", "Gene_ID"), drop = FALSE]
    }
    if (units == "log2RPM") {
      totals <- colSums(counts)
      if (any(!is.finite(totals)) || any(totals <= 0)) {
        stop("log2RPM requires finite, positive sample totals; zero total counts are undefined.",
             call. = FALSE)
      }
      counts <- log2(sweep(counts, 2L, totals, "/") * 1e6 + 1)
    }
    list(keys = keys, features = features, counts = counts, isoforms = isoforms)
  })
  names(tables) <- types
  keys <- unique(unlist(lapply(tables, `[[`, "keys"), use.names = FALSE))
  if (any(vapply(tables, function(table) !setequal(table$keys, keys), logical(1)))) {
    warning("Feature sets differ across assignment files; missing quantitations are NA.",
            call. = FALSE)
  }
  features <- do.call(rbind, lapply(tables, `[[`, "features"))
  all_keys <- unlist(lapply(tables, `[[`, "keys"), use.names = FALSE)
  features <- features[match(keys, all_keys), , drop = FALSE]
  if (level == "gene") {
    number <- rep(NA_integer_, length(keys))
    for (table in tables) {
      if (is.null(table$isoforms)) next
      positions <- match(table$keys, keys)
      previous <- number[positions]
      previous[is.na(previous)] <- 0L
      number[positions] <- pmax(previous, table$isoforms)
    }
    features$Number_of_Isoforms <- number
  }
  result <- features[rep(seq_along(keys), times = length(samples)), , drop = FALSE]
  result$sample <- rep(samples, each = length(keys))
  for (type in types) {
    table <- tables[[type]]
    aligned <- matrix(NA_real_, nrow = length(keys), ncol = length(samples))
    if (length(table$keys)) aligned[match(table$keys, keys), ] <- table$counts
    result[[type]] <- as.vector(aligned)
  }
  rownames(result) <- NULL
  result
}
