library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- c("sample 1", "Gene_ID")
  make_table <- function(genes, counts, transcripts = NULL) {
    annotation <- data.frame(Gene_ID = genes, Gene_Name = genes,
      Chr = rep("1", length(genes)), Start = rep("1", length(genes)),
      End = rep("10", length(genes)), Strand = rep("+", length(genes)))
    if (!is.null(transcripts)) annotation <- cbind(Transcript_ID = transcripts, annotation)
    result <- cbind(annotation, counts)
    names(result) <- c(names(annotation), samples)
    result
  }
  write_table <- function(data, type) utils::write.table(data,
    file.path(folder, paste0(type, ".txt")), sep = "\t", quote = FALSE,
    row.names = FALSE)
  object <- function() suppressWarnings(read_nexons(folder, ""))
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
      grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  warned <- function(expr) {
    warnings <- character()
    result <- withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    stopifnot(length(warnings) == 1L,
      grepl("Feature sets differ", warnings, fixed = TRUE))
    result
  }
  gene <- make_table(c("001", "002"), matrix(c(10, 30, 20, 20), 2))
  partial <- make_table(c("001", "001", "002"),
    matrix(c(1.5, 2.5, 6, 0, 2, 8), 3), c("01", "02", "03"))
  unique <- make_table(c("002", "001", "001"),
    matrix(c(3, 1, 0, 4, 0, 1), 3), c("03", "01", "02"))
  write_table(gene, "gene")
  write_table(partial, "partial")
  # Different source sample order must not change the comparison alignment.
  reordered <- unique[, c(1:7, 9, 8)]
  names(reordered) <- names(unique)[c(1:7, 9, 8)]
  write_table(reordered, "unique")
  x <- object()
  raw <- x$get_data("partial")
  result <- x$compare_assignments()
  stopifnot(identical(names(result),
    c("Gene_ID", "Number_of_Isoforms", "sample", "gene", "partial", "unique")),
    identical(result$Gene_ID, rep(c("001", "002"), 2)),
    identical(result$sample, rep(samples, each = 2)),
    identical(result$Number_of_Isoforms, rep(c(2L, 1L), 2)),
    identical(result$gene, c(10, 30, 20, 20)),
    identical(result$partial, c(4, 6, 2, 8)),
    identical(result$unique, c(1, 3, 1, 4)))
  normalized <- x$compare_assignments(samples = rev(samples), units = "log2RPM")
  stopifnot(isTRUE(all.equal(normalized$partial,
    log2(c(2, 8, 4, 6) / 10 * 1e6 + 1))),
    isTRUE(all.equal(normalized$unique,
    log2(c(1 / 5, 4 / 5, 1 / 4, 3 / 4) * 1e6 + 1))),
    identical(x$get_data("partial"), raw))
  iso <- x$compare_assignments("isoform", samples = samples[1])
  stopifnot(identical(names(iso), c("Transcript_ID", "Gene_ID", "sample", "partial", "unique")),
    identical(iso$Transcript_ID, c("01", "02", "03")),
    identical(iso$partial, c(1.5, 2.5, 6)), identical(iso$unique, c(1, 0, 3)))
  iso_log <- x$compare_assignments("isoform", samples = samples[1], units = "log2RPM")
  stopifnot(isTRUE(all.equal(iso_log$partial, log2(c(1.5, 2.5, 6) / 10 * 1e6 + 1))))
  for (selection in list(character(), "unknown", rep(samples[1], 2), NA_character_, 1)) {
    error(x$compare_assignments(samples = selection), "`samples`")
  }
  error(x$compare_assignments(level = "transcript"), "`level`")
  error(x$compare_assignments(units = "RPM"), "`units`")

  # Different isoform counts for the same gene use the maximum without a
  # gene-level feature warning; isoform-level mismatches warn once.
  write_table(unique[-3, ], "unique")
  x <- object()
  stopifnot(identical(x$compare_assignments()$Number_of_Isoforms, rep(c(2L, 1L), 2)))
  mismatch <- warned(x$compare_assignments("isoform"))
  stopifnot(all(is.na(mismatch$unique[mismatch$Transcript_ID == "02"])))
  write_table(gene[1, ], "gene")
  mismatch <- warned(object()$compare_assignments())
  stopifnot(all(is.na(mismatch$gene[mismatch$Gene_ID == "002"])))

  # Every nonempty subset of available files has the corresponding columns.
  for (mask in 1:7) {
    types <- c("gene", "partial", "unique")
    present <- as.logical(intToBits(mask)[1:3])
    unlink(file.path(folder, paste0(types, ".txt")))
    fixtures <- list(gene, partial, unique)
    for (i in which(present)) write_table(fixtures[[i]], types[i])
    x <- object()
    result <- x$compare_assignments()
    stopifnot(identical(names(result),
      c("Gene_ID", "Number_of_Isoforms", "sample", types[present])))
    if (!any(present[2:3])) stopifnot(all(is.na(result$Number_of_Isoforms)))
    if (!all(present[2:3])) error(x$compare_assignments("isoform"), "requires both")
  }
  x <- object()
  x$compare_assignments("isoform")
  stopifnot(identical(x$loaded, c("partial", "unique")))

  write_table(partial[c(1, 1, 3), ], "partial")
  error(object()$compare_assignments(), "Duplicate feature keys")
  bad <- partial
  bad[1, 2] <- ""
  write_table(bad, "partial")
  error(object()$compare_assignments(), "identifiers")
  for (value in c(-1, Inf, NA_real_)) {
    bad <- partial
    bad[1, 8] <- value
    write_table(bad, "partial")
    error(object()$compare_assignments(units = "log2RPM"), "finite, non-negative")
  }
  bad <- partial
  bad[, 8] <- 0
  write_table(bad, "partial")
  error(object()$compare_assignments(units = "log2RPM"), "zero total counts")
  write_table(gene[FALSE, ], "gene")
  write_table(partial[FALSE, ], "partial")
  write_table(unique[FALSE, ], "unique")
  x <- object()
  result <- x$compare_assignments()
  stopifnot(nrow(result) == 0L, is.character(result$Gene_ID),
    is.integer(result$Number_of_Isoforms), is.character(result$sample),
    is.numeric(result$gene), nrow(x$compare_assignments("isoform")) == 0L)
  error(x$compare_assignments(units = "log2RPM"), "zero total counts")
})
