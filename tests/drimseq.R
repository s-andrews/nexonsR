library(nexonsR)

local({
  prepare <- getFromNamespace("prepare_drimseq", "nexonsR")
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- paste0("sample-", 1:6)
  annotation <- data.frame(Transcript_ID = paste0("t", 1:6),
    Gene_ID = rep(c("g2", "g1"), each = 3), Gene_Name = "gene",
    Chr = "1", Start = 1, End = 10, Strand = "+")
  counts <- matrix(11:46, nrow = 6, dimnames = list(NULL, samples))
  original <- cbind(annotation, counts)
  write_table <- function(data, file = "unique") utils::write.table(data,
    file.path(folder, paste0(file, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write_table(original[, c(1:7, 13:8)])
  write_table(original, "partial")
  object <- suppressWarnings(read_nexons(folder, ""))
  metadata <- object$metadata
  metadata$condition <- factor(rep(c("control", "treated"), each = 3),
                               levels = c("control", "treated", "unused"))
  metadata$batch <- factor(rep(c("A", "B", "A"), 2))
  object$set_metadata(metadata)
  selected <- samples[c(6, 2, 4, 1)]
  prep <- function(...) prepare(object, min_samps_feature_expr = 2, ...)
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  stopifnot(!length(object$loaded))
  for (file in c("partial", "unique")) {
    input <- prep(design = ~ condition, file = file, samples = selected)
    if (file == "partial") stopifnot(identical(object$loaded, "partial"))
    expected <- counts[, selected]
    storage.mode(expected) <- "integer"
    stopifnot(identical(as.matrix(input$counts[, -(1:2)]), expected),
      identical(input$counts$gene_id, annotation$Gene_ID),
      identical(input$counts$feature_id, annotation$Transcript_ID),
      identical(input$metadata$sample, selected),
      identical(input$metadata$sample_id, selected),
      identical(levels(input$metadata$condition), c("control", "treated")),
      identical(unname(input$design[, 2]), c(1, 0, 1, 0)),
      input$filter$min_samps_gene_expr == 4,
      input$filter$min_samps_feature_prop == 2)
  }
  stopifnot(identical(object$metadata, metadata))
  stopifnot(ncol(prep(design = ~1)$design) == 1L)
  matrix_design <- stats::model.matrix(~ batch + condition, droplevels(metadata))
  stopifnot(identical(prep(design = ~ batch + condition)$design,
                     prep(design = matrix_design)$design))
  continuous <- metadata
  continuous$age <- c(20, 30, 45, 25, 35, 50)
  object$set_metadata(continuous)
  stopifnot(ncol(prep(design = ~age * condition)$design) == 4L)
  object$set_metadata(metadata)
  error(object$run_drimseq(), "design")
  error(object$run_drimseq(~condition), "min_samps_feature_expr")
  for (file in list("gene", "other", NA_character_, NULL, c("partial", "unique"))) {
    error(prep(design = ~condition, file = file), "`file`")
  }
  for (selection in list(character(), "unknown", c(samples[1], samples[1]), NA_character_)) {
    error(prep(design = ~condition, samples = selection), "`samples`")
  }
  for (design in list("condition", matrix(1, 5, 1), matrix(NA_real_, 6, 1),
                     matrix(Inf, 6, 1), matrix(numeric(), 6, 0))) {
    error(prep(design = design), "`design`")
  }
  error(prep(design = matrix(1, 6, 2)), "rank")
  error(prep(design = diag(6)), "residual")
  error(prep(design = condition ~ batch), "one-sided")
  error(prep(design = ~unknown_column), "unknown_column")
  bad_metadata <- metadata
  bad_metadata$batch[1] <- NA
  object$set_metadata(bad_metadata)
  error(prep(design = ~batch + condition), "missing")
  object$set_metadata(metadata)
  for (name in c("min_samps_feature_expr", "min_samps_feature_prop", "min_samps_gene_expr",
                 "min_feature_expr", "min_feature_prop", "min_gene_expr")) {
    for (value in list(-1, NA_real_, Inf, "3", numeric(), c(1, 2))) {
      args <- list(object = object, design = ~condition, min_samps_feature_expr = 2)
      args[[name]] <- value
      error(do.call(prepare, args), name)
    }
  }
  for (name in c("min_samps_feature_expr", "min_samps_feature_prop", "min_samps_gene_expr")) {
    for (value in c(1.5, 7)) {
      args <- list(object = object, design = ~condition, min_samps_feature_expr = 2)
      args[[name]] <- value
      error(do.call(prepare, args), name)
    }
  }
  error(prep(design = ~condition, min_feature_prop = 1.1), "min_feature_prop")
  for (value in list(NA, 1, NULL, c(TRUE, FALSE))) {
    error(prep(design = ~condition, run_gene_twice = value), "run_gene_twice")
  }
  for (column in c("Gene_ID", "Transcript_ID")) {
    for (value in c(NA_character_, "", " ")) {
      invalid <- original
      invalid[1, column] <- value
      write_table(invalid)
      object$clear_cache()
      error(prep(design = ~condition), column)
    }
  }
  invalid <- original
  invalid$Transcript_ID[2] <- invalid$Transcript_ID[1]
  write_table(invalid)
  object$clear_cache()
  error(prep(design = ~condition), "unique identifiers")
  for (value in c(-1, 0.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    invalid <- original
    invalid[1, 8] <- value
    write_table(invalid)
    object$clear_cache()
    error(prep(design = ~condition),
          if (is.finite(value) && value > .Machine$integer.max) "integer range" else "whole-number")
  }
  write_table(original)
  object$clear_cache()
  zero_filter <- prepare(object, ~condition, min_samps_feature_expr = 0,
    min_samps_feature_prop = 0, min_samps_gene_expr = 0,
    min_feature_expr = 0, min_feature_prop = 0, min_gene_expr = 0)
  stopifnot(all(unlist(zero_filter$filter) == 0))
  write_table(original[FALSE, ])
  object$clear_cache()
  error(prep(design = ~condition), "no rows")
  write_table(original)
  object$clear_cache()
  # Check the proportion filter separately from the expression filter.
  filter_counts <- data.frame(gene_id = rep(c("keep", "drop"), each = 3),
    feature_id = paste0("f", 1:6), matrix(rep(c(50, 50, 1, 100, 1, 1), 6), ncol = 6))
  names(filter_counts)[-(1:2)] <- samples
  filtered <- DRIMSeq::dmFilter(DRIMSeq::dmDSdata(filter_counts,
    data.frame(sample_id = samples)), min_samps_feature_prop = 3,
    min_feature_prop = 0.1)
  stopifnot(setequal(as.character(DRIMSeq::counts(filtered)$feature_id), c("f1", "f2")))
  # Many replicated genes support DRIMSeq's default precision estimation.
  set.seed(12)
  counts <- matrix(stats::rnbinom(1200, mu = 100, size = 20), ncol = 6,
                   dimnames = list(NULL, samples))
  annotation <- data.frame(Transcript_ID = paste0("t", 1:200),
    Gene_ID = rep(paste0("g", 1:100), each = 2), Gene_Name = "gene",
    Chr = "1", Start = 1, End = 10, Strand = "+")
  counts[1, ] <- 1 # Removing this weak feature also removes its singleton gene.
  data <- cbind(annotation, counts)
  write_table(data)
  write_table(data, "partial")
  object$clear_cache()
  fits <- list()
  for (file in c("unique", "partial")) {
    set.seed(123)
    fits[[file]] <- object$run_drimseq(
      if (file == "unique") ~batch + condition else matrix_design,
      file = file, min_samps_feature_expr = 3)
    stopifnot(inherits(fits[[file]], "dmDSfit"),
      !"g1" %in% DRIMSeq::counts(fits[[file]])$gene_id,
      identical(as.character(DRIMSeq::samples(fits[[file]])$sample_id), samples))
  }
  stopifnot(isTRUE(all.equal(DRIMSeq::proportions(fits$unique),
                            DRIMSeq::proportions(fits$partial))))
  tested <- DRIMSeq::dmTest(fits$unique, coef = "conditiontreated")
  stopifnot(nrow(DRIMSeq::results(tested)) == 99L,
    nrow(DRIMSeq::results(tested, level = "feature")) == 198L,
    identical(object$metadata, metadata))
  error(object$run_drimseq(~condition, min_samps_feature_expr = 3,
                          min_gene_expr = 1e9), "no analysable genes")
  unlink(file.path(folder, "unique.txt"))
  missing_file <- suppressWarnings(read_nexons(folder, ""))
  error(missing_file$run_drimseq(~1, min_samps_feature_expr = 2), "No unique output")
})
