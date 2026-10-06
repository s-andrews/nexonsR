library(nexonsR)

local({
  prepare <- getFromNamespace("prepare_deseq2", "nexonsR")
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- paste0("sample-", 1:6)
  annotation <- data.frame(Gene_ID = c("g2", "g1", "g2"),
    Gene_Name = "gene", Chr = "1", Start = 1, End = 10, Strand = "+")
  counts <- matrix(seq_len(18), nrow = 3, dimnames = list(NULL, samples))
  isoforms <- cbind(Transcript_ID = c("t1", "t2", "t3"), annotation, counts)
  write_table <- function(data, file) utils::write.table(data,
    file.path(folder, paste0(file, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write_table(isoforms, "partial")
  write_table(isoforms[, c(1:7, 13:8)], "unique")
  genes <- cbind(annotation[1:2, ], counts[1:2, ])
  write_table(genes, "gene")
  object <- read_nexons(folder, "")
  metadata <- object$metadata
  metadata$condition <- factor(rep(c("control", "treated"), each = 3))
  object$set_metadata(metadata)
  selected <- samples[c(6, 2, 4, 1)]
  for (file in c("gene", "partial", "unique")) {
    input <- prepare(object, file, "gene", selected)
    expected <- if (file == "gene") counts[1:2, selected] else
      rbind(counts[1, selected] + counts[3, selected], counts[2, selected])
    rownames(expected) <- c("g2", "g1")
    stopifnot(identical(input$counts, expected),
      identical(rownames(input$metadata), selected),
      identical(input$metadata$sample, selected),
      identical(as.character(input$metadata$condition), c("treated", "control", "treated", "control")))
  }
  for (file in c("partial", "unique")) {
    input <- prepare(object, file, "isoform", NULL)
    expected <- counts
    rownames(expected) <- c("t1", "t2", "t3")
    stopifnot(identical(input$counts, expected))
  }
  stopifnot(ncol(prepare(object, "gene", "gene", samples[1])$counts) == 1L,
            identical(object$metadata, metadata))
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  error(object$run_deseq2(), "design")
  for (file in c("g", "other")) error(object$run_deseq2(~condition, file), "`file`")
  error(object$run_deseq2(~condition, level = "other"), "`level`")
  error(object$run_deseq2(~condition, level = "isoform"), "requires")
  for (selection in list(character(), "unknown", c(samples[1], samples[1]), NA_character_)) {
    error(object$run_deseq2(~condition, samples = selection), "`samples`")
  }
  for (value in c(-1, 0.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    invalid <- genes
    invalid[1, 7] <- value
    write_table(invalid, "gene")
    object$clear_cache()
    error(prepare(object, "gene", "gene", NULL),
          if (is.finite(value) && value > .Machine$integer.max) "integer range" else "whole-number")
  }
  write_table(genes, "gene")
  unlink(file.path(folder, "unique.txt"))
  missing_file <- suppressWarnings(read_nexons(folder, ""))
  error(missing_file$run_deseq2(~1, file = "unique"), "No unique output")

  # Exercise real model fitting on replicated, overdispersed counts.
  set.seed(123)
  counts <- matrix(stats::rnbinom(6000, mu = 100, size = 5), ncol = 6,
                   dimnames = list(NULL, samples))
  annotation <- data.frame(Gene_ID = paste0("g", rep(1:500, each = 2)),
    Gene_Name = "gene", Chr = "1", Start = 1, End = 10, Strand = "+")
  write_table(cbind(Transcript_ID = paste0("t", 1:1000), annotation, counts), "partial")
  annotation$Gene_ID <- paste0("g", 1:1000)
  write_table(cbind(annotation, counts), "gene")
  object <- suppressWarnings(read_nexons(folder, ""))
  metadata <- object$metadata
  metadata$condition <- factor(rep(c("control", "treated"), each = 3))
  object$set_metadata(metadata)
  for (level in c("gene", "isoform")) {
    model <- object$run_deseq2(~condition, file = "partial", level = level,
                             samples = selected)
    stopifnot(inherits(model, "DESeqDataSet"),
      identical(colnames(model), selected),
      nrow(DESeq2::results(model)) == if (level == "gene") 500L else 1000L)
  }
  model <- object$run_deseq2(~condition)
  stopifnot(inherits(model, "DESeqDataSet"), length(DESeq2::sizeFactors(model)) == 6L)
})
