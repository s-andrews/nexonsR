library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  annotation <- data.frame(Gene_ID = c("g1", "g1", "g2", "g3"),
    Gene_Name = "gene", Chr = "1", Start = 1, End = 10, Strand = "+")
  counts <- data.frame(c(0.5, 0.5, 0, 2), c(0, 0, 1, 0), check.names = FALSE)
  names(counts) <- c("sample-B", "Gene_ID")
  isoforms <- cbind(Transcript_ID = paste0("t", 1:4), annotation, counts)
  genes <- cbind(annotation[c(1, 3, 4), ],
                 data.frame(c(1, 0, 2), c(0, 1, 0)))
  names(genes)[7:8] <- names(counts)
  write_table <- function(data, file) utils::write.table(data,
    file.path(folder, paste0(file, ".txt")), sep = "\t", quote = FALSE,
    row.names = FALSE)
  write_table(genes, "gene")
  write_table(isoforms, "partial")
  reordered <- isoforms[, c(1:7, 9, 8)]
  names(reordered) <- names(isoforms)[c(1:7, 9, 8)]
  write_table(reordered, "unique")
  x <- read_nexons(folder, "")
  expected <- data.frame(sample = names(counts), gene = c(2, 1))
  stopifnot(length(x$loaded) == 0L, identical(x$coverage(), expected),
            identical(x$loaded, "gene"))
  for (file in c("gene", "partial", "unique")) {
    stopifnot(identical(x$coverage(file), expected),
      identical(x$coverage(file, threshold = 2)$gene, c(1, 0)),
      identical(x$coverage(file, threshold = 3)$gene, c(0, 0)))
  }
  for (file in c("partial", "unique")) {
    stopifnot(identical(x$coverage(file, "isoform"),
      data.frame(sample = names(counts), isoform = c(1, 1))),
      identical(x$coverage(file, "isoform", 0.5)$isoform, c(3, 1)))
  }
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  error(x$coverage(level = "isoform"), "requires")
  for (value in list("other", NA_character_, NULL, c("gene", "partial"))) {
    error(x$coverage(file = value), "`file`")
    error(x$coverage(level = value), "`level`")
  }
  for (value in list(0, -1, NA_real_, Inf, NaN, "1", TRUE, NULL, c(1, 2))) {
    error(x$coverage(threshold = value), "`threshold`")
  }
  # Metadata is retained and aligned even when source sample order differs.
  original_metadata <- x$get_metadata()
  metadata <- original_metadata[2:1, , drop = FALSE]
  metadata$condition <- factor(c("treated", "control"))
  metadata$batch <- c(2L, 1L)
  x$set_metadata(metadata)
  saved_metadata <- x$get_metadata()
  for (file in c("gene", "partial", "unique")) {
    annotated <- saved_metadata
    annotated$gene <- c(2, 1)
    stopifnot(identical(x$coverage(file), annotated))
  }
  annotated <- saved_metadata
  annotated$isoform <- c(1, 1)
  stopifnot(identical(x$coverage("unique", "isoform"), annotated),
            identical(x$get_metadata(), saved_metadata))
  # A metadata column with the same name must not be overwritten.
  metadata <- saved_metadata
  metadata$gene <- c("a", "b")
  metadata$gene.1 <- c("c", "d")
  x$set_metadata(metadata)
  annotated <- metadata
  annotated$gene.2 <- c(2, 1)
  stopifnot(identical(x$coverage(), annotated),
            identical(x$get_metadata(), metadata))
  x$set_metadata(original_metadata)
  # Cached counts survive source changes until explicitly cleared.
  write_table(genes[FALSE, ], "gene")
  stopifnot(identical(x$coverage(), expected))
  x$clear_cache()
  stopifnot(identical(x$coverage()$gene, c(0, 0)))
  for (file in c("partial", "unique")) {
    empty <- isoforms[FALSE, ]
    names(empty) <- names(isoforms)
    write_table(empty, file)
    stopifnot(identical(x$coverage(file)$gene, c(0, 0)),
              identical(x$coverage(file, "isoform")$isoform, c(0, 0)))
  }
  for (value in c(-1, NA_real_, Inf)) {
    invalid <- genes
    invalid[1, 7] <- value
    write_table(invalid, "gene")
    x$clear_cache()
    error(x$coverage(), "finite, non-negative")
  }
  write_table(genes[, 1:7], "gene")
  unlink(file.path(folder, c("partial.txt", "unique.txt")))
  single <- suppressWarnings(read_nexons(folder, ""))
  stopifnot(identical(single$coverage(), expected[1, , drop = FALSE]))
  error(single$coverage("unique"), "No unique output file")
})
