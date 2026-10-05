library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  annotation <- data.frame(Gene_ID = c("002", "001", "002", "002"),
    Gene_Name = c("second", "first", "second", "alias"),
    Chr = "1", Start = 1, End = 10, Strand = "+", sample = 5)
  write_table <- function(data, file) utils::write.table(data,
    file.path(folder, paste0(file, ".txt")), sep = "\t", quote = FALSE,
    row.names = FALSE)
  write_table(annotation, "gene")
  for (file in c("partial", "unique")) {
    write_table(cbind(Transcript_ID = paste0(file, 1:4), annotation), file)
  }
  object <- read_nexons(folder, "")
  stopifnot(identical(object$gene_metadata(), annotation[, 1:2]),
            identical(object$loaded, "gene"))
  expected <- annotation[c(1, 2, 4), 1:2]
  rownames(expected) <- NULL
  for (file in c("partial", "unique")) {
    stopifnot(identical(object$gene_metadata(file), expected),
      identical(object$transcript_metdata(file),
        cbind(Transcript_ID = paste0(file, 1:4), annotation[, 1:2])))
  }
  stopifnot(identical(object$transcript_metdata(),
                      object$transcript_metdata("unique")))
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  error(object$transcript_metdata("gene"), "`file`")
  for (file in list("other", "g", "", NA_character_, NULL, c("partial", "unique"))) {
    error(object$gene_metadata(file), "`file`")
    error(object$transcript_metdata(file), "`file`")
  }
  # Empty inputs still return the expected annotation columns.
  for (file in c("gene", "partial", "unique")) {
    empty <- object$get_data(file)[FALSE, ]
    write_table(empty, file)
  }
  object$clear_cache()
  stopifnot(identical(dim(object$gene_metadata()), c(0L, 2L)),
            identical(dim(object$gene_metadata("partial")), c(0L, 2L)),
            identical(dim(object$transcript_metdata()), c(0L, 3L)))
  unlink(file.path(folder, "unique.txt"))
  missing_file <- suppressWarnings(read_nexons(folder, ""))
  error(missing_file$gene_metadata("unique"), "No unique output file")
  error(missing_file$transcript_metdata(), "No unique output file")
})
