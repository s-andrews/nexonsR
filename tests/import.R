library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  types <- c("gene", "partial", "unique")
  paths <- file.path(folder, paste0("nexons_output_", types, ".txt"))
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  error(read_nexons(folder), "No nexons output files")
  file.create(file.path(folder, "nexons_output_other.txt"))
  dir.create(paths[1])
  error(read_nexons(folder), "No nexons output files")
  unlink(paths[1], recursive = TRUE)
  # Invalid contents must not be read during discovery or printing.
  writeLines("not a count table", paths[1])
  warnings <- character()
  x <- withCallingHandlers(read_nexons(folder), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  stopifnot(inherits(x, "nexonsR"), length(warnings) == 1L,
    grepl(basename(paths[2]), warnings, fixed = TRUE),
    grepl(basename(paths[3]), warnings, fixed = TRUE),
    identical(names(x$files), types), sum(is.na(x$files)) == 2L,
    length(x$loaded) == 0L)
  invisible(capture.output(print(x)))
  stopifnot(length(x$loaded) == 0L)
  error(x$get_data("gene"), "Invalid header")
  error(x$get_data("partial"), "No partial output file")
  header <- "Gene_ID\tGene_Name\tChr\tStart\tEnd\tStrand\t01 sample-A\tsample.B"
  row <- "001\tO'Brien#gene\t01\t1\t20\t+\t3\t0"
  writeLines(c(header, row), paths[1])
  for (path in paths[-1]) {
    writeLines(c(paste0("Transcript_ID\t", header), paste0("002\t", row)), path)
  }
  x <- read_nexons(folder)
  gene <- x$get_data("gene")
  stopifnot(identical(gene$Gene_ID, "001"),
    identical(gene$Gene_Name, "O'Brien#gene"),
    identical(names(gene)[7:8], c("01 sample-A", "sample.B")),
    identical(gene[[7]], 3), identical(x$loaded, "gene"))
  stopifnot(identical(x$get_data("partial")$Transcript_ID, "002"),
    identical(x$get_data("unique")$Transcript_ID, "002"))
  y <- read_nexons(folder)
  stopifnot(length(y$loaded) == 0L)
  unlink(paths[1])
  stopifnot(identical(x$get_data("gene"), gene))
  clone <- x$clone(deep = TRUE)
  x$clear_cache()
  stopifnot(length(x$loaded) == 0L, length(clone$loaded) == 3L)
  suppressWarnings(error(x$get_data("gene"), "cannot open"))
  # Literal metacharacters in prefixes, relative paths, and empty prefixes.
  custom <- "run[1].+_"
  for (type in types) file.create(file.path(folder, paste0(custom, type, ".txt")))
  custom_result <- read_nexons(folder, custom)
  stopifnot(all(!is.na(custom_result$files)), length(custom_result$loaded) == 0L)
  for (type in types) file.create(file.path(folder, paste0(type, ".txt")))
  stopifnot(all(!is.na(read_nexons(folder, "")$files)))
  error(read_nexons(NA_character_), "`folder`")
  error(read_nexons(c(folder, folder)), "`folder`")
  error(read_nexons(file.path(folder, "absent")), "Directory does not exist")
  error(read_nexons(folder, NA_character_), "`prefix`")
  error(read_nexons(folder, "../"), "path separators")
})
