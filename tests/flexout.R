library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- c("PBM75668_pass_L001_barcode01.bam", "sample.bam.middle.bam")
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  for (prefix in c("nexons_output_", "run[1].+_", "")) {
    write_headers <- function(names) {
      header <- paste(c("Gene_ID", "Gene_Name", "Chr", "Start", "End", "Strand",
                        names), collapse = "\t")
      for (type in c("gene", "partial", "unique")) {
        writeLines(if (type == "gene") header else paste0("Transcript_ID\t", header),
                   file.path(folder, paste0(prefix, type, ".txt")))
      }
    }
    write_headers(samples)
    paths <- setNames(file.path(folder, paste0(prefix,
      sub("\\.bam$", "", samples), "_flexout.txt.gz")), samples)
    stopifnot(all(is.na(read_nexons(folder, prefix)$flexout_files)))
    # Invalid gzip content proves discovery does not open these files.
    for (path in paths) writeLines("not gzip data", path)
    x <- read_nexons(folder, prefix)
    stopifnot(identical(x$flexout_files, paths), length(x$loaded) == 0L)
    x$get_data("gene")
    x$clear_cache()
    stopifnot(identical(x$flexout_files, paths), length(x$loaded) == 0L)
    error(x$flexout_files <- character(), "unused argument")
    unlink(paths[2])
    dir.create(paths[2])
    warning_text <- character()
    incomplete <- withCallingHandlers(read_nexons(folder, prefix), warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    stopifnot(length(warning_text) == 1L, is.na(incomplete$flexout_files[2]),
              identical(x$clone(deep = TRUE)$flexout_files, paths))
    extra <- file.path(folder, paste0(prefix, "unknown_flexout.txt.gz"))
    file.create(extra)
    error(read_nexons(folder, prefix), "Flexout sample names do not match")
    unlink(extra)
    write_headers(c("sample", "sample.bam"))
    error(read_nexons(folder, prefix), "Ambiguous flexout sample names")
    unlink(list.files(folder, full.names = TRUE), recursive = TRUE)
  }
})
