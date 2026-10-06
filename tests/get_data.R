library(nexonsR)

local({
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- c("sample 1", "sample-2", "sample.3")
  annotation <- data.frame(Gene_ID = c("g1", "g2"), Gene_Name = c("one", "two"),
    Chr = "1", Start = c(1, 20), End = c(10, 30), Strand = "+")
  counts <- matrix(c(10, 30, 0, 20, 5, 5), nrow = 2,
                   dimnames = list(NULL, samples))
  write_table <- function(data, type) utils::write.table(data,
    file.path(folder, paste0(type, ".txt")), sep = "\t", quote = FALSE,
    row.names = FALSE)
  write_table(cbind(annotation, counts), "gene")
  for (type in c("partial", "unique")) {
    write_table(cbind(Transcript_ID = c("t1", "t2"), annotation, counts), type)
  }
  object <- read_nexons(folder, "")
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }

  raw <- object$get_data("gene")
  selected <- object$get_data("gene", samples = samples[c(3, 1)])
  stopifnot(identical(names(selected), c(names(annotation), samples[c(3, 1)])),
            identical(selected[[samples[3]]], counts[, 3]),
            identical(selected[[samples[1]]], counts[, 1]),
            identical(object$get_data("gene", units = "counts"), raw),
            identical(object$loaded, "gene"))

  transformed <- object$get_data("gene", samples = samples[c(2, 1)], units = "log2RPM")
  expected <- log2(sweep(counts[, c(2, 1), drop = FALSE], 2,
                         colSums(counts[, c(2, 1), drop = FALSE]), "/") * 1e6 + 1)
  stopifnot(identical(as.matrix(transformed[, samples[c(2, 1)]]), expected),
            identical(transformed[, names(annotation)], raw[, names(annotation)]),
            identical(object$get_data("gene"), raw))

  for (type in c("partial", "unique")) {
    result <- object$get_data(type, samples = "sample-2", units = "log2RPM")
    stopifnot(ncol(result) == 8L, identical(names(result)[8], "sample-2"),
              isTRUE(all.equal(result[[8]], as.numeric(expected[, 1]))))
  }
  for (selection in list(character(), "unknown", c(samples[1], samples[1]),
                         NA_character_, 1)) {
    error(object$get_data("gene", samples = selection), "`samples`")
  }
  error(object$get_data("gene", units = "rpm"), "arg")

  zero_counts <- counts
  zero_counts[, 2] <- 0
  write_table(cbind(annotation, zero_counts), "gene")
  object$clear_cache()
  error(object$get_data("gene", samples = "sample-2", units = "log2RPM"),
        "zero total counts")
  stopifnot(identical(object$get_data("gene", samples = "sample-2")[[7]], c(0, 0)))

  invalid_counts <- counts
  invalid_counts[1, 1] <- -1
  write_table(cbind(annotation, invalid_counts), "gene")
  object$clear_cache()
  error(object$get_data("gene", samples = "sample 1", units = "log2RPM"),
        "finite, non-negative")
})
