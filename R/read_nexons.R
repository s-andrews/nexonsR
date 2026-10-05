#' Discover nexons results without reading count tables
#'
#' @param folder Path to a directory containing nexons output files.
#' @param prefix Literal filename prefix. Defaults to `"nexons_output_"`.
#' @return An R6 object of class `nexonsR`. Read-only fields `folder`, `prefix`,
#'   and `files` describe the input. `files` is a named character vector with
#'   entries `gene`, `partial`, and `unique`; missing files have `NA` paths.
#'   `loaded` lists the tables currently cached in memory.
#' @details
#' The prefix is followed by `gene.txt`, `partial.txt`, or
#' `unique.txt`. No other files or subdirectories are searched.
#' Missing files produce a warning; if all three are missing, an error is raised.
#' Discovery reads only the first line of each available file. Sample names
#' must be non-empty, unique, and identical across files, regardless of order.
#' Count rows are not read until requested.
#'
#' `x$metadata` or `x$get_metadata()` returns a data frame initially containing
#' one character column, `sample`. Samples follow the first available file's
#' order (gene, partial, then unique). Add annotation columns and pass the
#' modified data frame to `x$set_metadata(metadata)`. It must contain each
#' original sample exactly once; rows are restored to the original order.
#' The `metadata` field is read-only; use the setter to replace it.
#'
#' Use `x$get_data("gene")` (or `"partial"`, `"unique"`) to read and cache a
#' table on first access. Required annotation headers are validated then.
#' Sample names are preserved verbatim. Cached tables are reused, even if the
#' source file changes. `x$clear_cache()` releases all cached tables so the next
#' access reads from disk again. Keep source files in place until read.
#'
#' `x$gene_metadata(file = "gene")` returns `Gene_ID` and `Gene_Name`.
#' Gene rows are returned directly; `"partial"` and `"unique"` return distinct
#' ID/name pairs in first-occurrence order.
#' `x$transcript_metdata(file = "unique")` returns `Transcript_ID`, `Gene_ID`,
#' and `Gene_Name`, retaining all rows. Only `"partial"` and `"unique"` are
#' allowed. Both methods preserve source column names and use cached tables.
#'
#' Run `x$run_deseq2(design, file = "gene", level = "gene", samples = NULL)`
#' to fit and return a DESeq2 model. See [run_deseq2] for method arguments.
#'
#' R6 objects have reference semantics: assigning `y <- x` shares the object.
#' Use `x$clone(deep = TRUE)` for an independent copy.
#' @examples
#' folder <- tempfile()
#' dir.create(folder)
#' header <- "Gene_ID\tGene_Name\tChr\tStart\tEnd\tStrand\tsample-1"
#' writeLines(c(header, "gene1\tExample\tchr1\t1\t100\t+\t5"),
#'            file.path(folder, "demo_gene.txt"))
#' # Warns because the two isoform files are absent:
#' x <- suppressWarnings(read_nexons(folder, prefix = "demo_"))
#' x$files
#' x$get_data("gene")
#' unlink(folder, recursive = TRUE)
#' @export
#' @importFrom R6 R6Class
#' @importFrom stats setNames
#' @importFrom utils read.delim
read_nexons <- function(folder, prefix = "nexons_output_") {
  NexonsResults$new(folder, prefix)
}

validate_string <- function(value, name, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      (!allow_empty && !nzchar(value))) {
    stop("`", name, "` must be a single ",
         if (allow_empty) "non-missing" else "non-empty",
         " character string.", call. = FALSE)
  }
}

read_nexons_header <- function(path, type) {
  annotation <- c("Gene_ID", "Gene_Name", "Chr", "Start", "End", "Strand")
  if (type != "gene") annotation <- c("Transcript_ID", annotation)
  line <- readLines(path, n = 1L, warn = FALSE)
  # Appending a sentinel preserves empty trailing fields when splitting.
  header <- if (length(line)) {
    fields <- strsplit(paste0(line, "\t."), "\t", fixed = TRUE)[[1L]]
    fields[-length(fields)]
  } else character()
  if (length(header) <= length(annotation) ||
      !identical(header[seq_along(annotation)], annotation)) {
    stop("Invalid header in ", basename(path), ": expected ",
         paste(annotation, collapse = ", "),
         ", followed by sample columns.", call. = FALSE)
  }
  samples <- header[-seq_along(annotation)]
  if (any(!nzchar(trimws(samples))) || anyDuplicated(samples)) {
    stop("Sample names must be non-empty and unique in ", basename(path),
         ".", call. = FALSE)
  }
  list(annotation = annotation, samples = samples)
}

NexonsResults <- R6::R6Class(
  "nexonsR",
  public = list(
    initialize = function(folder, prefix = "nexons_output_") {
      validate_string(folder, "folder")
      validate_string(prefix, "prefix", allow_empty = TRUE)
      if (grepl("/", prefix, fixed = TRUE) ||
          grepl("\\", prefix, fixed = TRUE)) {
        stop("`prefix` must be a filename prefix without path separators.",
             call. = FALSE)
      }
      if (!dir.exists(folder)) {
        stop("Directory does not exist: ", folder, call. = FALSE)
      }
      private$.folder <- normalizePath(folder, winslash = "/", mustWork = TRUE)
      private$.prefix <- prefix
      types <- c("gene", "partial", "unique")
      paths <- stats::setNames(file.path(private$.folder,
        paste0(prefix, types, ".txt")), types)
      found <- file.exists(paths) & !dir.exists(paths)
      missing_files <- paste(basename(paths[!found]), collapse = ", ")
      if (!any(found)) {
        stop("No nexons output files found in ", private$.folder,
             ". Expected: ", missing_files, call. = FALSE)
      }
      paths[!found] <- NA_character_
      private$.files <- paths
      private$.cache <- list()
      if (!all(found)) {
        warning("Missing nexons output files: ", missing_files, call. = FALSE)
      }
      samples <- NULL
      for (type in types[found]) {
        current <- read_nexons_header(paths[[type]], type)$samples
        if (is.null(samples)) samples <- current
        if (!setequal(samples, current)) {
          stop("Sample names do not match in ", basename(paths[[type]]),
               ". Missing: ", paste(setdiff(samples, current), collapse = ", "),
               "; unexpected: ", paste(setdiff(current, samples), collapse = ", "),
               call. = FALSE)
        }
      }
      private$.metadata <- data.frame(sample = samples, stringsAsFactors = FALSE)
    },
    get_metadata = function() {
      private$.metadata
    },
    set_metadata = function(metadata) {
      if (!is.data.frame(metadata) || anyDuplicated(names(metadata)) ||
          !"sample" %in% names(metadata)) {
        stop("`metadata` must be a data frame with unique column names and a `sample` column.",
             call. = FALSE)
      }
      samples <- metadata[["sample"]]
      if (!(is.character(samples) || is.factor(samples)) ||
          anyNA(samples) || anyDuplicated(samples) ||
          !setequal(as.character(samples), private$.metadata$sample)) {
        stop("`metadata$sample` must contain each original sample exactly once.",
             call. = FALSE)
      }
      metadata <- as.data.frame(metadata)
      metadata$sample <- as.character(samples)
      metadata <- metadata[match(private$.metadata$sample, metadata$sample), , drop = FALSE]
      rownames(metadata) <- NULL
      private$.metadata <- metadata
      invisible(self)
    },
    get_data = function(type = c("gene", "partial", "unique")) {
      type <- match.arg(type)
      if (!is.null(private$.cache[[type]])) return(private$.cache[[type]])
      path <- private$.files[[type]]
      if (is.na(path)) {
        stop("No ", type, " output file was found for this object.", call. = FALSE)
      }
      header <- read_nexons_header(path, type)
      annotation <- header$annotation
      if (!setequal(header$samples, private$.metadata$sample)) {
        stop("Sample names no longer match metadata in ", basename(path),
             ".", call. = FALSE)
      }
      # Keep identifiers (including numeric-looking IDs) as text.
      classes <- c(rep("character", length(annotation)),
                   rep("numeric", length(header$samples)))
      data <- utils::read.delim(path, check.names = FALSE, quote = "",
        comment.char = "", stringsAsFactors = FALSE, row.names = NULL,
        colClasses = classes, fill = FALSE)
      private$.cache[[type]] <- data
      data
    },
    gene_metadata = function(file = "gene") {
      validate_string(file, "file")
      if (!file %in% c("gene", "partial", "unique")) {
        stop("`file` must be 'gene', 'partial', or 'unique'.", call. = FALSE)
      }
      metadata <- self$get_data(file)[, c("Gene_ID", "Gene_Name"), drop = FALSE]
      if (file != "gene") metadata <- unique(metadata)
      rownames(metadata) <- NULL
      metadata
    },
    transcript_metdata = function(file = "unique") {
      validate_string(file, "file")
      if (!file %in% c("partial", "unique")) {
        stop("`file` must be 'partial' or 'unique'; gene files have no transcript metadata.",
             call. = FALSE)
      }
      self$get_data(file)[, c("Transcript_ID", "Gene_ID", "Gene_Name"), drop = FALSE]
    },
    clear_cache = function() {
      private$.cache <- list()
      invisible(self)
    },
    run_deseq2 = function(design, file = "gene", level = "gene", samples = NULL) {
      if (missing(design)) stop("`design` must be supplied.", call. = FALSE)
      input <- prepare_deseq2(self, file, level, samples)
      if (!requireNamespace("DESeq2", quietly = TRUE)) {
        stop("Install DESeq2 to use this method: BiocManager::install(\"DESeq2\").",
             call. = FALSE)
      }
      dataset <- DESeq2::DESeqDataSetFromMatrix(
        countData = input$counts, colData = input$metadata, design = design)
      DESeq2::DESeq(dataset)
    },
    print = function(...) {
      cat("<nexonsR>\n", "Folder: ", private$.folder, "\n",
          "Prefix: ", private$.prefix, "\n", sep = "")
      for (type in names(private$.files)) {
        status <- if (is.na(private$.files[[type]])) "missing" else
          if (type %in% names(private$.cache)) "cached" else "available (not loaded)"
        cat("  ", type, ": ", status, "\n", sep = "")
      }
      invisible(self)
    }
  ),
  active = list(
    folder = function() private$.folder,
    prefix = function() private$.prefix,
    files = function() private$.files,
    metadata = function() private$.metadata,
    loaded = function() names(private$.cache)
  ),
  private = list(.folder = NULL, .prefix = NULL, .files = NULL, .cache = NULL,
                 .metadata = NULL)
)
