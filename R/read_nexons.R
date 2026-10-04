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
#' Discovery does not open the files or validate their contents.
#'
#' Use `x$get_data("gene")` (or `"partial"`, `"unique"`) to read and cache a
#' table on first access. Required annotation headers are validated then.
#' Sample names are preserved verbatim. Cached tables are reused, even if the
#' source file changes. `x$clear_cache()` releases all cached tables so the next
#' access reads from disk again. Keep source files in place until read.
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
    },
    get_data = function(type = c("gene", "partial", "unique")) {
      type <- match.arg(type)
      if (!is.null(private$.cache[[type]])) return(private$.cache[[type]])
      path <- private$.files[[type]]
      if (is.na(path)) {
        stop("No ", type, " output file was found for this object.", call. = FALSE)
      }
      annotation <- c("Gene_ID", "Gene_Name", "Chr", "Start", "End", "Strand")
      if (type != "gene") annotation <- c("Transcript_ID", annotation)
      header <- names(utils::read.delim(path, nrows = 0L, check.names = FALSE,
        quote = "", comment.char = "", row.names = NULL))
      if (length(header) <= length(annotation) ||
          !identical(header[seq_along(annotation)], annotation)) {
        stop("Invalid header in ", basename(path), ": expected ",
             paste(annotation, collapse = ", "),
             ", followed by sample columns.", call. = FALSE)
      }
      # Keep identifiers (including numeric-looking IDs) as text.
      classes <- c(rep("character", length(annotation)),
                   rep("numeric", length(header) - length(annotation)))
      data <- utils::read.delim(path, check.names = FALSE, quote = "",
        comment.char = "", stringsAsFactors = FALSE, row.names = NULL,
        colClasses = classes, fill = FALSE)
      private$.cache[[type]] <- data
      data
    },
    clear_cache = function() {
      private$.cache <- list()
      invisible(self)
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
    loaded = function() names(private$.cache)
  ),
  private = list(.folder = NULL, .prefix = NULL, .files = NULL, .cache = NULL)
)
