#' Parse and query flexout frequency tables
#'
#' @name parse_flexout
#' @aliases get_flexout
#' @param samples Unique sample names including their original `.bam` suffix.
#'   `NULL` parses all samples with discovered flexout files.
#' @param chunk_size Maximum number of data lines read at once (default 100000).
#' @param transcript_id A single transcript ID.
#' @param end Either `"Start_Flex"` or `"End_Flex"`.
#' @details `x$parse_flexout()` reads gzipped files in bounded chunks and emits
#'   a progress message for each uncached sample. Completed samples are cached;
#'   a failed sample is not cached. Repeated calls skip completed samples.
#'   `x$flexout_loaded` lists cached samples. `x$clear_cache()` releases both
#'   count tables and parsed flexout data, retaining discovered file paths.
#'
#'   Each sample has a transcript index with separate sparse frequency tables
#'   for the two ends. Only observed values are stored, as signed R integers;
#'   counts use doubles to avoid 32-bit count overflow. Memory scales with
#'   distinct transcript/end/value combinations plus the input chunk, rather
#'   than the number of reads or the width of the value range. Start/end pairing
#'   is intentionally discarded. Values must fit R's non-missing integer range
#'   (-2147483647 to 2147483647).
#'
#'   `x$get_flexout(transcript_id, end, samples)` queries already parsed samples;
#'   `samples = NULL` selects all parsed samples. It never triggers parsing.
#'   An absent transcript returns an empty frequency table for that sample.
#' @return `parse_flexout()` returns the object invisibly. `get_flexout()` returns
#'   a named list of data frames (one per sample), with integer `value` and
#'   numeric `count` columns sorted by value, without expanding repeated reads.
#' @rawRd \usage{\special{x$parse_flexout(samples = NULL, chunk_size = 100000L)
#' x$get_flexout(transcript_id, end = c("Start_Flex", "End_Flex"), samples = NULL)}}
NULL

flexout_samples <- function(samples, available) {
  if (is.null(samples)) return(available)
  if (!is.character(samples) || !length(samples) || anyNA(samples) ||
      anyDuplicated(samples) || !all(samples %in% available)) {
    stop("`samples` must contain unique available sample names.", call. = FALSE)
  }
  samples
}

empty_flexout_counts <- function() data.frame(value = integer(), count = numeric())

add_flexout_counts <- function(previous, values) {
  observed <- sort(unique(values))
  counts <- as.double(tabulate(match(values, observed), nbins = length(observed)))
  if (is.null(previous)) return(data.frame(value = observed, count = counts))
  positions <- match(observed, previous$value)
  found <- !is.na(positions)
  previous$count[positions[found]] <- previous$count[positions[found]] + counts[found]
  if (any(!found)) {
    previous <- rbind(previous, data.frame(value = observed[!found], count = counts[!found]))
  }
  previous
}

parse_flexout_file <- function(path, chunk_size) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  header <- readLines(connection, n = 1L, warn = FALSE)
  if (!identical(header, "Transcript_ID\tStart_Flex\tEnd_Flex")) {
    stop("Invalid flexout header in ", basename(path), ".", call. = FALSE)
  }
  transcripts <- new.env(hash = TRUE, parent = emptyenv())
  offset <- 1
  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (!length(lines)) break
    # Sentinel ensures an empty trailing field is not silently discarded.
    fields <- strsplit(paste0(lines, "\t."), "\t", fixed = TRUE)
    valid <- lengths(fields) == 4L
    if (!all(valid)) {
      stop("Expected three columns in ", basename(path), " at line ",
           offset + which(!valid)[1L], ".", call. = FALSE)
    }
    fields <- matrix(unlist(fields, use.names = FALSE), nrow = 4L)
    ids <- fields[1L, ]
    start <- suppressWarnings(as.integer(fields[2L, ]))
    end <- suppressWarnings(as.integer(fields[3L, ]))
    valid <- nzchar(trimws(ids)) & grepl("^[+-]?[0-9]+$", fields[2L, ]) &
      grepl("^[+-]?[0-9]+$", fields[3L, ]) & !is.na(start) & !is.na(end)
    if (!all(valid)) {
      stop("Invalid transcript ID or integer flex value in ", basename(path),
           " at line ", offset + which(!valid)[1L], ".", call. = FALSE)
    }
    groups <- split(seq_along(ids), ids)
    for (id in names(groups)) {
      rows <- groups[[id]]
      previous <- transcripts[[id]]
      transcripts[[id]] <- list(
        Start_Flex = add_flexout_counts(previous$Start_Flex, start[rows]),
        End_Flex = add_flexout_counts(previous$End_Flex, end[rows]))
    }
    offset <- offset + length(lines)
  }
  # Store lists, not environments, so deep clones and returned tables remain independent.
  as.list(transcripts, all.names = TRUE)
}
