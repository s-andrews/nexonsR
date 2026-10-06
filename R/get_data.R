#' Extract nexons counts or log2 reads per million
#'
#' Read and cache a nexons count table, optionally selecting samples and
#' converting their raw counts to log2 reads per million (log2RPM).
#'
#' @name get_data
#' @param type Exactly one of `"gene"` (default), `"partial"`, or `"unique"`.
#' @param samples Optional non-empty character vector of unique sample names.
#'   Only these sample columns are returned, in the supplied order. `NULL`
#'   returns all samples in the source table's order.
#' @param units Either `"counts"` (the default) for raw counts or `"log2RPM"`.
#' @return A data frame containing the annotation columns for `type`, followed
#'   by the selected sample columns. Annotation values and row order are
#'   preserved. Sample columns contain raw counts when `units = "counts"` and
#'   log2RPM values otherwise.
#' @details Call this method on an object returned by [read_nexons()]. The full
#'   raw table is read and cached on first access. Sample selection and log2RPM
#'   conversion are applied to the returned copy and do not alter the cache.
#'
#'   For each selected sample and row, log2RPM is calculated as
#'   `log2(count / sum(counts) * 1e6 + 1)`, where the sum is over all rows of
#'   the selected `type` table. Counts must be finite and non-negative, and
#'   every selected sample must have a positive total count.
#' @rawRd \usage{\special{x$get_data(type = c("gene", "partial", "unique"),
#'   samples = NULL, units = c("counts", "log2RPM"))}}
#' @seealso [read_nexons()]
#' @examples
#' \dontrun{
#' x <- read_nexons("results")
#' x$get_data("gene")
#' x$get_data("partial", samples = c("sample-2", "sample-1"))
#' x$get_data("unique", samples = "sample-1", units = "log2RPM")
#' }
NULL
