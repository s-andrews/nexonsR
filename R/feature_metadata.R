#' Extract gene annotations
#'
#' Return gene annotations from a nexons result object's selected count table.
#'
#' @name gene_metadata
#' @param file Exactly one of `"gene"` (default), `"partial"`, or `"unique"`.
#' @return A data frame containing `Gene_ID` and `Gene_Name`, preserving the
#'   original column names and character identifiers.
#' @details Call this method on an object returned by [read_nexons()].
#'   Gene files retain all rows. Partial and unique files return distinct
#'   ID/name pairs in first-occurrence order. The selected table is loaded
#'   on first access and cached in the object. Missing files raise an error.
#' @rawRd \usage{\special{x$gene_metadata(file = "gene")}}
#' @seealso [transcript_metadata()], [read_nexons()]
#' @examples
#' \dontrun{
#' x <- read_nexons("results")
#' x$gene_metadata()
#' x$gene_metadata(file = "partial")
#' }
NULL

#' Extract transcript annotations
#'
#' Return transcript and gene annotations from a nexons result object's
#' selected isoform count table.
#'
#' @name transcript_metadata
#' @param file Exactly one of `"unique"` (default) or `"partial"`.
#'   Specifying `"gene"` raises an error.
#' @return A data frame containing `Transcript_ID`, `Gene_ID`, and `Gene_Name`
#'   in that order, preserving original column names and character identifiers.
#' @details Call this method on an object returned by [read_nexons()].
#'   All rows are retained in source order. The selected table is loaded on
#'   first access and cached in the object. Missing files raise an error.
#'   The older spelling `x$transcript_metdata()` remains a compatibility alias.
#' @rawRd \usage{\special{x$transcript_metadata(file = "unique")}}
#' @seealso [gene_metadata()], [read_nexons()]
#' @examples
#' \dontrun{
#' x <- read_nexons("results")
#' x$transcript_metadata()
#' x$transcript_metadata(file = "partial")
#' }
NULL
