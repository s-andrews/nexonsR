\name{compare_assignments}
\alias{compare_assignments}
\title{Compare quantitations across assignment files}
\description{Compare gene or isoform quantitations from available nexons
assignment files using a nexonsR object's method.}
\usage{
\special{x$compare_assignments(level = "gene", samples = NULL, units = "counts")}
}
\arguments{
\item{level}{Either \code{"gene"} (default) or \code{"isoform"}.}
\item{samples}{Optional non-empty character vector of unique known sample
names, in the desired order. \code{NULL} selects all samples in metadata order.}
\item{units}{Either \code{"counts"} (default) or \code{"log2RPM"}.}
}
\value{A data frame with one row per feature and sample. Gene output contains
\code{Gene_ID}, \code{Number_of_Isoforms}, \code{sample}, followed by available
\code{gene}, \code{partial}, and \code{unique} columns. Isoform output contains
\code{Transcript_ID}, \code{Gene_ID}, \code{sample}, \code{partial}, and
\code{unique}. Rows form sample blocks in selected sample order, with features
in first-occurrence order across gene, partial, then unique files.
Empty tables return typed zero-row output for counts.}
\details{
Gene quantitations sum raw isoform counts by gene before log2RPM conversion.
Each file uses its own sample totals and the formula
\code{log2(count / sample_total * 1e6 + 1)}. Normalization requires finite,
non-negative counts and positive sample totals, as in \code{\link{get_data}}.
Isoform counts include zero-count transcripts and use the larger distinct
transcript count from partial and unique, or the sole available count.
They are \code{NA} if neither isoform file contains the gene.

Gene output omits columns for unavailable files. Isoform output requires
both partial and unique files. Features are aligned by gene ID or the
transcript/gene pair. All features are retained; missing quantitations are
\code{NA}. Differing feature sets produce one warning per call. Missing or
empty identifiers and duplicate feature keys within a file raise errors.
Raw tables are cached through \code{\link{get_data}} and remain unchanged.
}
\examples{
\dontrun{
x <- read_nexons("results")
x$compare_assignments()
x$compare_assignments("isoform", samples = "sample-1", units = "log2RPM")
}
}
\seealso{\code{\link{read_nexons}}, \code{\link{get_data}}}
