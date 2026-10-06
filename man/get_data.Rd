\name{get_data}
\alias{get_data}
\title{Extract nexons counts or log2 reads per million}
\description{Read and cache a nexons count table, optionally selecting samples
and converting their raw counts to log2 reads per million (log2RPM).}
\usage{
\special{x$get_data(type = c("gene", "partial", "unique"),
  samples = NULL, units = c("counts", "log2RPM"))}
}
\arguments{
\item{type}{Exactly one of \code{"gene"} (default), \code{"partial"}, or
\code{"unique"}.}
\item{samples}{Optional non-empty character vector of unique sample names.
Only these sample columns are returned, in the supplied order. \code{NULL}
returns all samples in the source table's order.}
\item{units}{Either \code{"counts"} (the default) for raw counts or
\code{"log2RPM"}.}
}
\value{A data frame containing the annotation columns for \code{type}, followed
by the selected sample columns. Annotation values and row order are preserved.
Sample columns contain raw counts when \code{units = "counts"} and log2RPM
values otherwise.}
\details{Call this method on an object returned by \code{read_nexons()}.
The full raw table is read and cached on first access. Sample selection and
log2RPM conversion are applied to the returned copy and do not alter the cache.

For each selected sample and row, log2RPM is calculated as
\code{log2(count / sum(counts) * 1e6 + 1)}, where the sum is over all rows of
the selected \code{type} table. Counts must be finite and non-negative, and
every selected sample must have a positive total count.}
\examples{
\dontrun{
x <- read_nexons("results")
x$get_data("gene")
x$get_data("partial", samples = c("sample-2", "sample-1"))
x$get_data("unique", samples = "sample-1", units = "log2RPM")
}
}
\seealso{\code{\link{read_nexons}}}
