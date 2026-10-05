\name{coverage}
\alias{coverage}
\title{Count observed genes or isoforms per sample}
\description{Count features meeting a minimum count threshold using the
\code{coverage} method of a \code{nexonsR} object.}
\usage{
\special{x$coverage(file = "gene", level = "gene", threshold = 1)}
}
\arguments{
\item{file}{One of \code{"gene"} (default), \code{"partial"}, or \code{"unique"}.}
\item{level}{One of \code{"gene"} (default) or \code{"isoform"}.
Isoform coverage requires a partial or unique file.}
\item{threshold}{A single finite positive number, defaulting to 1.
Counts greater than or equal to this value count as observed.}
}
\value{A data frame containing all sample metadata columns, followed by a
\code{gene} or \code{isoform} column containing the number of observed features.
Metadata column types and row and column order are preserved; the object's
metadata is not modified. If the count column name already exists, a numeric
suffix is added to the new column (for example, \code{gene.1}).
Empty tables return zero counts.}
\details{
The selected table is read on first access and cached in the object.
At gene level, partial and unique counts are summed by \code{Gene_ID}
before applying the threshold. At isoform level, each \code{Transcript_ID}
is counted separately. Counts must be finite and non-negative; fractional
counts are allowed. Identifiers must be non-empty and non-missing.
Duplicate identifiers are allowed only when aggregating isoforms by gene.
}
\examples{
\dontrun{
x <- read_nexons("results")
x$coverage()
x$coverage(file = "partial", threshold = 10)
x$coverage(file = "unique", level = "isoform")
}
}
\seealso{\code{\link{read_nexons}}}
