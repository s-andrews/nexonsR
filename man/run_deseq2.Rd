\name{run_deseq2}
\alias{run_deseq2}
\title{Fit a DESeq2 model from nexons counts and metadata}
\description{The \code{run_deseq2} method of a \code{nexonsR} object prepares
counts and sample metadata, constructs a DESeq2 dataset, and runs
\code{DESeq2::DESeq()} with its default settings.}
\usage{
\special{x$run_deseq2(design, file = "gene", level = "gene", samples = NULL)}
}
\arguments{
\item{design}{The design passed to \code{DESeq2::DESeqDataSetFromMatrix},
usually a formula such as \code{~ condition} or \code{~ batch + condition}.
Variables refer to columns in the object's metadata. If supplying a model
matrix, its rows must already correspond to the selected sample order.}
\item{file}{Exactly one of \code{"gene"} (default), \code{"partial"}, or
\code{"unique"}. The selected file must have been found at object creation.}
\item{level}{Exactly \code{"gene"} (default) or \code{"isoform"}.
Gene-level analysis of partial or unique tables sums counts by \code{Gene_ID}.
Isoform analysis uses \code{Transcript_ID} and requires partial or unique data.}
\item{samples}{Optional non-empty character vector of unique sample names.
Both counts and metadata are selected and ordered by this vector.
\code{NULL} uses all samples in metadata order.}
}
\value{A fitted \code{DESeqDataSet}, suitable for \code{DESeq2::results()}.
The model is not stored in the nexonsR object.}
\details{
Call this method on an object returned by \code{read_nexons()}.
DESeq2 is an optional dependency; install it using
\code{BiocManager::install("DESeq2")} before running the method.
Only the selected table is lazily read and cached. Metadata is not modified.
Unused factor levels are dropped from the selected metadata.
Count columns and metadata rows are aligned by sample name.

Counts must be finite, non-negative whole numbers and fit in R's integer
range, including after aggregation. Counts are never rounded or normalized
before fitting. Feature identifiers must be non-missing and non-empty;
duplicates are allowed only when aggregating isoforms by gene.
DESeq2 performs its own design and model checks and propagates fitting errors.
}
\examples{
\dontrun{
x <- read_nexons("results")
metadata <- x$get_metadata()
# Supply annotations appropriate to your experiment:
metadata$condition <- factor(c("control", "control", "treated", "treated"))
x$set_metadata(metadata)
model <- x$run_deseq2(~ condition)
result <- DESeq2::results(model)
isoform_model <- x$run_deseq2(~ condition, file = "unique", level = "isoform")
}
}
\seealso{\code{\link{read_nexons}}}
