\name{read_nexons}
\alias{read_nexons}
\title{Discover nexons results without reading count tables}
\usage{read_nexons(folder, prefix = "nexons_output_")}
\arguments{
\item{folder}{Path to a directory containing nexons output files.}
\item{prefix}{Literal filename prefix. Defaults to \code{"nexons_output_"}.}
}
\description{Locate gene and isoform count tables and create an R6 object
by reading only the first line of each available file. Count rows remain unloaded.}
\value{An R6 object of class \code{nexonsR}, with read-only fields
\code{folder}, \code{prefix}, \code{files}, \code{metadata}, and \code{loaded}.
\code{files} is a named character vector of absolute paths with entries
\code{gene}, \code{partial}, and \code{unique}; missing files have \code{NA}
paths. \code{loaded} lists the tables currently cached in memory.}
\details{
The prefix is followed by \code{gene.txt}, \code{partial.txt},
or \code{unique.txt}. Other files and subdirectories are ignored.
Missing files produce a warning; if all three are missing, an error is raised.
Sample names must be non-empty, unique, and identical across available files,
regardless of order. Invalid headers or mismatched sample sets raise an error.

\code{x$get_metadata()} or \code{x$metadata} returns a data frame initially
containing one character column, \code{sample}. Samples follow the first
available file's order (gene, partial, then unique).
Add annotation columns and call \code{x$set_metadata(metadata)} to replace it.
The argument must be a data frame with unique column names and a \code{sample}
column containing each original sample exactly once (character or factor).
Rows are restored to the original sample order and annotation columns are
preserved. The setter returns the object invisibly. Invalid replacements
leave the existing metadata unchanged. The metadata field is read-only.
Count tables retain their own column order; match annotations by sample name.

An empty prefix is allowed. Prefixes cannot contain path separators.

Use \code{x$get_data("gene")} (or \code{"partial"}, \code{"unique"}) to read
and cache a table on first access. Annotation headers are validated at that
point. Annotation columns are read as character vectors and sample counts
as numeric vectors. Column names are preserved verbatim.
Cached tables are reused even if the source file changes.
\code{x$clear_cache()} releases all cached tables, so the next access reads
from disk again. Keep source files in place until read.

\code{x$gene_metadata(file = "gene")} returns a data frame containing
\code{Gene_ID} and \code{Gene_Name}. Gene rows are returned directly;
\code{"partial"} and \code{"unique"} return distinct ID/name pairs in
first-occurrence order.
\code{x$transcript_metdata(file = "unique")} returns a data frame containing
\code{Transcript_ID}, \code{Gene_ID}, and \code{Gene_Name}, retaining all rows
in source order. Only \code{"partial"} and \code{"unique"} are allowed;
\code{"gene"} raises an error. Both methods preserve source column names
and use cached tables, reading the selected file on first access.

R6 objects have reference semantics: \code{y <- x} shares the object.
Use \code{x$clone(deep = TRUE)} for an independent copy.

Use \code{x$run_deseq2(design, file = "gene", level = "gene", samples = NULL)}
to fit and return a DESeq2 model; see \code{\link{run_deseq2}}.
Count observed features with \code{x$coverage()}; see \code{\link{coverage}}.
}
\examples{
folder <- tempfile()
dir.create(folder)
header <- "Gene_ID\tGene_Name\tChr\tStart\tEnd\tStrand\tsample-1"
writeLines(c(header, "gene1\tExample\tchr1\t1\t100\t+\t5"),
           file.path(folder, "demo_gene.txt"))
x <- suppressWarnings(read_nexons(folder, prefix = "demo_"))
x$files
x$get_data("gene")
unlink(folder, recursive = TRUE)
}
