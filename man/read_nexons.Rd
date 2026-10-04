\name{read_nexons}
\alias{read_nexons}
\title{Discover nexons results without reading count tables}
\usage{read_nexons(folder, prefix = "nexons_output_")}
\arguments{
\item{folder}{Path to a directory containing nexons output files.}
\item{prefix}{Literal filename prefix. Defaults to \code{"nexons_output_"}.}
}
\description{Locate gene and isoform count tables and create an R6 object
without opening the files or reading their contents.}
\value{An R6 object of class \code{nexonsR}, with read-only fields
\code{folder}, \code{prefix}, \code{files}, and \code{loaded}.
\code{files} is a named character vector of absolute paths with entries
\code{gene}, \code{partial}, and \code{unique}; missing files have \code{NA}
paths. \code{loaded} lists the tables currently cached in memory.}
\details{
The prefix is followed by \code{gene.txt}, \code{partial.txt},
or \code{unique.txt}. Other files and subdirectories are ignored.
Missing files produce a warning; if all three are missing, an error is raised.
An empty prefix is allowed. Prefixes cannot contain path separators.

Use \code{x$get_data("gene")} (or \code{"partial"}, \code{"unique"}) to read
and cache a table on first access. Annotation headers are validated at that
point. Annotation columns are read as character vectors and sample counts
as numeric vectors. Column names are preserved verbatim.
Cached tables are reused even if the source file changes.
\code{x$clear_cache()} releases all cached tables, so the next access reads
from disk again. Keep source files in place until read.

R6 objects have reference semantics: \code{y <- x} shares the object.
Use \code{x$clone(deep = TRUE)} for an independent copy.
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
