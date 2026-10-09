\name{parse_flexout}
\alias{parse_flexout}
\alias{get_flexout}
\title{Parse and query flexout frequency tables}
\description{Read gzipped flexout files in bounded chunks and cache separate
sparse frequency tables for start and end values, indexed by sample and transcript.}
\usage{
\special{x$parse_flexout(samples = NULL, chunk_size = 100000L)
x$get_flexout(transcript_id, end = c("Start_Flex", "End_Flex"), samples = NULL)}
}
\arguments{
\item{samples}{Unique original sample names, including any trailing \code{.bam}.
For parsing, \code{NULL} selects all samples with discovered files. For queries,
\code{NULL} selects all parsed samples. Explicit unavailable samples raise an error.}
\item{chunk_size}{Maximum number of data lines read at once; a positive integer.}
\item{transcript_id}{One non-empty transcript ID.}
\item{end}{Either \code{"Start_Flex"} or \code{"End_Flex"}.}
}
\value{\code{parse_flexout()} returns the object invisibly.
\code{get_flexout()} returns a named list of data frames, one per sample, with
integer \code{value} and numeric \code{count} columns sorted by value.
An absent transcript returns an empty table. Queries never trigger parsing.}
\details{
Files must have exactly three tab-separated columns with the header
\code{Transcript_ID}, \code{Start_Flex}, \code{End_Flex}, in that order.
Transcript IDs must be non-empty. Flex values must be signed integers within
R's non-missing integer range (-2147483647 to 2147483647).
Only observed values are stored. Counts use doubles to avoid 32-bit count
overflow. Memory scales with distinct transcript/end/value combinations plus
the input chunk, rather than the number of reads or width of the value range.
Start/end pairing is intentionally discarded.

A progress message is emitted for each uncached sample. Completed samples are
cached; failed samples are not. Repeated calls skip completed samples.
Read-only \code{x$flexout_loaded} lists cached samples.
\code{x$clear_cache()} releases both count tables and parsed flexout data,
retaining discovered paths. Deep clones have independent caches.
}
\examples{
\dontrun{
x <- read_nexons("results")
x$parse_flexout()
x$get_flexout("transcript-1", "Start_Flex")
x$get_flexout("transcript-1", "End_Flex", samples = "sample-1.bam")
}
}
\seealso{\code{\link{read_nexons}}}
