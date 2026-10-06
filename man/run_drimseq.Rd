\name{run_drimseq}
\alias{run_drimseq}
\title{Fit a DRIMSeq isoform usage model from nexons counts}
\description{The \code{run_drimseq} method filters isoform counts, estimates
precision and fits a DRIMSeq model. Hypothesis testing is performed separately.}
\usage{
\special{x$run_drimseq(design, file = "unique", samples = NULL,
  min_samps_feature_expr, min_feature_expr = 10,
  min_samps_feature_prop = min_samps_feature_expr,
  min_feature_prop = 0.1, min_samps_gene_expr = NULL,
  min_gene_expr = 10, run_gene_twice = FALSE)}
}
\arguments{
\item{design}{A one-sided formula using sample metadata, or a numeric model
matrix with rows already in selected sample order. Designs must be finite,
have full column rank, and leave residual degrees of freedom. Unnamed matrix
columns are named \code{coef1}, \code{coef2}, and so on.}
\item{file}{Exactly \code{"unique"} (default) or \code{"partial"}.
Both contain uniquely assigned reads and are handled identically.}
\item{samples}{Optional non-empty character vector of unique known sample names.
Counts and metadata are selected and ordered accordingly. \code{NULL} uses all
samples in metadata order. Unused factor levels are dropped.}
\item{min_samps_feature_expr}{Required number of samples meeting the isoform
count threshold. For a simple group comparison, use the smallest group size.
Groups are not inferred from the design.}
\item{min_feature_expr}{Minimum isoform count, default 10.}
\item{min_samps_feature_prop}{Number of samples meeting the isoform proportion
threshold. Defaults to \code{min_samps_feature_expr}.}
\item{min_feature_prop}{Minimum within-gene isoform proportion, default 0.1.
Must be between zero and one. Set to zero to disable the positive proportion
threshold, or lower it to retain minor isoforms.}
\item{min_samps_gene_expr}{Number of samples meeting the gene count threshold.
\code{NULL} (default) means all selected samples.}
\item{min_gene_expr}{Minimum summed isoform count per gene, default 10.}
\item{run_gene_twice}{Whether to repeat gene filtering after isoform removal.
Defaults to \code{FALSE}.}
}
\value{A native fitted \code{dmDSfit} object, suitable for
\code{DRIMSeq::dmTest()}. It is not stored in the nexons object.}
\details{
DRIMSeq is optional; install it with \code{BiocManager::install("DRIMSeq")}.
Only the selected table is lazily read and cached. Object metadata is unchanged.
Counts must be finite, non-negative whole numbers within R's integer range;
they are never rounded or normalized. Gene and transcript IDs must be nonempty
and non-missing, and transcript IDs must be unique.

Sample thresholds must be whole numbers from zero through the selected sample
count. Count thresholds must be finite non-negative numbers. DRIMSeq removes
unobserved isoforms and genes with fewer than two retained isoforms. An error
is raised if no analysable genes remain. For large cohorts, consider lowering
\code{min_samps_gene_expr} to tolerate occasional poorly observed samples.
The defaults are starting points, not nanopore-specific calibration.

The workflow calls \code{dmDSdata()}, \code{dmFilter()}, \code{dmPrecision()}
and \code{dmFit()}, leaving precision and fitting settings at DRIMSeq defaults.
Call \code{set.seed()} before fitting for reproducibility. Choose the hypothesis
after fitting using \code{dmTest()}'s coefficient, contrast or reduced-design
arguments, then extract gene or feature results with \code{DRIMSeq::results()}.

Returned sample metadata retains original names in \code{sample};
\code{sample_id} is reserved for DRIMSeq's sample identifiers. If sample names
collide with \code{gene_id} or \code{feature_id}, internal sample identifiers
are generated to avoid collisions in the DRIMSeq count table.
}
\examples{
\dontrun{
x <- read_nexons("results")
metadata <- x$get_metadata()
# Example annotations for six samples:
metadata$condition <- factor(rep(c("control", "treated"), each = 3))
metadata$batch <- factor(rep(c("A", "B", "A"), 2))
x$set_metadata(metadata)
set.seed(123)
fit <- x$run_drimseq(~ batch + condition, min_samps_feature_expr = 3)
tested <- DRIMSeq::dmTest(fit, coef = "conditiontreated")
genes <- DRIMSeq::results(tested)
isoforms <- DRIMSeq::results(tested, level = "feature")
partial_fit <- x$run_drimseq(~ batch + condition, file = "partial",
                           min_samps_feature_expr = 3)
}
}
\seealso{\code{\link{read_nexons}}, \code{\link{run_deseq2}}}
