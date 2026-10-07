# nexonsR
An R package to analyse nanopore RNA-Seq quantitated with nexons

```r
library(nexonsR)
results <- read_nexons("path/to/results", prefix = "nexons_output_")
results             # Shows file availability without reading tables
results$files       # Absolute paths; missing files have NA paths
genes <- results$get_data("gene") # Reads and caches only the gene table
isoforms <- results$get_data("partial")
selected <- results$get_data("gene", samples = c("sample-2", "sample-1"))
log2rpm <- results$get_data("unique", units = "log2RPM")
results$clear_cache()
```

The reader uses an R6 object and looks for the literal prefix followed by
`gene.txt`, `partial.txt`, and `unique.txt`. Thus the
default gene filename is `nexons_output_gene.txt`. Use `prefix = ""`
for files named simply `gene.txt`, etc.

`partial` counts include complete and partial reads; `unique` counts include
only complete reads. Other files are ignored. Missing result files generate
a warning; finding none generates an error. Discovery reads only the first line
of each available file and checks that all files have the same unique sample
names, regardless of order. Count rows are read only when accessed. Sample names and
annotation identifiers are preserved. Tables remain cached until cleared;
source files must remain available for uncached reads.

`get_data()` returns raw counts by default. Its `samples` argument returns only
the named samples in the requested order while retaining every annotation
column. With `units = "log2RPM"`, each selected sample is converted using
`log2(count / sum(counts) * 1e6 + 1)`. Selection and conversion do not change
the cached raw table.

The package provides import, lazy access, and DESeq2 and DRIMSeq model fitting.
DESeq2 and DRIMSeq are required dependencies. When installing from a repository,
configure both CRAN and Bioconductor repositories so the installer can resolve
them. Local tarball installs with `repos = NULL` require dependencies to be
installed beforehand, for example with
`BiocManager::install(c("DESeq2", "DRIMSeq"))`.

Sample annotations are stored in `metadata`, initially with a single `sample`
column. Access and update them using:

```r
metadata <- results$get_metadata() # Also available as results$metadata
metadata$condition <- rep("control", nrow(metadata))
results$set_metadata(metadata)
```

Updates must retain each original sample exactly once. You can reorder rows;
the setter restores the order from the first available file (gene, partial,
then unique), keeping annotations attached to their samples. Additional columns
are preserved. Count tables retain their own sample column order; use sample
names when matching them to metadata.

Extract gene and transcript annotations as data frames:

```r
genes <- results$gene_metadata() # Gene_ID and Gene_Name from the gene file
genes <- results$gene_metadata(file = "partial") # Distinct ID/name pairs
transcripts <- results$transcript_metadata() # Defaults to the unique file
transcripts <- results$transcript_metadata(file = "partial")
```

Column names retain their original spelling: `Gene_ID`, `Gene_Name`, and
`Transcript_ID`. Gene metadata preserves all rows from the gene file and
returns distinct pairs from partial/unique files in first-occurrence order.
Transcript metadata returns `Transcript_ID`, `Gene_ID`, and `Gene_Name` in
that order, preserving all rows; `file = "gene"` raises an error.

Count observed genes or isoforms per sample:

```r
results$coverage()
results$coverage(file = "partial", threshold = 10)
results$coverage(file = "unique", level = "isoform")
```

The returned data frame contains all sample metadata columns followed by
`gene` or `isoform`, preserving metadata order and column types. Existing
metadata is unchanged. If the count column name already exists, the new column
receives a numeric suffix (such as `gene.1`).
Counts at or above the positive threshold (default 1)
are observed. Gene coverage from partial/unique files sums counts by `Gene_ID`
before applying the threshold. Isoform coverage requires partial/unique data.
The selected table is loaded and cached as needed.

Compare quantitations across assignment files:

```r
results$compare_assignments()
results$compare_assignments(level = "isoform", samples = "sample-1",
                            units = "log2RPM")
```

Gene output contains `Gene_ID`, `Number_of_Isoforms`, `sample`, and available
`gene`, `partial`, and `unique` quantitations. Isoform counts use the larger
count from partial and unique, or `NA` when neither provides the gene.
Isoform output contains `Transcript_ID`, `Gene_ID`, `sample`, `partial`, and
`unique` and requires both isoform files. Raw isoform counts are summed by gene
before log2RPM conversion, using each file's own sample totals. All features
are retained; differing feature sets produce one warning and missing values
are `NA`. Samples follow the supplied order, or metadata order by default.

Fit a DESeq2 model after adding the design variables to metadata:

```r
model <- results$run_deseq2(~ condition)
DESeq2::results(model)

# Sum isoform counts by Gene_ID for gene-level analysis:
model <- results$run_deseq2(~ condition, file = "partial", level = "gene")

# Analyse individual isoforms in selected samples:
model <- results$run_deseq2(~ condition, file = "unique", level = "isoform",
                           samples = c("control1", "control2", "treated1", "treated2"))
```

`design` is the only required argument. The method aligns counts and metadata
by sample name and returns the fitted model without storing it in the nexonsR
object. Gene files cannot be used for isoform analysis. Counts must be
non-negative whole numbers; fractional counts are rejected without rounding.

Fit differential isoform usage models with DRIMSeq:

```r
# Add condition and batch columns to the object's metadata first.
set.seed(123)
fit <- results$run_drimseq(~ batch + condition, filter_group = "condition")
tested <- DRIMSeq::dmTest(fit, coef = "conditiontreated")
DRIMSeq::results(tested)                    # Gene-level usage changes
DRIMSeq::results(tested, level = "feature") # Individual isoforms

partial_fit <- results$run_drimseq(~ batch + condition, file = "partial",
                                  filter_group = "condition")

# Infer the smallest group size from a simple categorical design:
simple_fit <- results$run_drimseq(~ condition)
# Or supply an explicit threshold:
explicit_fit <- results$run_drimseq(~ batch + condition, min_samps_feature_expr = 3)
```

Both `partial` and `unique` (the default) contain uniquely assigned reads and
receive identical handling. Supply a formula or a numeric design matrix whose
rows match the selected sample order. The method returns a fitted `dmDSfit`;
choose coefficients, contrasts or a reduced design in `DRIMSeq::dmTest()`.

`min_samps_feature_expr = NULL` uses the smallest observed group size after
sample selection, ignoring unused factor levels. A formula with one bare
categorical predictor identifies groups automatically. For more complex
formulas or matrix designs, supply `filter_group` naming a categorical metadata
column, or an explicit numeric threshold. Filtering groups must have no missing
values and at least two observed groups. `filter_group` controls filtering only;
it does not alter the model and is ignored with an explicit threshold.
Omitted or `NULL` `min_samps_feature_prop` inherits the resolved sample threshold.
This permits isoforms expressed only in one condition to survive filtering.

Defaults require isoform counts of at least 10 and
within-gene proportions of at least 0.1 in that many samples, plus summed gene
counts of at least 10 in all selected samples. Lower `min_feature_prop` (or set
it to zero) to retain minor isoforms; lower `min_samps_gene_expr` to tolerate
poorly observed samples in larger cohorts. These are starting thresholds, not
nanopore-specific calibration. See `?run_drimseq` for all filtering arguments.
