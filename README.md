# nexonsR
An R package to analyse nanopore RNA-Seq quantitated with nexons

```r
library(nexonsR)
results <- read_nexons("path/to/results", prefix = "nexons_output_")
results             # Shows file availability without reading tables
results$files       # Absolute paths; missing files have NA paths
genes <- results$get_data("gene") # Reads and caches only the gene table
isoforms <- results$get_data("partial")
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

The package provides import, lazy access, and DESeq2 and DRIMSeq model fitting.

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

Fit a DESeq2 model after adding the design variables to metadata:

```r
# Install the optional dependency first: BiocManager::install("DESeq2")
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
# Install the optional dependency first: BiocManager::install("DRIMSeq")
# Add condition and batch columns to the object's metadata first.
set.seed(123)
fit <- results$run_drimseq(~ batch + condition, min_samps_feature_expr = 3)
tested <- DRIMSeq::dmTest(fit, coef = "conditiontreated")
DRIMSeq::results(tested)                    # Gene-level usage changes
DRIMSeq::results(tested, level = "feature") # Individual isoforms

partial_fit <- results$run_drimseq(~ batch + condition, file = "partial",
                                  min_samps_feature_expr = 3)
```

Both `partial` and `unique` (the default) contain uniquely assigned reads and
receive identical handling. Supply a formula or a numeric design matrix whose
rows match the selected sample order. The method returns a fitted `dmDSfit`;
choose coefficients, contrasts or a reduced design in `DRIMSeq::dmTest()`.

`min_samps_feature_expr` must be explicit; use the smallest group size for a
simple group comparison. Defaults require isoform counts of at least 10 and
within-gene proportions of at least 0.1 in that many samples, plus summed gene
counts of at least 10 in all selected samples. Lower `min_feature_prop` (or set
it to zero) to retain minor isoforms; lower `min_samps_gene_expr` to tolerate
poorly observed samples in larger cohorts. These are starting thresholds, not
nanopore-specific calibration. See `?run_drimseq` for all filtering arguments.
