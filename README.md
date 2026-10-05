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

The package provides import, lazy access, and DESeq2 model fitting.
DRIMSeq support will be added separately.

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
