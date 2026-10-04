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

This initial version provides import and lazy access. Conversion to DESeq2
and DRIMSeq inputs will be added separately.

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
