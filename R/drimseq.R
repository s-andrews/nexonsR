# Prepare and validate inputs independently of DRIMSeq model fitting.
prepare_drimseq <- function(object, design, file = "unique", samples = NULL,
                            min_samps_feature_expr = NULL, min_feature_expr = 10,
                            min_samps_feature_prop = min_samps_feature_expr,
                            min_feature_prop = 0.1, min_samps_gene_expr = NULL,
                            min_gene_expr = 10, run_gene_twice = FALSE,
                            filter_group = NULL) {
  validate_string(file, "file")
  if (!file %in% c("partial", "unique")) {
    stop("`file` must be 'partial' or 'unique'.", call. = FALSE)
  }
  input <- prepare_model_counts(object, file, "isoform", samples, "DRIMSeq")
  metadata <- input$metadata
  n <- nrow(metadata)
  genes <- object$get_data(file)$Gene_ID
  if (anyNA(genes) || any(!nzchar(trimws(genes)))) {
    stop("Gene_ID must not contain missing or empty identifiers.", call. = FALSE)
  }

  if (is.null(min_samps_feature_expr)) {
    min_samps_feature_expr <- drimseq_group_size(design, metadata, filter_group)
  }
  if (is.null(min_samps_feature_prop)) {
    min_samps_feature_prop <- min_samps_feature_expr
  }

  if (inherits(design, "formula")) {
    if (length(design) != 2L) {
      stop("`design` must be a one-sided formula.", call. = FALSE)
    }
    frame <- stats::model.frame(design, data = metadata, na.action = stats::na.fail)
    design <- stats::model.matrix(design, data = frame)
  }
  if (!is.matrix(design) || !is.numeric(design) || nrow(design) != n ||
      !ncol(design) || any(!is.finite(design))) {
    stop("`design` must produce a finite numeric matrix with one row per selected sample and at least one column.",
         call. = FALSE)
  }
  if (qr(design)$rank != ncol(design)) {
    stop("`design` must have full column rank.", call. = FALSE)
  }
  if (ncol(design) >= n) {
    stop("`design` must leave residual degrees of freedom.", call. = FALSE)
  }
  # Matrix designs are positional, just as in run_deseq2.
  rownames(design) <- metadata$sample
  if (is.null(colnames(design))) colnames(design) <- paste0("coef", seq_len(ncol(design)))

  if (is.null(min_samps_gene_expr)) min_samps_gene_expr <- n
  filter <- list(min_samps_feature_expr = min_samps_feature_expr,
    min_feature_expr = min_feature_expr,
    min_samps_feature_prop = min_samps_feature_prop,
    min_feature_prop = min_feature_prop,
    min_samps_gene_expr = min_samps_gene_expr, min_gene_expr = min_gene_expr,
    run_gene_twice = run_gene_twice)
  for (name in names(filter)[names(filter) != "run_gene_twice"]) {
    value <- filter[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0) {
      stop("`", name, "` must be a finite non-negative number.", call. = FALSE)
    }
    if (startsWith(name, "min_samps_") && (value != floor(value) || value > n)) {
      stop("`", name, "` must be a whole number from zero through the selected sample count.",
           call. = FALSE)
    }
  }
  if (min_feature_prop > 1) {
    stop("`min_feature_prop` must be between zero and one.", call. = FALSE)
  }
  if (!is.logical(run_gene_twice) || length(run_gene_twice) != 1L || is.na(run_gene_twice)) {
    stop("`run_gene_twice` must be a single non-missing logical value.", call. = FALSE)
  }

  # Internal count names avoid collisions with DRIMSeq's identifier columns.
  # The original sample names remain in metadata$sample for design variables.
  sample_ids <- metadata$sample
  if (any(sample_ids %in% c("gene_id", "feature_id"))) {
    sample_ids <- paste0("nexons_sample_", seq_len(n))
  }
  colnames(input$counts) <- sample_ids
  metadata$sample_id <- sample_ids
  counts <- data.frame(gene_id = genes, feature_id = rownames(input$counts),
    input$counts, check.names = FALSE, row.names = NULL)
  list(counts = counts, metadata = metadata, design = design, filter = filter)
}

# Only a bare categorical predictor identifies groups unambiguously.
drimseq_group_size <- function(design, metadata, filter_group) {
  guidance <- "Supply `filter_group` naming a categorical metadata column or an explicit `min_samps_feature_expr`."
  if (is.null(filter_group)) {
    if (inherits(design, "formula") && length(design) == 2L) {
      terms <- tryCatch(stats::terms(design), error = function(e) NULL)
      variables <- as.list(attr(terms, "variables"))[-1L]
      bare_rhs <- function(x) {
        if (is.symbol(x)) return(as.character(x) != ".")
        if (is.numeric(x)) return(length(x) == 1L && x %in% c(0, 1))
        is.call(x) && as.character(x[[1L]]) %in% c("+", "-") &&
          all(vapply(as.list(x)[-1L], bare_rhs, logical(1)))
      }
      if (bare_rhs(design[[2L]]) &&
          length(variables) == 1L && is.symbol(variables[[1L]]) &&
          length(attr(terms, "term.labels")) == 1L &&
          identical(attr(terms, "order"), 1L)) {
        filter_group <- as.character(variables[[1L]])
      }
    }
    if (is.null(filter_group)) {
      stop("Cannot infer filtering groups from this design. ", guidance, call. = FALSE)
    }
  }
  validate_string(filter_group, "filter_group")
  if (!filter_group %in% names(metadata)) {
    stop("`filter_group` must name a metadata column. ", guidance, call. = FALSE)
  }
  group <- metadata[[filter_group]]
  if (!(is.factor(group) || is.character(group) || is.logical(group)) ||
      !is.null(dim(group)) || anyNA(group)) {
    stop("Filtering groups must be categorical and have no missing values. ",
         guidance, call. = FALSE)
  }
  sizes <- table(group)
  sizes <- sizes[sizes > 0L]
  if (length(sizes) < 2L) {
    stop("Filtering requires at least two observed groups in the selected samples. ",
         guidance, call. = FALSE)
  }
  min(as.integer(sizes))
}
