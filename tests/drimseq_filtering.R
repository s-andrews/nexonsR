library(nexonsR)

local({
  prepare <- getFromNamespace("prepare_drimseq", "nexonsR")
  folder <- tempfile()
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  samples <- paste0("s", 1:8)
  counts <- rbind(c(rep(100, 3), rep(0, 5)), rep(100, 8))
  colnames(counts) <- samples
  data <- data.frame(Transcript_ID = c("specific", "shared"), Gene_ID = "gene",
    Gene_Name = "gene", Chr = "1", Start = 1, End = 10, Strand = "+", counts)
  for (file in c("partial", "unique")) {
    utils::write.table(data, file.path(folder, paste0(file, ".txt")),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  object <- suppressWarnings(read_nexons(folder, ""))
  metadata <- object$metadata
  metadata$condition <- factor(rep(c("A", "B"), c(3, 5)), levels = c("A", "B", "unused"))
  metadata$batch <- factor(rep(c("X", "Y"), 4))
  metadata$age <- 1:8
  metadata$character_group <- as.character(metadata$condition)
  metadata$logical_group <- metadata$condition == "A"
  metadata$three_groups <- factor(rep(c("A", "B", "C"), c(2, 3, 3)))
  metadata$missing_group <- c(NA_character_, rep("B", 7))
  metadata$one_group <- "A"
  object$set_metadata(metadata)
  error <- function(expr, pattern) {
    result <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"),
              grepl(pattern, conditionMessage(result), fixed = TRUE))
  }
  for (design in list(~condition, ~0 + condition, ~condition - 1,
                     ~character_group, ~logical_group)) {
    input <- prepare(object, design)
    stopifnot(input$filter$min_samps_feature_expr == 3,
              input$filter$min_samps_feature_prop == 3,
              input$filter$min_samps_gene_expr == 8)
  }
  stopifnot(prepare(object, ~three_groups)$filter$min_samps_feature_expr == 2)
  selected <- samples[c(8, 1, 4, 2)]
  subset <- prepare(object, ~condition, samples = selected)
  stopifnot(subset$filter$min_samps_feature_expr == 2,
            subset$filter$min_samps_feature_prop == 2,
            subset$filter$min_samps_gene_expr == 4,
            identical(subset$metadata$sample, selected))
  stopifnot(prepare(object, ~condition, samples = samples[c(1, 4:8)])$filter$min_samps_feature_expr == 1)
  # Drop the unused condition level before building a full-rank matrix.
  matrix_design <- stats::model.matrix(~batch + condition, droplevels(metadata))
  for (design in list(~batch + condition, ~batch * condition, ~age,
                     ~factor(condition), ~I(age), ~condition^2, ~1, ~., matrix_design)) {
    error(prepare(object, design), "filter_group")
  }
  for (design in list(~batch + condition, ~batch * condition, ~age, matrix_design)) {
    for (group in c("condition", "character_group", "logical_group")) {
      input <- prepare(object, design, filter_group = group)
      stopifnot(input$filter$min_samps_feature_expr == 3)
    }
  }
  for (group in list("unknown", "", NA_character_, 1, c("condition", "batch"))) {
    error(prepare(object, ~condition, filter_group = group), "filter_group")
  }
  for (group in c("age", "missing_group", "one_group")) {
    error(prepare(object, ~condition, filter_group = group), "filter_group")
  }
  error(prepare(object, ~1, samples = samples[1:3], filter_group = "condition"), "two observed groups")
  error(prepare(object, ~age), "categorical")
  stopifnot(prepare(object, ~age, min_samps_feature_expr = 2,
                    filter_group = "unknown")$filter$min_samps_feature_prop == 2,
            prepare(object, ~condition, min_samps_feature_prop = 1)$filter$min_samps_feature_prop == 1,
            prepare(object, ~condition, min_samps_feature_prop = NULL)$filter$min_samps_feature_prop == 3,
            prepare(object, ~age, min_samps_feature_expr = 2,
                    min_samps_feature_prop = NULL)$filter$min_samps_feature_prop == 2)
  apply_filter <- function(input) {
    dataset <- DRIMSeq::dmDSdata(input$counts, input$metadata)
    do.call(DRIMSeq::dmFilter, c(list(x = dataset), input$filter))
  }
  for (file in c("partial", "unique")) {
    inferred <- prepare(object, ~condition, file = file)
    explicit <- prepare(object, ~condition, file = file, min_samps_feature_expr = 3)
    stopifnot(isTRUE(all.equal(inferred, explicit)))
    retained <- DRIMSeq::counts(apply_filter(inferred))
    stopifnot(setequal(as.character(retained$feature_id), c("specific", "shared")),
              identical(retained, DRIMSeq::counts(apply_filter(explicit))))
  }
  stopifnot(identical(object$metadata, metadata))
  # Exercise inference through the public method without fitting a tiny dataset.
  error(object$run_drimseq(~condition, min_gene_expr = 1e9), "no analysable genes")
  error(object$run_drimseq(~batch + condition, filter_group = "condition",
                          min_gene_expr = 1e9), "no analysable genes")
})
