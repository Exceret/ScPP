#' @title Find Binary Feature Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, column1 are sample names (colname is "Sample") and column2 are feature (colname is "Feature") labels of each sample.
#' @param ref_group A character to indicate which feature is the control group.
#' @param Log2FC_cutoff Absolute cutoff value of fold change, default is 0.585.
#'
#' @return A gene list of feature markers with two binary groups and ranked genes based on t-test's statistic.
#' @export
#' @family scPP
#'
marker_Binary.optimized <- function(
  bulk_data,
  features,
  ref_group,
  Log2FC_cutoff = 0.585
) {
  # if (missing(ref_group)) {
  #     cli::cli_abort("{.arg ref_group }is missing or incorrect.")
  # }
  # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
  #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
  # }
  # if (missing(features) || !inherits(features, c("matrix", "data.frame"))) {
  #     cli::cli_abort("{.arg features} is missing or incorrect.")
  # }

  if (!is.matrix(bulk_data)) {
    bulk_data <- as.matrix(bulk_data)
  }

  features_dt <- data.table::as.data.table(features)
  ref <- features_dt[Feature == ref_group, Sample]
  tes <- features_dt[Feature != ref_group, Sample]

  ref_pos <- which(colnames(bulk_data) %chin% ref)
  tes_pos <- which(colnames(bulk_data) %chin% tes)

  log2FCs <- rowMeans(bulk_data[, tes_pos, drop = FALSE]) -
    rowMeans(bulk_data[, ref_pos, drop = FALSE])

  if (requireNamespace("matrixTests", quietly = TRUE)) {
    row_t_welch <- getExportedValue("matrixTests", "row_t_welch")
    t_results <- row_t_welch(
      bulk_data[, tes_pos, drop = FALSE],
      bulk_data[, ref_pos, drop = FALSE]
    )
    pvalues <- t_results$pvalue
    statistics <- t_results$statistic
    names(pvalues) <- rownames(bulk_data)
    names(statistics) <- rownames(bulk_data)
  } else {
    n_genes <- nrow(bulk_data)
    pvalues <- numeric(n_genes)
    statistics <- numeric(n_genes)
    names(pvalues) <- rownames(bulk_data)
    names(statistics) <- rownames(bulk_data)

    for (i in seq_len(n_genes)) {
      test_result <- rlang::try_fetch(
        stats::t.test(bulk_data[i, tes_pos], bulk_data[i, ref_pos]),
        error = function(e) list(p.value = NA, statistic = NA)
      )
      pvalues[i] <- test_result$p.value
      statistics[i] <- test_result$statistic
    }
  }

  genes_sort <- sort(statistics[!is.na(statistics)], decreasing = TRUE)

  res <- data.table::data.table(
    gene = names(pvalues),
    pvalue = pvalues,
    log2FC = log2FCs
  )

  res[, fdr := stats::p.adjust(pvalue, method = "fdr")]

  gene_pos <- res[pvalue < 0.05 & log2FC > Log2FC_cutoff, gene]
  gene_neg <- res[pvalue < 0.05 & log2FC < -Log2FC_cutoff, gene]

  geneList <- list(
    gene_pos = gene_pos,
    gene_neg = gene_neg,
    genes_sort = genes_sort
  )

  has_pos <- length(gene_pos) > 0
  has_neg <- length(gene_neg) > 0

  if (has_pos && has_neg) {
    return(geneList)
  } else if (!has_pos) {
    cli::cli_warn(
      "There are no genes positively correlated with the given feature in this bulk dataset."
    )
    return(list(gene_neg = gene_neg))
  } else {
    cli::cli_warn(
      "There are no genes negatively correlated with the given feature in this bulk dataset."
    )
    return(list(gene_pos = gene_pos))
  }
}
