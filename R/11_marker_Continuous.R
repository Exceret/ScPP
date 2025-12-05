#' @title Find Continuous Feature Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, such as TMB or CNA values of each sample.
#' @param estimate_cutoff Absolute cutoff value of correlation coefficient, default is 0.2.
#' @param method Method uses for cor.test, default is "spearman", another choice is "pearson".
#'
#' @return A gene list of feature markers and ranked genes based on correlation coefficients.
#' @export
#' @family scPP
#'
marker_Continuous.optimized <- function(
    bulk_data,
    features,
    method = "spearman",
    estimate_cutoff = 0.2
) {
    # # 输入验证 (上游已检查)
    # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
    # }
    # if (missing(features) || !is.numeric(features)) {
    #     cli::cli_abort("{.arg features} is missing or incorrect.")
    # }

    if (!is.matrix(bulk_data)) {
        bulk_data <- as.matrix(bulk_data)
    }

    features_log <- log2(as.numeric(features) + 1)
    n_genes <- nrow(bulk_data)

    pvalues <- numeric(n_genes)
    estimates <- numeric(n_genes)
    names(pvalues) <- rownames(bulk_data)
    names(estimates) <- rownames(bulk_data)

    for (i in seq_len(n_genes)) {
        cor_result <- rlang::try_fetch(
            stats::cor.test(
                bulk_data[i, ],
                features_log,
                method = method
            ),
            error = function(e) list(p.value = NA, estimate = NA)
        )
        pvalues[i] <- cor_result$p.value
        estimates[i] <- cor_result$estimate
    }

    genes_sort <- sort(estimates[!is.na(estimates)], decreasing = TRUE)

    res <- data.table::data.table(
        gene = names(pvalues),
        pvalue = pvalues,
        estimate = estimates
    )

    data.table::setorder(res, pvalue)
    res[, fdr := stats::p.adjust(pvalue, method = "fdr")]

    gene_pos <- res[fdr < 0.05 & estimate > estimate_cutoff, gene]
    gene_neg <- res[fdr < 0.05 & estimate < -estimate_cutoff, gene]

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
