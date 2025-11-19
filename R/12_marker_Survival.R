#' @title Find Survival-Associated Markers
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param survival_data Survival data with time in column1 and status in column2. Rownames are sample name.
#'
#' @return A gene list of survival-associated markers and ranked genes based on Cox hazard ratios.
#' @export
#' @family scPP
#'
marker_Survival2 <- function(bulk_data, survival_data) {
    # # 输入验证 (上游已检查)
    # if (missing(bulk_data) || !inherits(bulk_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg bulk_data} is missing or incorrect.")
    # }
    # if (missing(survival_data) || !inherits(survival_data, c("matrix", "data.frame"))) {
    #     cli::cli_abort("{.arg survival_data} is missing or incorrect.")
    # }

    SurvivalData <- data.frame(cbind(survival_data, Matrix::t(bulk_data)))
    colnames(SurvivalData) = make.names(colnames(SurvivalData))
    var <- make.names(rownames(bulk_data))

    Model_Formula <- sapply(var, function(x) {
        stats::as.formula(paste("survival::Surv(time, status) ~", x))
    })

    Model_all <- lapply(Model_Formula, function(x) {
        survival::coxph(x, data = SurvivalData)
    })

    res <- lapply(seq_along(Model_all), function(i) {
        coef_summary <- Matrix::summary(Model_all[[i]])$coefficients
        data.frame(
            variable = var[i],
            pvalue = coef_summary[, 5],
            coef = coef_summary[, 2]
        )
    }) %>%
        dplyr::bind_rows()

    genes_sort <- res %>%
        dplyr::arrange(dplyr::desc(coef)) %>%
        dplyr::pull(coef, name = variable)

    res <- res[order(res$pvalue), ]
    res$fdr <- stats::p.adjust(res$pvalue, method = "fdr")

    gene_pos <- res %>%
        dplyr::filter(fdr < 0.05, coef > 1) %>%
        dplyr::pull(variable) # correalted with worse survival
    gene_neg <- res %>%
        dplyr::filter(fdr < 0.05, coef < 1) %>%
        dplyr::pull(variable) # correlated with better survival

    geneList <- list(
        gene_pos = gene_pos,
        gene_neg = gene_neg,
        genes_sort = genes_sort
    )

    # ? It's confusing but it is
    if (length(gene_pos) > 0 & length(gene_neg) > 0) {
        return(geneList)
    } else if (length(gene_pos) == 0) {
        cli::cli_warn(
            "There are no genes negatively correlated with patients' prognosis in this bulk dataset."
        )
        return(list(gene_pos = gene_neg))
    } else if (length(gene_neg) == 0) {
        cli::cli_warn(
            "There are no genes positively correlated with patients' prognosis in this bulk dataset."
        )
        return(list(gene_neg = gene_pos))
    }
}
