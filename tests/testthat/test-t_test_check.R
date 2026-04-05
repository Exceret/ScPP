test_that("t test works", {
  set.seed(123L)
  mat <- matrix(rnbinom(100, size = 2, prob = 0.5), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  ref_pos <- c(1, 2, 6, 7, 10)
  tes_pos <- c(3, 4, 5, 8, 9)

  log2FCs <- rowMeans(mat[, tes_pos, drop = FALSE]) -
    rowMeans(mat[, ref_pos, drop = FALSE])

  # basic
  n_genes <- nrow(mat)
  pvalues <- numeric(n_genes)
  statistics <- numeric(n_genes)
  names(pvalues) <- rownames(mat)
  names(statistics) <- rownames(mat)

  for (i in seq_len(n_genes)) {
    test_result <- rlang::try_fetch(
      stats::t.test(mat[i, tes_pos], mat[i, ref_pos]),
      error = function(e) {
        cli::cli_warn(
          "t.test failed for gene {.pkg {i}}, return NA: {e$message}"
        )
        list(p.value = NA, statistic = NA)
      }
    )
    pvalues[i] <- test_result$p.value
    statistics[i] <- test_result$statistic
  }

  res1 <- dplyr::bind_cols(pvalues = pvalues, statistics = statistics)

  # matrixTest
  t_results <- matrixTests::row_t_welch(
    mat[, tes_pos, drop = FALSE],
    mat[, ref_pos, drop = FALSE]
  )
  pvalues <- t_results$pvalue
  statistics <- t_results$statistic
  names(pvalues) <- rownames(mat)
  names(statistics) <- rownames(mat)

  res2 <- dplyr::bind_cols(pvalues = pvalues, statistics = statistics)

  expect_lt(
    max(abs(c(res1$pvalues - res2$pvalues, res1$statistics - res2$statistics))),
    1e-6
  )

  # genefilter
  # ! result doesn't match
  #   grp <- factor(c(rep("tes", length(tes_pos)), rep("ref", length(ref_pos))))
  #   combined_mat <- mat[, c(tes_pos, ref_pos), drop = FALSE]

  #   res <- genefilter::rowttests(combined_mat, grp)

  #   pvalues <- setNames(res$p.value, rownames(mat))
  #   statistics <- setNames(res$statistic, rownames(mat))

  #     res3 <- dplyr::bind_cols(pvalues = pvalues, statistics = statistics)

  # viper
  t_results <- viper::rowTtest(
    mat[, tes_pos, drop = FALSE],
    mat[, ref_pos, drop = FALSE]
  )
  pvalues <- t_results$p.value
  statistics <- t_results$statistic
  names(pvalues) <- rownames(mat)
  names(statistics) <- rownames(mat)

  res4 <- dplyr::bind_cols(pvalues = pvalues, statistics = statistics)
})

test_that("t test better", {
  set.seed(123L)
  mat <- matrix(rnbinom(100, size = 2, prob = 0.5), nrow = 10)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("sample", 1:10)

  ref_pos <- c(1, 2, 6, 7, 10)
  tes_pos <- c(3, 4, 5, 8, 9)

  log2FCs <- rowMeans(mat[, tes_pos, drop = FALSE]) -
    rowMeans(mat[, ref_pos, drop = FALSE])

  time <- microbenchmark::microbenchmark(
    basic = {
      n_genes <- nrow(mat)
      pvalues <- numeric(n_genes)
      statistics <- numeric(n_genes)
      names(pvalues) <- rownames(mat)
      names(statistics) <- rownames(mat)

      for (i in seq_len(n_genes)) {
        test_result <-
          stats::t.test(mat[i, tes_pos], mat[i, ref_pos])

        pvalues[i] <- test_result$p.value
        statistics[i] <- test_result$statistic
      }
    },
    matrixTest = {
      t_results <- matrixTests::row_t_welch(
        mat[, tes_pos, drop = FALSE],
        mat[, ref_pos, drop = FALSE]
      )
      pvalues <- t_results$pvalue
      statistics <- t_results$statistic
      names(pvalues) <- rownames(mat)
      names(statistics) <- rownames(mat)
    },
    viper = {
      t_results <- viper::rowTtest(
        mat[, tes_pos, drop = FALSE],
        mat[, ref_pos, drop = FALSE]
      )
      pvalues <- t_results$p.value
      statistics <- t_results$statistic
      names(pvalues) <- rownames(mat)
      names(statistics) <- rownames(mat)
    }
  )
})
