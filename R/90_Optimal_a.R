# ? This function was integrated into ScPP()

# #' Title
# #'
# #' @param sc_dataset A seurat object of single cell RNA sequencing data.
# #' @param geneList A gene list correlated with interested features.
# #' @param probs Cutoff value of ScPP, default is seq(0.2, 0.45, by = 0.05).
# #'
# #' @return The optimal prob based on NES.
# #' @export
# #'
# #' @examples
# optimal_a = function(sc_dataset, geneList, probs = seq(0.2, 0.45, by = 0.05)) {
#     set.seed(123)
#     geneList_AUC = geneList[names(geneList) != "genes_sort"]
#     genes_sort = geneList$genes_sort

#     if (length(geneList_AUC) != 2) {
#         stop(
#             "This gene list do not have enough information correlated with interested feature."
#         )
#     }

#     if (missing(sc_dataset) || class(sc_dataset) != "Seurat") {
#         stop("'sc_dataset' is missing or not a seurat object.")
#     }

#     library(AUCell)
#     seurat_version <- packageVersion("Seurat")
#     if (!is.na(seurat_version) && seurat_version >= "5.0.0") {
#         cellrankings = AUCell_buildRankings(
#             sc_dataset@assays$RNA$data,
#             plotStats = FALSE
#         )
#     } else {
#         cellrankings = AUCell_buildRankings(
#             sc_dataset@assays$RNA@data,
#             plotStats = FALSE
#         )
#     }
#     cellAUC = AUCell_calcAUC(geneList_AUC, cellrankings)

#     metadata = as.data.frame(sc_dataset@meta.data)
#     metadata$AUCup <- as.numeric(getAUC(cellAUC)["gene_pos", ])
#     metadata$AUCdown <- as.numeric(getAUC(cellAUC)["gene_neg", ])

#     NES_dif_res = c()
#     for (i in 1:length(probs)) {
#         downcells1 = rownames(metadata)[which(
#             metadata$AUCup <= quantile(metadata$AUCup, probs = probs[i])
#         )]
#         upcells1 = rownames(metadata)[which(
#             metadata$AUCup >= quantile(metadata$AUCup, probs = (1 - probs[i]))
#         )]
#         downcells2 = rownames(metadata)[which(
#             metadata$AUCdown >=
#                 quantile(metadata$AUCdown, probs = (1 - probs[i]))
#         )]
#         upcells2 = rownames(metadata)[which(
#             metadata$AUCdown <= quantile(metadata$AUCdown, probs = probs[i])
#         )]

#         ScPP_neg = intersect(downcells1, downcells2)
#         ScPP_pos = intersect(upcells1, upcells2)

#         metadata$ScPP <- ifelse(
#             rownames(metadata) %in% ScPP_pos,
#             "Phenotype+",
#             "Background"
#         )
#         metadata$ScPP <- ifelse(
#             rownames(metadata) %in% ScPP_neg,
#             "Phenotype-",
#             metadata$ScPP
#         )

#         sc_dataset$ScPP = metadata$ScPP
#         Idents(sc_dataset) = "ScPP"

#         markers <- FindMarkers(
#             sc_dataset,
#             ident.1 = "Phenotype+",
#             ident.2 = "Phenotype-"
#         )

#         genes_pos <- rownames(markers[
#             which(markers$avg_log2FC > 1 & markers$p_val_adj < 0.05),
#         ])
#         genes_neg <- rownames(markers[
#             which(markers$avg_log2FC < -1 & markers$p_val_adj < 0.05),
#         ])

#         if (length(genes_pos) == 0 | length(genes_neg) == 0) {
#             NES_dif = cbind(probs[i], NA)
#             NES_dif_res = rbind(NES_dif_res, NES_dif)
#         }

#         if (length(genes_pos) > 0 & length(genes_neg) > 0) {
#             res <- list(Genes_pos = genes_pos, Genes_neg = genes_neg)

#             library(fgsea)
#             fgseaRes <- fgsea(pathways = res, stats = genes_sort)
#             NES_dif = fgseaRes$NES[fgseaRes$pathway == "Genes_pos"] -
#                 fgseaRes$NES[fgseaRes$pathway == "Genes_neg"]
#             NES_dif = cbind(probs[i], NES_dif)
#             NES_dif_res = rbind(NES_dif_res, NES_dif)
#         }
#     }
#     colnames(NES_dif_res) = c("prob", "NES_dif")
#     NES_dif_res = as.data.frame(NES_dif_res)
#     opt = NES_dif_res$prob[which.max(na.omit(NES_dif_res$NES_dif))]
#     return(opt)
# }
