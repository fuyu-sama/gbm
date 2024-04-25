#! /usr/bin/env Rscript

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#  Codes are far away from bugs with Buddha's bless
#

# %% environment config
library <- function(...) suppressMessages(base::library(...))
library('parallel')
library('dplyr')
library('xlsx')
library('ggplot2')
library('ggvenn')
library('ComplexHeatmap')
library('Seurat')
library('clusterProfiler')
library('org.Hs.eg.db')

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
source(fs::path(WORKDIR, "src", "R", "utils.R"))

# %% read data
readData <- function(idx) {
    read.dir <- fs::path(WORKDIR, "spaceranger", idx, "outs")
    filename <- "filtered_feature_bc_matrix.h5"
    if (idx == "21B-603-5") filename <- "raw_feature_bc_matrix.h5"
    image <- Read10X_Image(fs::path(read.dir, "spatial"), filter.matrix = FALSE)
    seurat.obj <- Load10X_Spatial(
        read.dir,
        filename = filename,
        filter.matrix = FALSE,
        image = image
    )
    names(seurat.obj@images) <- idx
    return(seurat.obj)
}
seurat.list <- mclapply(idx.list, readData)

# %% preprocessing
preProcessing <- function(seurat.obj) {
    seurat.obj <- SCTransform(
        seurat.obj,
        assay = "Spatial",
        variable.features.n = 9000,
        return.only.var.genes = FALSE,
        verbose = FALSE
    )
    seurat.obj <- RunPCA(
        seurat.obj, assay = "SCT", verbose = FALSE)
    seurat.obj <- FindNeighbors(
        seurat.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    seurat.obj <- RunUMAP(
        seurat.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    seurat.obj <- RunTSNE(
        seurat.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, preProcessing)

# %% batch clustering
batchCluster <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    save.dir <- fs::path(save.dirs[[idx]], "cluster")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    for (resolution in seq(0.06, 0.5, 0.01)) {
        seurat.obj <- FindClusters(
            seurat.obj, verbose = FALSE, resolution = resolution)
        p1 <- DimPlot(
            seurat.obj, reduction = "tsne", label = TRUE)
        p2 <- SpatialDimPlot(
            seurat.obj, label = TRUE, label.size = 3, alpha = 0.6)
        save.path <- fs::path(save.dir, paste0(resolution, ".cluster.pdf"))
        ggsave(save.path, p1 + p2, width = 14, height = 7)
    }
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, batchCluster)

# %% final cluster
finalCluster <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    resolutions <- list(
        "21B-603-5" = 0.15,
        "22F-10823-3" = 0.18,
        #"22F-21576-1" = 0.14,
        #"22F-23738-2" = 0.11
        "22F-21576-1" = 0.2,
        "22F-23738-2" = 0.2
    )
    seurat.obj <- FindClusters(
        seurat.obj, verbose = FALSE, resolution = resolutions[[idx]])
    p1 <- DimPlot(
        seurat.obj, reduction = "tsne", label = TRUE)
    p2 <- SpatialDimPlot(
        seurat.obj, label = TRUE, label.size = 3, alpha = 0.6)
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".cluster.pdf"))
    ggsave(save.path, p1 + p2, width = 14, height = 7)
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, finalCluster)
saveSeuratList()

# %% save expression matrix
saveExpression <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    DefaultAssay(seurat.obj) <- "Spatial"
    expression.df <- GetAssayData(
        seurat.obj, slot = "count", assay = "Spatial"
    ) %>% t()
    write.csv(
        expression.df,
        fs::path(WORKDIR, "Data", "counts", paste0(idx, ".csv"))
    )
}
mclapply(seurat.list, saveExpression)

# %% annotation
regionAnnotation <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    region.annotations <- list()
    region.annotations[["21B-603-5"]] <- c(
        "4" = "Blood vessel rich tumor area",
        "2" = "Tumor cell densely populated area",
        "0" = "Tumor area 1",
        "1" = "Junction area",
        "3" = "Parancerous area"
    )
    region.annotations[["22F-10823-3"]] <- c(
        "2" = "Blood vessel rich tumor area",
        #"0" = "Tumor cell densely populated area",
        #"5" = "Tumor area 2",
        "5" = "Tumor cell densely populated area",
        "0" = "Tumor area 2",
        "3" = "Junction area 2",
        "1" = "Junction area",
        "4" = "Parancerous area"
    )
    region.annotations[["22F-21576-1"]] <- c(
        "0" = "Blood vessel rich tumor area",
        "2" = "Tumor cell densely populated area",
        "3" = "Tumor area 3",
        "1" = "Tumor area 4"
    )
    region.annotations[["22F-23738-2"]] <- c(
        "2" = "Blood vessel rich tumor area",
        "0" = "Tumor cell densely populated area",
        "1" = "Tumor area 5",
        "4" = "Tumor area 6",
        "3" = "Tumor area 7",
        "5" = "Tumor area 8"
    )
    seurat.obj <- RenameIdents(seurat.obj, region.annotations[[idx]])
    p1 <- DimPlot(seurat.obj, reduction = "tsne", label = FALSE) +
        theme(legend.text = element_text(size = 20))
    p2 <- SpatialDimPlot(seurat.obj, alpha = 0.6, label = FALSE) +
        theme(legend.text = element_text(size = 20))
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区.pdf"))
    ggsave(save.path, p1 + p2, width = 25, height = 7)
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, regionAnnotation)

# %% DE 1
differentialExpression1 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    markers <- FindAllMarkers(
        seurat.obj, logfc.threshold = 0.1, verbose = FALSE
    )
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区域差异表达.csv"))
    write.csv(markers, save.path)

    tops <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
    spots <- names(sort(Idents(seurat.obj)))
    draw.matrix <- GetAssayData(
        seurat.obj, assay = "SCT", slot = "scale.data")[tops$gene, spots]
    save.path <- fs::path(
        save.dirs[[idx]], paste0(idx, ".分区域差异表达热图.pdf")
    )
    p <- DoHeatmap(seurat.obj, features = tops$gene) + NoLegend()
    ggsave(save.path, p, height = 20, width = 15)

    return(markers)
}
markers.list1 <- mclapply(seurat.list, differentialExpression1)

# %% GO 1
enrichment1 <- function(idx) {
    clusters <- unique(markers.list1[[idx]]$cluster)
    for (cluster in clusters) {
        upscale.table <- markers.list1[[idx]] %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC > 0)
        downscale.table <- markers.list1[[idx]] %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC < 0)
        upscale.genes <- upscale.table[upscale.table$cluster == cluster, ]$gene
        upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
        downscale.genes <- downscale.table[downscale.table$cluster == cluster, ]$gene
        downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
        if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
            print(paste(idx, cluster))
            next
        }

        save.dir <- fs::path(save.dirs[[idx]], "GO")
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
        upscale.ego <- enrichGO(
            upscale.genes,
            OrgDb = OrgDb,
            keyType = "SYMBOL",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        downscale.ego <- enrichGO(
            downscale.genes,
            OrgDb = OrgDb,
            keyType = "SYMBOL",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        p1 <- dotplot(upscale.ego, split = "ONTOLOGY") +
            facet_grid(ONTOLOGY~., scale = "free") +
            ggtitle("GO Enrichment for log2FC > 0 genes")
        p2 <- dotplot(downscale.ego, split = "ONTOLOGY") +
            facet_grid(ONTOLOGY~., scale = "free") +
            ggtitle("GO Enrichment for log2FC < 0 genes")
        cluster <- stringr::str_replace_all(cluster, " ", "_")
        ggsave(
            fs::path(save.dir, paste(idx, cluster, "GO.pdf", sep = ".")),
            p1 + p2,
            width = 14,
            height = 15
        )
        write.xlsx2(
            as.data.frame(upscale.ego),
            fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC > 0"
        )
        write.xlsx2(
            as.data.frame(downscale.ego),
            fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )

        save.dir <- fs::path(save.dirs[[idx]], "KEGG")
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
        upscale.ids <- mapIds(
            OrgDb, keys = upscale.genes, column = "ENTREZID", keytype = "SYMBOL")
        downscale.ids <- mapIds(
            OrgDb, keys = downscale.genes, column = "ENTREZID", keytype = "SYMBOL")
        upscale.kk <- enrichKEGG(
            upscale.ids,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        downscale.kk <- enrichKEGG(
            downscale.ids,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        if (dim(upscale.kk)[1] > 0) {
            p1 <- dotplot(upscale.kk) +
                ggtitle("KEGG Enrichment for log2FC > 0 genes")
            write.xlsx2(
                as.data.frame(upscale.kk),
                fs::path(
                    save.dir,
                    paste(idx, cluster, "KEGG.xlsx", sep = ".")
                    ),
                sheetName = "avg_log2FC > 0"
            )
        } else {
            p1 <- NULL
        }
        if (dim(downscale.kk)[1] > 0) {
            p2 <- dotplot(downscale.kk) +
                ggtitle("KEGG Enrichment for log2FC < 0 genes")
            write.xlsx2(
                as.data.frame(downscale.kk),
                fs::path(
                    save.dir,
                    paste(idx, cluster, "KEGG.xlsx", sep = ".")
                    ),
                sheetName = "avg_log2FC < 0",
                append = TRUE
            )
        } else {
            p2 <- NULL
        }
        ggsave(
            fs::path(save.dir, paste(idx, cluster, "KEGG.pdf", sep = ".")),
            p1 + p2,
            width = 14
        )
    }
}
sapply(idx.full, enrichment1)

# %% DE 2 tumor(tumor + blood + cell densly) vs para (para + junc)
differentialExpression2 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    compare.ident.1 <- list(
        "21B-603-5" = c(
            "Blood vessel rich tumor area",
            "Tumor cell densely populated area",
            "Tumor area 1"
            ),
        "22F-10823-3" = c(
            "Blood vessel rich tumor area",
            "Tumor cell densely populated area",
            "Tumor area 2"
        )
    )
    compare.ident.2 <- list(
        "21B-603-5" = c("Junction area", "Parancerous area"),
        "22F-10823-3" = c("Junction area", "Parancerous area", "Junction area 2")
    )
    markers <- FindMarkers(
        seurat.obj,
        logfc.threshold = 0.1,
        ident.1 = compare.ident.1[[idx]],
        ident.2 = compare.ident.2[[idx]],
        verbose = FALSE
    )
    save.dir <- fs::path(save.dirs[[idx]], "肿瘤vs瘤旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, paste0(idx, ".差异表达.csv")))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    bottoms <- markers %>% top_n(n = -50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(seurat.obj, features = c(tops, bottoms)) + NoLegend()
    save.path <- fs::path(save.dir, paste0(idx, ".差异表达热图.pdf"))
    ggsave(save.path, p, height = 15, width = 15)

    return(markers)
}
markers.list2 <- mclapply(
    seurat.list[c("21B-603-5", "22F-10823-3")],
    differentialExpression2
)

# %% GO 2
enrichment2 <- function(idx) {
    upscale.table <- markers.list2[[idx]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- markers.list2[[idx]] %>%
        filter(p_val_adj < 0.05, avg_log2FC < 0)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
    if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
        print(paste(idx))
        next
    }

    save.dir <- fs::path(save.dirs[[idx]], "肿瘤vs瘤旁", "GO")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    upscale.ego <- enrichGO(
        upscale.genes,
        OrgDb = OrgDb,
        keyType = "SYMBOL",
        ont = "ALL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.1
    )
    downscale.ego <- enrichGO(
        downscale.genes,
        OrgDb = OrgDb,
        keyType = "SYMBOL",
        ont = "ALL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.1
    )
    p1 <- dotplot(upscale.ego, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY~., scale = "free") +
        ggtitle("GO Enrichment for log2FC > 0 genes")
    p2 <- dotplot(downscale.ego, split = "ONTOLOGY") +
        facet_grid(ONTOLOGY~., scale = "free") +
        ggtitle("GO Enrichment for log2FC < 0 genes")
    ggsave(
        fs::path(save.dir, paste(idx, "GO.pdf", sep = ".")),
        p1 + p2,
        width = 14,
        height = 15
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir, paste(idx, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir, paste(idx, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )

    save.dir <- fs::path(save.dirs[[idx]], "肿瘤vs瘤旁", "KEGG")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    upscale.ids <- mapIds(
        OrgDb, keys = upscale.genes, column = "ENTREZID", keytype = "SYMBOL")
    downscale.ids <- mapIds(
        OrgDb, keys = downscale.genes, column = "ENTREZID", keytype = "SYMBOL")
    upscale.kk <- enrichKEGG(
        upscale.ids,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.1
    )
    downscale.kk <- enrichKEGG(
        downscale.ids,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.1
    )
    if (dim(upscale.kk)[1] > 0) {
        p1 <- dotplot(upscale.kk) +
            ggtitle("KEGG Enrichment for log2FC > 0 genes")
        write.xlsx2(
            as.data.frame(upscale.kk),
            fs::path(save.dir, paste(idx, "KEGG.xlsx", sep = ".")),
            sheetName = "avg_log2FC > 0"
        )
    } else {
        p1 <- NULL
    }
    if (dim(downscale.kk)[1] > 0) {
        p2 <- dotplot(downscale.kk) +
            ggtitle("KEGG Enrichment for log2FC < 0 genes")
        write.xlsx2(
            as.data.frame(downscale.kk),
            fs::path(save.dir, paste(idx, "KEGG.xlsx", sep = ".")),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )
    } else {
        p2 <- NULL
    }
    ggsave(
        fs::path(save.dir, paste(idx, "KEGG.pdf", sep = ".")),
        p1 + p2,
        width = 14
    )
}
sapply(c("21B-603-5", "22F-10823-3"), enrichment2)

# %% DE 3 inter tumor
differentialExpression3 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    names(seurat.obj@images) <- NULL
    compare.idents <- list()
    compare.idents[["21B-603-5"]] <- c(
        "Blood vessel rich tumor area",
        "Tumor cell densely populated area",
        "Tumor area 1"
    )
    compare.idents[["22F-10823-3"]] <- c(
        "Blood vessel rich tumor area",
        "Tumor cell densely populated area",
        "Tumor area 2"
    )
    seurat.obj <- subset(seurat.obj, idents = compare.idents[[idx]])

    save.dir <- fs::path(save.dirs[[idx]], "肿瘤各分区差异表达")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)

    markers <- FindAllMarkers(
        seurat.obj, logfc.threshold = 0.1, verbose = FALSE
    )
    write.csv(markers, fs::path(save.dir, paste0(idx, ".差异表达.csv")))

    tops <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
    p <- DoHeatmap(seurat.obj, features = tops$gene) + NoLegend()
    save.path <- fs::path(
        save.dir, paste0(idx, ".差异表达热图.pdf")
    )
    ggsave(save.path, p, height = 20, width = 15)

    return(markers)
}
markers.list3 <- mclapply(
    seurat.list[c("21B-603-5", "22F-10823-3")],
    differentialExpression3
)

# %% GO 3
enrichment3 <- function(idx) {
    clusters <- unique(markers.list3[[idx]]$cluster)
    for (cluster in clusters) {
        upscale.table <- markers.list3[[idx]] %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC > 0)
        downscale.table <- markers.list3[[idx]] %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC < 0)
        upscale.genes <- upscale.table[upscale.table$cluster == cluster, ]$gene
        upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
        downscale.genes <- downscale.table[downscale.table$cluster == cluster, ]$gene
        downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
        if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
            print(paste(idx, cluster))
            next
        }

        save.dir <- fs::path(save.dirs[[idx]], "肿瘤各分区差异表达", "GO")
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
        upscale.ego <- enrichGO(
            upscale.genes,
            OrgDb = OrgDb,
            keyType = "SYMBOL",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        downscale.ego <- enrichGO(
            downscale.genes,
            OrgDb = OrgDb,
            keyType = "SYMBOL",
            ont = "ALL",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        p1 <- dotplot(upscale.ego, split = "ONTOLOGY") +
            facet_grid(ONTOLOGY ~., scale = "free") +
            ggtitle("GO Enrichment for log2FC > 0 genes")
        p2 <- dotplot(downscale.ego, split = "ONTOLOGY") +
            facet_grid(ONTOLOGY~., scale = "free") +
            ggtitle("GO Enrichment for log2FC < 0 genes")
        cluster <- stringr::str_replace_all(cluster, " ", "_")
        ggsave(
            fs::path(save.dir, paste(idx, cluster, "GO.pdf", sep = ".")),
            p1 + p2,
            width = 14,
            height = 15
        )
        write.xlsx2(
            as.data.frame(upscale.ego),
            fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC > 0"
        )
        write.xlsx2(
            as.data.frame(downscale.ego),
            fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )

        save.dir <- fs::path(save.dirs[[idx]], "肿瘤各分区差异表达", "KEGG")
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
        upscale.ids <- mapIds(
            OrgDb, keys = upscale.genes, column = "ENTREZID", keytype = "SYMBOL")
        downscale.ids <- mapIds(
            OrgDb, keys = downscale.genes, column = "ENTREZID", keytype = "SYMBOL")
        upscale.kk <- enrichKEGG(
            upscale.ids,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        downscale.kk <- enrichKEGG(
            downscale.ids,
            organism = "hsa",
            keyType = "ncbi-geneid",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            qvalueCutoff = 0.1
        )
        if (dim(upscale.kk)[1] > 0) {
            p1 <- dotplot(upscale.kk) +
                ggtitle("KEGG Enrichment for log2FC > 0 genes")
            write.xlsx2(
                as.data.frame(upscale.kk),
                fs::path(
                    save.dir,
                    paste(idx, cluster, "KEGG.xlsx", sep = ".")
                    ),
                sheetName = "avg_log2FC > 0"
            )
        } else {
            p1 <- NULL
        }
        if (dim(downscale.kk)[1] > 0) {
            p2 <- dotplot(downscale.kk) +
                ggtitle("KEGG Enrichment for log2FC < 0 genes")
            write.xlsx2(
                as.data.frame(downscale.kk),
                fs::path(
                    save.dir,
                    paste(idx, cluster, "KEGG.xlsx", sep = ".")
                    ),
                sheetName = "avg_log2FC < 0",
                append = TRUE
            )
        } else {
            p2 <- NULL
        }
        ggsave(
            fs::path(save.dir, paste(idx, cluster, "KEGG.pdf", sep = ".")),
            p1 + p2,
            width = 14
        )
    }
}
sapply(c("21B-603-5", "22F-10823-3"), enrichment3)

# %% DE 4 blood vessel vs. others
differentialExpression4 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    markers.list <- list()
    for (cluster in unique(Idents(seurat.obj))) {
        if (cluster == "Blood vessel rich tumor area") next
        subset.obj <- seurat.obj
        names(subset.obj@images) <- NULL
        subset.obj <- subset(
            subset.obj,
            idents = c("Blood vessel rich tumor area", cluster)
        )
        markers <- FindAllMarkers(
            subset.obj,
            logfc.threshold = 0.1,
            verbose = FALSE
        )
        save.dir <- fs::path(save.dirs[[idx]], "富血管区对比其他区")
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)

        save.path <- fs::path(
            save.dir, paste(idx, cluster, "区域差异表达.csv", sep = ".")
        )
        write.csv(markers, save.path)

        tops <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
        p <- DoHeatmap(subset.obj, features = tops$gene) + NoLegend()
        save.path <- fs::path(
            save.dir, paste(idx, cluster, "区域差异表达热图.pdf", sep = ".")
        )
        ggsave(save.path, p, height = 20, width = 15)
        markers.list[[cluster]] <- markers
    }

    return(markers.list)
}

markers.list4 <- mclapply(seurat.list, differentialExpression4)

# %%
mclapply(
    seurat.list, drawGenelist, gene.list = loadKinase(), save.name = "激酶表达"
)
mclapply(
    seurat.list, drawGenelist, gene.list = loadUbiquitin(), save.name = "泛素表达"
)
mclapply(
    seurat.list, drawGenelist, gene.list = loadTF(), save.name = "转录因子表达"
)
mclapply(
    seurat.list, drawGenelist, gene.list = loadRBP(), save.name = "RBP表达"
)
