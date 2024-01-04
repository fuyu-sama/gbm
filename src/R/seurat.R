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
library('ggsignif')
library('Seurat')
library('clusterProfiler')
library('org.Hs.eg.db')

# options
options(mc.cores = 8)
options(future.globals.maxSize = 500 * 1024 ^ 3)

# public variables
ggsave <- function(...) suppressMessages(ggplot2::ggsave(...))
WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
idx.full <- c("21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2")
idx.list <- as.list(idx.full)
names(idx.list) <- idx.full

save.dirs <- lapply(
    idx.list,
    function(idx) {
        save.dir <- fs::path(WORKDIR, "results", idx)
        if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
        return(save.dir)
    }
)

OrgDb <- org.Hs.eg.db

# public functions
saveSeuratList <- function() {
    saveRDS(seurat.list, fs::path(WORKDIR, "results", "seurat.list.rds"))
}

loadSeuratList <- function() {
    seurat.list <<- readRDS(fs::path(WORKDIR, "results", "seurat.list.rds"))
}

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
        "22F-21576-1" = 0.14,
        "22F-23738-2" = 0.11
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

# %% annotation
regionAnnotation <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    region.annotations <- list()
    region.annotations[["21B-603-5"]] <- c(
        "4" = "Blood vessel rich area",
        "2" = "Tumor cell densely populated area",
        "0" = "Tumor area 1",
        "1" = "Junction area",
        "3" = "Parancerous area"
    )
    region.annotations[["22F-10823-3"]] <- c(
        "2" = "Blood vessel rich area",
        "0" = "Tumor cell densely populated area",
        "5" = "Tumor area 2",
        "3" = "Junction area 2",
        "1" = "Junction area",
        "4" = "Parancerous area"
    )
    region.annotations[["22F-21576-1"]] <- c(
        "0" = "Blood vessel rich area",
        "2" = "Tumor cell densely populated area",
        "3" = "Tumor area 3",
        "1" = "Tumor area 4"
    )
    region.annotations[["22F-23738-2"]] <- c(
        "2" = "Blood vessel rich area",
        "0" = "Tumor cell densely populated area",
        "1" = "Tumor area 5",
        "4" = "Tumor area 6",
        "3" = "Tumor area 7",
        "5" = "Tumor area 8"
    )
    seurat.obj <- RenameIdents(seurat.obj, region.annotations[[idx]])
    p1 <- DimPlot(seurat.obj, reduction = "tsne", label = TRUE)
    p2 <- SpatialDimPlot(seurat.obj, label.size = 3, alpha = 0.6, label = TRUE)
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区.pdf"))
    ggsave(save.path, p1 + p2, width = 21, height = 7)
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, regionAnnotation)

# %% cell cycle
cellCycle <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    seurat.obj <- CellCycleScoring(
        seurat.obj,
        g2m.features = cc.genes$g2m.genes,
        s.features = cc.genes$s.genes
    )
    ggsave(
        fs::path(save.dirs[[idx]], paste0(idx, ".CellCycle.Score.pdf")),
        SpatialFeaturePlot(seurat.obj, c("S.Score", "G2M.Score")),
        width = 14
    )
    ggsave(
        fs::path(save.dirs[[idx]], paste0(idx, ".CellCycle.Phase.pdf")),
        SpatialDimPlot(seurat.obj, "Phase")
    )
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, cellCycle)

# %% cell cycle plot
drawCellCycleViolin <- function(seurat.obj) {
    y_position <- list(
        "21B-603-5" = c(6.5, 7.5, 7),
        "22F-10823-3" = c(8, 9, 8.5),
        "22F-21576-1" = c(6.5, 7.5, 7),
        "22F-23738-2" = c(8, 9, 8.5)
    )
    resolutions <- list(
        "21B-603-5" = "SCT_snn_res.0.1",
        "22F-10823-3" = "SCT_snn_res.0.1",
        "22F-21576-1" = "SCT_snn_res.0.08",
        "22F-23738-2" = "SCT_snn_res.0.06"
    )
    idx <- names(seurat.obj@images)
    pdim <- SpatialDimPlot(seurat.obj, group.by = resolutions[[idx]])
    names(seurat.obj@images) <- NULL
    DefaultAssay(seurat.obj) <- "Spatial"
    regions <- seurat.obj@meta.data[resolutions[[idx]]]
    g1.spots <- colnames(subset(seurat.obj, subset = Phase == "G1"))
    g2m.spots <- colnames(subset(seurat.obj, subset = Phase == "G2M"))
    s.spots <- colnames(subset(seurat.obj, subset = Phase == "S"))
    draw.df <- rbind(
        data.frame(
            "CD58" = GetAssayData(seurat.obj)["CD58", g1.spots],
            "Phase" = rep("G1", length(g1.spots))),
        data.frame(
            "CD58" = GetAssayData(seurat.obj)["CD58", g2m.spots],
            "Phase" = rep("G2M", length(g2m.spots))),
        data.frame(
            "CD58" = GetAssayData(seurat.obj)["CD58", s.spots],
            "Phase" = rep("S", length(s.spots)))
    )
    draw.df$Region <- regions[rownames(draw.df), resolutions[[idx]]]
    draw.df <- draw.df[draw.df$CD58 > 0, ]
    comparisons <- list(c("G1", "G2M"), c("G1", "S"), c("G2M", "S"))
    p <- ggplot(draw.df, aes(x = Phase, y = CD58)) +
        geom_violin() +
        facet_wrap(~Region) +
        geom_signif(comparisons = comparisons, y_position = y_position[[idx]], map_signif_level = TRUE)

    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".CD58.Violin.pdf"))
    ggsave(save.path, p + pdim, width = 17)

    p <- ggplot(draw.df, aes(x = Region, fill = Phase)) + geom_bar(position="fill")
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".Phase.stack.pdf"))
    ggsave(save.path, p + pdim, width = 17)
}
mclapply(seurat.list, drawCellCycleViolin)

# %% DE 1
differentialExpression1 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    markers <- FindAllMarkers(
        seurat.obj, logfc.threshold = 0.1, verbose = FALSE
    )
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区域差异表达.csv"))
    write.csv(markers, save.path)

    tops <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
    p <- DoHeatmap(seurat.obj, features = tops$gene) + NoLegend()
    save.path <- fs::path(
        save.dirs[[idx]], paste0(idx, ".分区域差异表达热图.pdf")
    )
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
            "Blood vessel rich area",
            "Tumor cell densely populated area",
            "Tumor area 1"
            ),
        "22F-10823-3" = c(
            "Blood vessel rich area",
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
sapply(c("21B-603-5", "22F-10823-3") , enrichment2)

# %% DE 3 inter tumor
differentialExpression3 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    names(seurat.obj@images) <- NULL
    compare.idents <- list()
    compare.idents[["21B-603-5"]] <- c(
        "Blood vessel rich area",
        "Tumor cell densely populated area",
        "Tumor area 1"
    )
    compare.idents[["22F-10823-3"]] <- c(
        "Blood vessel rich area",
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
        if (cluster == "Blood vessel rich area") next
        subset.obj <- seurat.obj
        names(subset.obj@images) <- NULL
        subset.obj <- subset(
            subset.obj,
            idents = c("Blood vessel rich area", cluster)
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

# %% read DE 1
differentialExpression1 <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    read.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区域差异表达.csv"))
    markers <- read.csv(read.path, row.names = 1)
    markers$cluster <- factor(
        markers$cluster,
        levels = levels(Idents(seurat.obj))
    )
    return(markers)
}
markers.list1 <- mclapply(seurat.list, differentialExpression1)

# %% kinase
kinase.df <- jsonlite::fromJSON(fs::path(WORKDIR, "Data", "klifs.net.json"))
drawKinase <- function(seurat.obj) {
    idx <- names(seurat.obj@images)

    marker.genes <- markers.list1[[idx]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% kinase.df$name)

    write.csv(
        marker.genes,
        fs::path(save.dirs[[idx]], paste(idx, "激酶表达.csv", sep = "."))
    )


    p <- DoHeatmap(seurat.obj, features = marker.genes$gene)
    save.path <- fs::path(
        save.dirs[[idx]], paste(idx, "激酶表达.pdf", sep = ".")
    )
    ggsave(save.path, p, height = 20, width = 15)
}
mclapply(seurat.list, drawKinase)

# %%
kinase.df <- jsonlite::fromJSON(fs::path(WORKDIR, "Data", "klifs.net.json"))
drawKinaseTogether <- function(idx, level) {
    marker.genes.1 <- markers.list1[[idx[1]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% kinase.df$name)
    marker.genes.2 <- markers.list1[[idx[2]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% kinase.df$name)
    marker.genes <- rbind(marker.genes.1, marker.genes.2) %>%
        group_by(cluster)
    draw.genes <- unique(marker.genes$gene)
    #draw.genes <- intersect(marker.genes.1$gene, marker.genes.2$gene)

    p1 <- DoHeatmap(seurat.list[[idx[1]]], features = draw.genes)
    p2 <- DoHeatmap(seurat.list[[idx[2]]], features = draw.genes)
    p <- p1 + p2
    save.path <- fs::path(WORKDIR, "results", paste0(level, "激酶表达.1.pdf"))
    ggsave(save.path, p, height = 20, width = 30)

    draw.seurat <- merge(seurat.list[[idx[1]]], seurat.list[[idx[2]]])
    p <- DoHeatmap(draw.seurat, features = draw.genes)
    save.path <- fs::path(WORKDIR, "results", paste0(level, "激酶表达.2.pdf"))
    ggsave(save.path, p, height = 20, width = 15)
}
drawKinaseTogether(c(idx.full[1], idx.full[2]), "IV")
drawKinaseTogether(c(idx.full[3], idx.full[4]), "GBM")

# %% ubiquitin
flag <- TRUE
for (ubi in c("E1", "E2", "E3")) {
    df.path <- fs::path(WORKDIR, "Data", "iuucd", paste0("ALL_", ubi, ".txt"))
    read.df <- read.table(df.path, sep = "\t")
    read.df$V6 <- sapply(strsplit(read.df$V4, ";  "), function(x) return(x[1]))
    if (flag) {
        ubiquitin.df <- read.df
        flag <- FALSE
    } else {
        ubiquitin.df <- rbind(ubiquitin.df, read.df)
    }
}
drawUbiquitin <- function(seurat.obj) {
    idx <- names(seurat.obj@images)

    marker.genes <- markers.list1[[idx]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% ubiquitin.df$V6)

    write.csv(
        marker.genes,
        fs::path(save.dirs[[idx]], paste(idx, "泛素表达.csv", sep = "."))
    )

    p <- DoHeatmap(seurat.obj, features = marker.genes$gene)
    save.path <- fs::path(
        save.dirs[[idx]], paste(idx, "泛素表达.pdf", sep = ".")
    )
    ggsave(save.path, p, height = 20, width = 15)
}
mclapply(seurat.list, drawUbiquitin)

# %%
drawUbiquitinTogether <- function(idx, level) {
    marker.genes.1 <- markers.list1[[idx[1]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% ubiquitin.df$V6)
    marker.genes.2 <- markers.list1[[idx[2]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% ubiquitin.df$V6)

    marker.genes <- rbind(marker.genes.1, marker.genes.2) %>%
        group_by(cluster)
    draw.genes <- unique(marker.genes$gene)
    #draw.genes <- intersect(marker.genes.1$gene, marker.genes.2$gene)

    p1 <- DoHeatmap(seurat.list[[idx[1]]], features = draw.genes)
    p2 <- DoHeatmap(seurat.list[[idx[2]]], features = draw.genes)
    p <- p1 + p2
    save.path <- fs::path(WORKDIR, "results", paste0(level, "泛素表达.1.pdf"))
    ggsave(save.path, p, height = 20, width = 30)

    draw.seurat <- merge(seurat.list[[idx[1]]], seurat.list[[idx[2]]])
    p <- DoHeatmap(draw.seurat, features = draw.genes)
    save.path <- fs::path(WORKDIR, "results", paste0(level, "泛素表达.2.pdf"))
    ggsave(save.path, p, height = 20, width = 15)
}
drawUbiquitinTogether(c(idx.full[1], idx.full[2]), "IV")
drawUbiquitinTogether(c(idx.full[3], idx.full[4]), "GBM")

# %% TF
tf.df <- read.csv(fs::path(WORKDIR, "Data", "DatabaseExtract_v_1.01.csv"))
tf.df <- filter(tf.df, TF.assessment == "Known motif")
drawTF <- function(seurat.obj) {
    idx <- names(seurat.obj@images)

    marker.genes <- markers.list1[[idx]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% tf.df$HGNC.symbol)

    write.csv(
        marker.genes,
        fs::path(save.dirs[[idx]], paste(idx, "转录因子表达.csv", sep = "."))
    )

    p <- DoHeatmap(seurat.obj, features = marker.genes$gene)
    save.path <- fs::path(
        save.dirs[[idx]], paste(idx, "转录因子表达.pdf", sep = ".")
    )
    ggsave(save.path, p, height = 20, width = 15)
}
mclapply(seurat.list, drawTF)

# %%
drawTFTogether <- function(idx, level) {
    marker.genes.1 <- markers.list1[[idx[1]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% tf.df$HGNC.symbol)
    marker.genes.2 <- markers.list1[[idx[2]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% tf.df$HGNC.symbol)

    marker.genes <- rbind(marker.genes.1, marker.genes.2) %>%
        group_by(cluster)
    draw.genes <- unique(marker.genes$gene)
    #draw.genes <- intersect(marker.genes.1$gene, marker.genes.2$gene)

    p1 <- DoHeatmap(seurat.list[[idx[1]]], features = draw.genes)
    p2 <- DoHeatmap(seurat.list[[idx[2]]], features = draw.genes)
    p <- p1 + p2
    save.path <- fs::path(
        WORKDIR, "results", paste0(level, "转录因子表达.1.pdf"))
    ggsave(save.path, p, height = 20, width = 30)

    draw.seurat <- merge(seurat.list[[idx[1]]], seurat.list[[idx[2]]])
    p <- DoHeatmap(draw.seurat, features = draw.genes)
    save.path <- fs::path(
        WORKDIR, "results", paste0(level, "转录因子表达.2.pdf"))
    ggsave(save.path, p, height = 20, width = 15)
}
drawTFTogether(c(idx.full[1], idx.full[2]), "IV")
drawTFTogether(c(idx.full[3], idx.full[4]), "GBM")

# %% RBP
rbp.df <- read.csv(fs::path(WORKDIR, "Data", "cancers-751631-suppl-final.csv"))
drawRBP <- function(seurat.obj) {
    idx <- names(seurat.obj@images)

    marker.genes <- markers.list1[[idx]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% rbp.df$ID)

    write.csv(
        marker.genes,
        fs::path(save.dirs[[idx]], paste(idx, "RBP表达.csv", sep = "."))
    )

    p <- DoHeatmap(seurat.obj, features = marker.genes$gene)
    save.path <- fs::path(
        save.dirs[[idx]], paste(idx, "RBP表达.pdf", sep = ".")
    )
    ggsave(save.path, p, height = 30, width = 15)
}
mclapply(seurat.list, drawRBP)

# %%
drawRBPTogether <- function(idx, level) {
    marker.genes.1 <- markers.list1[[idx[1]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% rbp.df$ID)
    marker.genes.2 <- markers.list1[[idx[2]]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% rbp.df$ID)

    marker.genes <- rbind(marker.genes.1, marker.genes.2) %>%
        group_by(cluster)
    draw.genes <- unique(marker.genes$gene)
    #draw.genes <- intersect(marker.genes.1$gene, marker.genes.2$gene)

    p1 <- DoHeatmap(seurat.list[[idx[1]]], features = draw.genes)
    p2 <- DoHeatmap(seurat.list[[idx[2]]], features = draw.genes)
    p <- p1 + p2
    save.path <- fs::path(WORKDIR, "results", paste0(level, "RBP表达.1.pdf"))
    ggsave(save.path, p, height = 20, width = 30)

    draw.seurat <- merge(seurat.list[[idx[1]]], seurat.list[[idx[2]]])
    p <- DoHeatmap(draw.seurat, features = draw.genes)
    save.path <- fs::path(WORKDIR, "results", paste0(level, "RBP表达.2.pdf"))
    ggsave(save.path, p, height = 20, width = 15)
}
drawUbiquitinTogether(c(idx.full[1], idx.full[2]), "IV")
drawUbiquitinTogether(c(idx.full[3], idx.full[4]), "GBM")
