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
save.dirs[["int"]] <- fs::path(WORKDIR, "results", "integrate")
OrgDb <- org.Hs.eg.db

# public functions
saveSeuratList <- function() {
    saveRDS(seurat.list, fs::path(WORKDIR, "results", "seurat.list.rds"))
}

loadSeuratList <- function() {
    seurat.list <<- readRDS(fs::path(WORKDIR, "results", "seurat.list.rds"))
}

saveIntegrate <- function() {
    saveRDS(integrate.obj, fs::path(WORKDIR, "results", "integrate.rds"))
}

loadIntegrate <- function() {
    integrate.obj <<- readRDS(fs::path(WORKDIR, "results", "integrate.rds"))
}

saveMerge <- function() {
    saveRDS(merge.obj, fs::path(WORKDIR, "results", "merge.rds"))
}

loadMerge <- function() {
    merge.obj <<- readRDS(fs::path(WORKDIR, "results", "merge.rds"))
}

saveMarkerList <- function() {
    saveRDS(markers.list, fs::path(WORKDIR, "results", "markers.list.rds"))
}

loadMarkerList <- function() {
    markers.list <<- readRDS(fs::path(WORKDIR, "results", "markers.list.rds"))
}

# %% read data
readData <- function(idx) {
    read.dir <- fs::path(WORKDIR, "spaceranger", idx, "outs")
    filename <- "filtered_feature_bc_matrix.h5"
    if (idx == "21B-603-5") filename <- "raw_feature_bc_matrix.h5"
    seurat.obj <- Load10X_Spatial(
        read.dir,
        filename = filename,
        filter.matrix = FALSE,
        image = Read10X_Image(fs::path(read.dir, "spatial"), filter.matrix = FALSE)
    )
    names(seurat.obj@images) <- idx
    return(seurat.obj)
}
seurat.list <- mclapply(idx.list, readData)

# %% merge data
mergeData <- function(seurat.list) {
    for (seurat.obj in seurat.list) {
        idx <- names(seurat.obj@images)
        names(seurat.obj@images) <- "slice1_"
        seurat.obj@meta.data["sample"] <- rep(idx, dim(seurat.obj)[2])
        seurat.obj <- RenameCells(seurat.obj, add.cell.id = idx)
        names(seurat.obj@images) <- idx
        seurat.list[[idx]] <- seurat.obj
    }
    merge.obj <- merge(seurat.list[[idx.full[1]]], seurat.list[[idx.full[2]]])
    merge.obj <- merge(merge.obj, seurat.list[[idx.full[3]]])
    merge.obj <- merge(merge.obj, seurat.list[[idx.full[4]]])
    return(merge.obj)
}
merge.obj <- mergeData(seurat.list)

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
        "21B-603-5" = 0.14,
        "22F-10823-3" = 0.18,
        "22F-21576-1" = 0.14,
        "22F-23738-2" = 0.12
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
        "1" = "Junction area",
        "2" = "Tumor cell densely populated area",
        "3" = "Parancerous area",
        "4" = "Blood vessel rich area",
        "0" = "21B-603-5_0"
    )
    region.annotations[["22F-10823-3"]] <- c(
        "0" = "Tumor cell densely populated area",
        "1" = "Blood vessel rich area",
        "2" = "Junction area",
        "3" = "Parancerous area",
        "4" = "22F-10823-3_4",
        "5" = "22F-10823-3_5"
    )
    region.annotations[["22F-21576-1"]] <- c(
        "0" = "Blood vessel rich area",
        "2" = "Tumor cell densely populated area",
        "1" = "22F-21576-1_1",
        "3" = "22F-21576-1_3"
    )
    region.annotations[["22F-23738-2"]] <- c(
        "0" = "Tumor cell densely populated area",
        "2" = "Blood vessel rich area",
        "1" = "22F-23738-2_1",
        "3" = "22F-23738-2_3",
        "4" = "22F-23738-2_4",
        "5" = "22F-23738-2_5"
    )
    seurat.obj <- RenameIdents(seurat.obj, region.annotations[[idx]])
    p1 <- DimPlot(seurat.obj, reduction = "tsne", label = TRUE)
    p2 <- SpatialDimPlot(seurat.obj, label.size = 3, alpha = 0.6, label = TRUE)
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".annotation.pdf"))
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

# %% de
differentialExpression <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    markers <- FindAllMarkers(
        seurat.obj, logfc.threshold = 0.1, verbose = FALSE
    )
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".AllMarkers.csv"))
    write.csv(markers, save.path)

    tops <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    p <- DoHeatmap(seurat.obj, features = tops$gene) + NoLegend()
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".AllMarkers.pdf"))
    ggsave(save.path, p)

    return(markers)
}
markers.list <- mclapply(seurat.list, differentialExpression)

# %% GO
goEnrichment <- function(idx) {
    clusters <- unique(markers.list[[idx]]$cluster)
    for (cluster in clusters) {
        upscale.table <- markers.list[[idx]] %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC > 0)
        downscale.table <- markers.list[[idx]] %>%
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
            fs::path(save.dirs[[idx]], paste(idx, cluster, "GO.pdf", sep = ".")),
            p1 + p2,
            width = 14,
            height = 15
        )
        write.xlsx2(
            as.data.frame(upscale.ego),
            fs::path(save.dirs[[idx]], paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC > 0"
        )
        write.xlsx2(
            as.data.frame(downscale.ego),
            fs::path(save.dirs[[idx]], paste(idx, cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )
    }
}
sapply(idx.full, goEnrichment)

# %% integrate data
integrateData <- function(seurat.list) {
    cons <- list(
        "21B-603-5" = "IV",
        "22F-10823-3" = "IV",
        "22F-21576-1" = "GBM",
        "22F-23738-2" = "GBM"
    )
    for (seurat.obj in seurat.list) {
        idx <- names(seurat.obj@images)
        names(seurat.obj@images) <- "slice1_"
        seurat.obj@meta.data["sample"] <- rep(idx, dim(seurat.obj)[2])
        seurat.obj@meta.data["level"] <- rep(cons[[idx]], dim(seurat.obj)[2])
        seurat.obj <- RenameCells(seurat.obj, add.cell.id = idx)
        names(seurat.obj@images) <- idx
        seurat.list[[idx]] <- seurat.obj
    }
    integrate.features <- SelectIntegrationFeatures(
        seurat.list,
        nfeatures = 5000,
        verbose = FALSE
    )
    seurat.list <- PrepSCTIntegration(
        seurat.list,
        anchor.features = integrate.features,
        verbose = FALSE
    )
    seurat.list <- lapply(
        seurat.list,
        FUN = RunPCA,
        features = integrate.features,
        verbose = FALSE
    )
    anchors <- FindIntegrationAnchors(
        seurat.list,
        anchor.features = integrate.features,
        normalization.method = "SCT",
        reduction = "rpca",
        k.anchor = 10,
        dims = 1:20,
        verbose = FALSE
    )
    integrate.obj <- IntegrateData(
        anchors,
        normalization.method = "SCT",
        dims = 1:20,
        verbose = FALSE
    )

    integrate.obj <- RunPCA(
        integrate.obj, assay = "integrated", verbose = FALSE)
    integrate.obj <- FindNeighbors(
        integrate.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    integrate.obj <- RunUMAP(
        integrate.obj, reduction = "pca", dims = 1:30, verbose = FALSE)
    integrate.obj <- RunTSNE(
        integrate.obj, reduction = "pca", dims = 1:30, verbose = FALSE)

    p1 <- DimPlot(integrate.obj, reduction = "tsne", cols = DiscretePalette(13))
    p2 <- DimPlot(integrate.obj, reduction = "tsne", group.by = "sample")
    ggsave(fs::path(save.dirs[["int"]], "integrate.tsne.pdf"), p1 + p2, width = 20)

    return(integrate.obj)
}
integrate.obj <- integrateData(seurat.list)

# %% integrate de 1
differentialExpressionInt1 <- function(compare.ident) {
    subset.obj <- integrate.obj
    names(subset.obj@images) <- NULL
    subset.obj <- subset(subset.obj, idents = compare.ident)
    subset.obj@images <- integrate.obj@images

    markers <- FindMarkers(
        subset.obj,
        "IV",
        "GBM",
        logfc.threshold = 0.1,
        group.by = "level",
        verbose = FALSE
    )
    compare.ident <- stringr::str_replace_all(compare.ident, " ", "_")
    file.name <- paste(compare.ident, "IV.vs.GBM", sep = ".")
    save.dir <- fs::path(save.dirs[["int"]], "iv-vs-gbm")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    save.path <- fs::path(save.dir, paste0(file.name, ".csv"))
    write.csv(markers, save.path)

    tops <- c(
        rownames(top_n(markers, n = 20, wt = avg_log2FC)),
        rownames(top_n(markers, n = -20, wt = avg_log2FC))
    )
    p <- DoHeatmap(subset.obj, group.by = "level", features = tops) +
        NoLegend()
    save.path <- fs::path(save.dir, paste0(file.name, ".pdf"))
    ggsave(save.path, p)

    return(markers)
}

named.regions <- list(
    "Tumor cell densely populated area" = "Tumor cell densely populated area",
    "Blood vessel rich area" = "Blood vessel rich area"
)
integrate.markers1 <- mclapply(named.regions, differentialExpressionInt1)

# %% enrichment for integrate 1
enrichmentInt1 <- function(region) {
    upscale.table <- integrate.markers[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC < 0)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
    if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
        print(region)
        next
    }

    region <- stringr::str_replace_all(region, " ", "_")
    save.dir <- fs::path(save.dirs[["int"]], "iv-vs-gbm")

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
        fs::path(save.dir, paste(region, "GO.pdf", sep = ".")),
        p1 + p2,
        width = 14,
        height = 15
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir, paste(region, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir, paste(region, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )

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
    p1 <- dotplot(upscale.kk) +
        ggtitle("KEGG Enrichment for log2FC > 0 genes")
    p2 <- dotplot(downscale.kk) +
        ggtitle("KEGG Enrichment for log2FC < 0 genes")
    ggsave(
        fs::path(save.dir, paste(region, "KEGG.pdf", sep = ".")),
        p1 + p2,
        width = 14
    )
    write.xlsx2(
        as.data.frame(upscale.kk),
        fs::path(save.dir, paste(region, "KEGG.xlsx", sep = ".")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.kk),
        fs::path(save.dir, paste(region, "KEGG.xlsx", sep = ".")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )
}
sapply(names(integrate.markers1), enrichmentInt1)

# %% integrate de 2
differentialExpressionInt2 <- function(compare.ident) {
    subset.obj <- integrate.obj
    names(subset.obj@images) <- NULL
    new.ident <- paste(subset.obj@meta.data$level, Idents(subset.obj), sep = ".")
    Idents(subset.obj) <- new.ident
    subset.obj <- subset(subset.obj, idents = c(compare.ident, "IV.Parancerous area"))

    markers <- FindMarkers(
        subset.obj,
        compare.ident,
        "IV.Parancerous area",
        logfc.threshold = 0.1,
        verbose = FALSE
    )
    compare.ident <- stringr::str_replace_all(compare.ident, " ", "_")
    file.name <- paste(compare.ident, "vs.Parancerous_area", sep = ".")
    save.dir <- fs::path(save.dirs[["int"]], "other-region-vs-parancerous")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    save.path <- fs::path(save.dir, paste0(file.name, ".csv"))
    write.csv(markers, save.path)

    tops <- c(
        rownames(top_n(markers, n = 20, wt = avg_log2FC)),
        rownames(top_n(markers, n = -20, wt = avg_log2FC))
    )
    p <- DoHeatmap(subset.obj, features = tops) + NoLegend()
    save.path <- fs::path(save.dir, paste0(file.name, ".pdf"))
    ggsave(save.path, p)

    return(markers)
}

named.regions <- list(
    "IV.Tumor cell densely populated area" = "IV.Tumor cell densely populated area",
    "IV.Blood vessel rich area" = "IV.Blood vessel rich area",
    "IV.Junction area" = "IV.Junction area",
    "GBM.Tumor cell densely populated area" = "GBM.Tumor cell densely populated area",
    "GBM.Blood vessel rich area" = "GBM.Blood vessel rich area"
)
integrate.markers2 <- mclapply(named.regions, differentialExpressionInt2)

# %% enrichment for integrate 2
enrichmentInt2 <- function(region) {
    upscale.table <- integrate.markers2[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers2[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC < 0)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
    if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
        print(region)
        next
    }

    region <- stringr::str_replace_all(region, " ", "_")
    save.dir <- fs::path(save.dirs[["int"]], "other-region-vs-parancerous")
    file.name <- paste(region, "vs.Parancerous_area", sep = ".")

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
        fs::path(save.dir, paste(file.name, "GO.pdf", sep = ".")),
        p1 + p2,
        width = 14,
        height = 15
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir, paste(file.name, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir, paste(file.name, "GO.xlsx", sep = ".")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )

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
    p1 <- dotplot(upscale.kk) +
        ggtitle("KEGG Enrichment for log2FC > 0 genes")
    p2 <- dotplot(downscale.kk) +
        ggtitle("KEGG Enrichment for log2FC < 0 genes")
    ggsave(
        fs::path(save.dir, paste(file.name, "KEGG.pdf", sep = ".")),
        p1 + p2,
        width = 14
    )
    write.xlsx2(
        as.data.frame(upscale.kk),
        fs::path(save.dir, paste(file.name, "KEGG.xlsx", sep = ".")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.kk),
        fs::path(save.dir, paste(file.name, "KEGG.xlsx", sep = ".")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )
}
sapply(names(integrate.markers2), enrichmentInt2)

# %% integrate de 3
differentialExpressionInt3 <- function() {
    subset.obj <- integrate.obj
    names(subset.obj@images) <- NULL
    idents <- Idents(subset.obj)
    idents <- idents[! idents %in% c("Parancerous area", "Blood vessel rich area", "Tumor cell densely populated area", "Junction area")]
    subset.obj <- subset(subset.obj, idents = idents)

    markers <- FindMarkers(
        subset.obj,
        "IV",
        "GBM",
        group.by = "level",
        logfc.threshold = 0.1,
        verbose = FALSE
    )
    file.name <- "Tumor_cell_non-dense_area.IV.vs.GBM"
    save.dir <- fs::path(save.dirs[["int"]], "tumor_cell_non-dense_area")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    save.path <- fs::path(save.dir, paste0(file.name, ".csv"))
    write.csv(markers, save.path)

    tops <- c(
        rownames(top_n(markers, n = 20, wt = avg_log2FC)),
        rownames(top_n(markers, n = -20, wt = avg_log2FC))
    )
    p <- DoHeatmap(subset.obj, features = tops, group.by = "level") + NoLegend()
    save.path <- fs::path(save.dir, paste0(file.name, ".pdf"))
    ggsave(save.path, p)

    return(markers)
}
integrate.markers3 <- differentialExpressionInt3()

# %% enrichment for integrate 3
enrichmentInt3 <- function() {
    upscale.table <- integrate.markers3 %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers3 %>%
        filter(p_val_adj < 0.05, avg_log2FC < 0)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]
    if ((length(upscale.genes) <= 0) || (length(downscale.genes) <= 0)) {
        print(region)
        next
    }

    save.dir <- fs::path(save.dirs[["int"]], "tumor_cell_non-dense_area")
    file.name <- "Tumor_cell_non-dense_area.IV.vs.GBM"

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
        fs::path(save.dir, paste0(file.name, ".GO.pdf")),
        p1 + p2,
        width = 14,
        height = 15
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir, paste0(file.name, ".GO.xlsx")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir, paste0(file.name, ".GO.xlsx")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )

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
    p1 <- dotplot(upscale.kk) +
        ggtitle("KEGG Enrichment for log2FC > 0 genes")
    p2 <- dotplot(downscale.kk) +
        ggtitle("KEGG Enrichment for log2FC < 0 genes")
    ggsave(
        fs::path(save.dir, paste0(file.name, ".KEGG.pdf")),
        p1 + p2,
        width = 14
    )
    write.xlsx2(
        as.data.frame(upscale.kk),
        fs::path(save.dir, paste0(file.name, ".KEGG.xlsx")),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.kk),
        fs::path(save.dir, paste0(file.name, ".KEGG.xlsx")),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )
}
enrichmentInt3()

# %% draw tfs
drawTFs <- function() {
    tf.path <- fs::path(WORKDIR, "Data", "DatabaseExtract_v_1.01.csv")
    tf.df <- read.csv(tf.path, row.names = 1)
    tf.df <- filter(tf.df, Is.TF. == "Yes", TF.assessment == "Known motif")
    save.dir <- fs::path(WORKDIR, "gene-expression", "TFs")
    tf.genes <- tf.df$HGNC.symbol

    for (idx in idx.full) {
        draw.df <- markers.list[[idx]] %>%
            group_by(cluster) %>%
            top_n(500, wt = avg_log2FC) %>%
            filter(p_val_adj < 0.05, gene %in% tf.genes)
        draw.seurat <- integrate.obj
        names(draw.seurat@images) <- NULL
        draw.spots <- draw.seurat@meta.data[draw.seurat@meta.data$sample == idx, ]
        draw.seurat <- subset(draw.seurat, cells = rownames(draw.spots))
        p <- DoHeatmap(draw.seurat, features = draw.df$gene)
        ggsave(
            fs::path(WORKDIR, "gene-expression", paste0(idx, "-TFs-heatmap.jpg")),
            p
        )
    }
}
drawTFs()

# %% draw rbps
drawRBPs <- function() {
    rbp.path <- fs::path(WORKDIR, "Data", "cancers-751631-suppl-final.csv")
    rbp.df <- read.csv(rbp.path, row.names = 1)
    rbp.genes <- rownames(rbp.df)

    for (idx in idx.full) {
        draw.df <- markers.list[[idx]] %>% filter(
            p_val_adj < 0.05,
            gene %in% rbp.genes
        )
        draw.seurat <- integrate.obj
        names(draw.seurat@images) <- NULL
        draw.spots <- draw.seurat@meta.data[draw.seurat@meta.data$sample == idx, ]
        draw.seurat <- subset(draw.seurat, cells = rownames(draw.spots))
        p <- DoHeatmap(draw.seurat, features = draw.df$gene)
        ggsave(
            fs::path(WORKDIR, "gene-expression", paste0(idx, "-RBPs-heatmap.jpg")),
            p
        )
    }
}
drawRBPs()
