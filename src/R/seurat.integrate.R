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
save.dirs[["int"]] <- fs::path(WORKDIR, "results", "integrate")
if (! fs::dir_exists(save.dirs[["int"]])) fs::dir_create(save.dirs[["int"]])

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

# %% load Seurat list
loadSeuratList()

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
        nfeatures = 9000,
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
    upscale.table <- integrate.markers1[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers1[[region]] %>%
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

# %% integrate de 4
differentialExpressionInt4 <- function(level, findmarker = TRUE) {
    subset.obj <- integrate.obj
    compare.ident <- c(
        "Tumor cell densely populated area", "Blood vessel rich area")
    names(subset.obj@images) <- NULL

    para.spots <- colnames(subset(subset.obj, idents = "Parancerous area"))
    compare.spots <- colnames(subset(
            subset.obj[, subset.obj@meta.data$level == level],
            idents = compare.ident
            ))

    subset.obj <- subset(subset.obj, cells = c(compare.spots, para.spots))
    group <- ifelse(
        colnames(subset.obj) %in% para.spots,
        "Parancerous area",
        paste0(level, ".Tumor")
    )
    subset.obj@meta.data$group <- group

    save.dir <- fs::path(save.dirs[["int"]], "tumor-vs-parancerous")
    file.name <- paste0(level, ".Tumor.vs.Parancerous_area")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    if (findmarker) {
        markers <- FindMarkers(
            subset.obj,
            paste0(level, ".Tumor"),
            logfc.threshold = 1.5,
            group.by = "group",
            verbose = FALSE
        )
        save.path <- fs::path(save.dir, paste0(file.name, ".csv"))
        write.csv(markers, save.path)
    } else {
        markers <- integrate.markers4[[level]]
    }

    tops <- c(
        rownames(top_n(markers, n = 20, wt = avg_log2FC)),
        rownames(top_n(markers, n = -20, wt = avg_log2FC))
    )
    p <- DoHeatmap(subset.obj, features = tops, group.by = "group") + NoLegend()
    save.path <- fs::path(save.dir, paste0(file.name, ".pdf"))
    ggsave(save.path, p)

    tf.path <- fs::path(WORKDIR, "Data", "DatabaseExtract_v_1.01.csv")
    tf.df <- read.csv(tf.path, row.names = 1)
    tf.df <- filter(tf.df, Is.TF. == "Yes", TF.assessment == "Known motif")
    tf.genes <- tf.df$HGNC.symbol
    draw.df <- markers %>%
        filter(p_val_adj < 0.05, rownames(markers) %in% tf.genes) %>%
        arrange(desc(avg_log2FC))
    write.csv(draw.df, fs::path(save.dir, paste0(file.name, ".TFs.csv")))
    if (dim(draw.df)[1] > 100) {
        draw.df <- rbind(
            top_n(draw.df, 50, wt = avg_log2FC),
            top_n(draw.df, -50, wt = avg_log2FC)
        )
    }
    p <- DoHeatmap(subset.obj, features = rownames(draw.df), group.by = "group")
    ggsave(fs::path(save.dir, paste0(level, "-TFs-heatmap.pdf")), p, height = 13)

    rbp.path <- fs::path(WORKDIR, "Data", "cancers-751631-suppl-final.csv")
    rbp.df <- read.csv(rbp.path, row.names = 1)
    rbp.genes <- rownames(rbp.df)
    draw.df <- markers %>%
        filter(p_val_adj < 0.05, rownames(markers) %in% rbp.genes) %>%
        arrange(desc(avg_log2FC))
    write.csv(draw.df, fs::path(save.dir, paste0(file.name, ".RBPs.csv")))
    if (dim(draw.df)[1] > 100) {
        draw.df <- rbind(
            top_n(draw.df, 50, wt = avg_log2FC),
            top_n(draw.df, -50, wt = avg_log2FC)
        )
    }
    p <- DoHeatmap(subset.obj, features = rownames(draw.df), group.by = "group")
    ggsave(fs::path(save.dir, paste0(level, "-RBPs-heatmap.pdf")), p, height = 13)

    return(markers)
}

integrate.markers4 <- mclapply(
    c("IV", "GBM"),
    differentialExpressionInt4,
    findmarker = TRUE
)
names(integrate.markers4) <- c("IV", "GBM")

# %% integrate de 4 draw
differentialExpressionInt4Multi <- function(level) {
    source(fs::path(WORKDIR, "src", "R", "utils.R"))
    subset.obj <- integrate.obj
    compare.ident <- c(
        "Tumor cell densely populated area", "Blood vessel rich area")
    names(subset.obj@images) <- NULL

    para.spots <- colnames(subset(subset.obj, idents = "Parancerous area"))
    compare.spots <- colnames(subset(
            subset.obj[, subset.obj@meta.data$level == level],
            idents = compare.ident
            ))

    subset.obj <- subset(subset.obj, cells = c(compare.spots, para.spots))
    group <- ifelse(
        colnames(subset.obj) %in% para.spots,
        "Parancerous area",
        paste0(level, ".Tumor")
    )
    subset.obj@meta.data$group <- group

    save.dir <- fs::path(save.dirs[["int"]], "tumor-vs-parancerous")
    file.name <- paste0(level, ".Tumor.vs.Parancerous_area")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    markers <- integrate.markers4[[level]]

    tf.path <- fs::path(WORKDIR, "Data", "DatabaseExtract_v_1.01.csv")
    tf.df <- read.csv(tf.path, row.names = 1)
    tf.df <- filter(tf.df, Is.TF. == "Yes", TF.assessment == "Known motif")
    tf.genes <- tf.df$HGNC.symbol
    draw.df <- markers %>%
        filter(p_val_adj < 0.05, rownames(markers) %in% tf.genes) %>%
        arrange(desc(avg_log2FC))
    if (dim(draw.df)[1] > 100) {
        draw.df <- rbind(
            top_n(draw.df, 50, wt = avg_log2FC),
            top_n(draw.df, -50, wt = avg_log2FC)
        )
    }
    p <- DoMultiBarHeatmap(
        subset.obj,
        features = rownames(draw.df),
        group.by = "group",
        additional.group.by = c("sample")
    )
    ggsave(fs::path(save.dir, paste0(level, "-TFs-heatmap.1.pdf")), p, height = 13)

    rbp.path <- fs::path(WORKDIR, "Data", "cancers-751631-suppl-final.csv")
    rbp.df <- read.csv(rbp.path, row.names = 1)
    rbp.genes <- rownames(rbp.df)
    draw.df <- markers %>%
        filter(p_val_adj < 0.05, rownames(markers) %in% rbp.genes) %>%
        arrange(desc(avg_log2FC))
    if (dim(draw.df)[1] > 100) {
        draw.df <- rbind(
            top_n(draw.df, 50, wt = avg_log2FC),
            top_n(draw.df, -50, wt = avg_log2FC)
        )
    }
    p <- DoMultiBarHeatmap(
        subset.obj,
        features = rownames(draw.df),
        group.by = "group",
        additional.group.by = c("sample")
    )
    ggsave(fs::path(save.dir, paste0(level, "-RBPs-heatmap.1.pdf")), p, height = 13)
}

mclapply(c("IV", "GBM"), differentialExpressionInt4Multi)
