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
library('Seurat')
library('clusterProfiler')
library('org.Hs.eg.db')

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
source(fs::path(WORKDIR, "src", "R", "utils.R"))

# %% integrate data
loadSeuratList()
integrateData <- function(seurat.list) {
    cons <- list(
        "21B-603-5" = "IV",
        "22F-10823-3" = "IV",
        "22F-21576-1" = "GBM",
        "22F-23738-2" = "GBM"
    )
    all.genes <- rownames(seurat.list[[1]])
    for (seurat.obj in seurat.list) {
        idx <- names(seurat.obj@images)
        names(seurat.obj@images) <- "slice1_"
        seurat.obj@meta.data["sample"] <- rep(idx, dim(seurat.obj)[2])
        seurat.obj@meta.data["level"] <- rep(cons[[idx]], dim(seurat.obj)[2])
        seurat.obj <- RenameCells(seurat.obj, add.cell.id = idx)
        names(seurat.obj@images) <- idx
        seurat.list[[idx]] <- seurat.obj
        all.genes <- intersect(all.genes, rownames(seurat.obj))
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
    anchors <- FindIntegrationAnchors(
        seurat.list,
        anchor.features = integrate.features,
        normalization.method = "SCT",
        k.anchor = 5,
        dims = 1:20,
        verbose = FALSE
    )
    integrate.obj <- IntegrateData(
        anchors,
        normalization.method = "SCT",
        features.to.integrate = all.genes,
        dims = 1:50,
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

    return(integrate.obj)
}
integrate.obj <- integrateData(seurat.list)
#saveIntegrate()

# annotation
regionAnnotation <- function(seurat.obj) {
    region.annotations <- list(
        "21B-603-5.4" = "Blood vessel rich tumor area",
        "22F-10823-3.2" = "Blood vessel rich tumor area",
        "22F-21576-1.0" = "Blood vessel rich tumor area",
        "22F-23738-2.2" = "Blood vessel rich tumor area",
        "21B-603-5.2" = "Tumor cell densely populated area",
        #"22F-10823-3.0" = "Tumor cell densely populated area",
        "22F-10823-3.5" = "Tumor cell densely populated area",
        "22F-21576-1.2" = "Tumor cell densely populated area",
        "22F-23738-2.0" = "Tumor cell densely populated area",
        "21B-603-5.0" = "Tumor area 1",
        #"22F-10823-3.5" = "Tumor area 2",
        "22F-10823-3.0" = "Tumor area 2",
        "22F-21576-1.3" = "Tumor area 3",
        "22F-21576-1.1" = "Tumor area 4",
        "22F-23738-2.1" = "Tumor area 5",
        "22F-23738-2.4" = "Tumor area 6",
        "22F-23738-2.3" = "Tumor area 7",
        "22F-23738-2.5" = "Tumor area 8",
        "21B-603-5.1" = "Junction area",
        "22F-10823-3.1" = "Junction area",
        "22F-10823-3.3" = "Junction area 2",
        "21B-603-5.3" = "Parancerous area",
        "22F-10823-3.4" = "Parancerous area"
    )
    Idents(seurat.obj) <- paste(
        seurat.obj@meta.data$sample, Idents(seurat.obj), sep = ".")
    seurat.obj <- RenameIdents(seurat.obj, region.annotations)

    regions <- unique(as.character(region.annotations))
    cols <- DiscretePalette(length(regions))
    names(cols) <- regions
    cols["Tumor area 2"] <- "#8fb2c9"

    p <- SpatialDimPlot(
        seurat.obj, label.size = 3, alpha = 1, label = TRUE, cols = cols
    )
    save.path <- fs::path(save.dirs[["int"]], "分区.pdf")
    ggsave(save.path, p, width = 28, height = 7)

    p1 <- DimPlot(seurat.obj, reduction = "tsne", cols = cols)
    p2 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "sample")
    ggsave(fs::path(save.dirs[["int"]], "integrate.tsne.pdf"), p1 + p2, width = 20)
    return(seurat.obj)
}
integrate.obj <- regionAnnotation(integrate.obj)

# %%
differentialExpression1 <- function(integrate.obj) {
    markers <- FindAllMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]])
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 20, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(integrate.obj, features = tops)
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list1 <- differentialExpression1(integrate.obj)

# %% 2.IV整合对比
differentialExpression2 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    iv.obj@images <- integrate.obj@images[1:2]
    markers <- FindAllMarkers(iv.obj, logfc.threshold = 0.1, verbose = FALSE)

    save.dir <- fs::path(save.dirs[["int"]], "IV整合对比")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "IV整合差异表格.csv"))

    tops <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
    p <- DoHeatmap(iv.obj, features = tops$gene)
    ggsave(fs::path(save.dir, "IV整合差异热图.pdf"), p)
    return(markers)
}
markers.list2 <- differentialExpression2(integrate.obj)
sapply(
    unique(markers.list2$cluster),
    enrichmentFindAllMarkers,
    markers = markers.list2,
    save.dir = fs::path(save.dirs[["int"]], "IV整合对比")
)

# %% 3.1 IV 富血管 vs. GBM 富血管
differentialExpression3.1 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    iv.blood <- subset(iv.obj, idents = "Blood vessel rich tumor area")
    iv.blood <- colnames(iv.blood)

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "GBM"]
    gbm.blood <- subset(gbm.obj, idents = "Blood vessel rich tumor area")
    gbm.blood <- colnames(gbm.blood)

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = iv.blood,
        ident.2 = gbm.blood,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "IV富血管vsGBM富血管")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(integrate.obj, features = tops, cells = c(iv.blood, gbm.blood))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list3.1 <- differentialExpression3.1(integrate.obj)
enrichmentFindMarkers(
    markers.list3.1, fs::path(save.dirs[["int"]], "IV富血管vsGBM富血管")
)

# %% 3.2 IV 密集 vs. GBM 密集
differentialExpression3.2 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    iv.densly <- subset(iv.obj, idents = "Tumor cell densely populated area")
    iv.densly <- colnames(iv.densly)

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "GBM"]
    gbm.densly <- subset(gbm.obj, idents = "Tumor cell densely populated area")
    gbm.densly <- colnames(gbm.densly)

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = iv.densly,
        ident.2 = gbm.densly,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "IV密集vsGBM密集")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(iv.densly, gbm.densly))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list3.2 <- differentialExpression3.2(integrate.obj)
enrichmentFindMarkers(
    markers.list3.2, fs::path(save.dirs[["int"]], "IV密集vsGBM密集")
)

# %% 4.1 GBM 密集 vs IV 癌旁
differentialExpression4.1 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "GBM"]
    densly.cells <- colnames(
        subset(gbm.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "GBM密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.1 <- differentialExpression4.1(integrate.obj)
enrichmentFindMarkers(
    markers.list4.1, fs::path(save.dirs[["int"]], "GBM密集vsIV癌旁")
)

# %% 4.2 GBM 22F-21576-1 密集 vs IV 癌旁
differentialExpression4.2 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "GBM"]
    gbm.obj <- gbm.obj[, gbm.obj$sample == "22F-21576-1"]
    densly.cells <- colnames(
        subset(gbm.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "GBM22F-21576-1密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.2 <- differentialExpression4.2(integrate.obj)

# %% 4.3 GBM 22F-23738-2 密集 vs IV 癌旁
differentialExpression4.3 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "GBM"]
    gbm.obj <- gbm.obj[, gbm.obj$sample == "22F-23738-2"]
    densly.cells <- colnames(
        subset(gbm.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "GBM22F-23738-2密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.3 <- differentialExpression4.3(integrate.obj)

# %% 4.4 IV 密集 vs IV 癌旁
differentialExpression4.4 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    gbm.obj <- integrate.obj
    names(gbm.obj@images) <- NULL
    gbm.obj <- gbm.obj[, gbm.obj$level == "IV"]
    densly.cells <- colnames(
        subset(gbm.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "GBM密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.4 <- differentialExpression4.4(integrate.obj)
enrichmentFindMarkers(
    markers.list4.4, fs::path(save.dirs[["int"]], "GBM密集vsIV癌旁")
)

# %% 4.5 IV 21B-603-5 密集 vs IV 癌旁
differentialExpression4.5 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    iv.obj <- iv.obj[, iv.obj$sample == "21B-603-5"]
    densly.cells <- colnames(
        subset(iv.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "IV21B-603-5密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.5 <- read.csv(
    fs::path(save.dirs[["int"]], "IV21B-603-5密集vsIV癌旁", "差异表格.csv"),
    row.names = 1
)
markers.list4.5 <- differentialExpression4.5(integrate.obj)

# %% 4.6 IV 22F-10823-3 密集 vs IV 癌旁
differentialExpression4.6 <- function(integrate.obj) {
    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    para.cells <- colnames(subset(iv.obj, idents = "Parancerous area"))

    iv.obj <- integrate.obj
    names(iv.obj@images) <- NULL
    iv.obj <- iv.obj[, iv.obj$level == "IV"]
    iv.obj <- iv.obj[, iv.obj$sample == "22F-10823-3"]
    densly.cells <- colnames(
        subset(iv.obj, idents = "Tumor cell densely populated area"))

    markers <- FindMarkers(
        integrate.obj,
        logfc.threshold = 0.1,
        ident.1 = densly.cells,
        ident.2 = para.cells,
        verbose = FALSE
    )

    save.dir <- fs::path(save.dirs[["int"]], "IV22F-10823-3密集vsIV癌旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "差异表格.csv"))

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(
        integrate.obj, features = tops, cells = c(densly.cells, para.cells))
    ggsave(fs::path(save.dir, "差异热图.pdf"), p)
    return(markers)
}
markers.list4.6 <- read.csv(
    fs::path(save.dirs[["int"]], "IV22F-10823-3密集vsIV癌旁", "差异表格.csv"),
    row.names = 1
)
markers.list4.6 <- differentialExpression4.6(integrate.obj)

# %% set list
acquireSetList <- function() {
    set.dir <- fs::path(WORKDIR, "results", "integrate", "交集基因")
    # set 1: IV 密集 vs IV 癌旁
    set1.up <- intersect(
        rownames(filter(markers.list4.5, p_val_adj < 0.05, avg_log2FC > 0.5)),
        rownames(filter(markers.list4.6, p_val_adj < 0.05, avg_log2FC > 0.5))
    )
    set1.down <- intersect(
        rownames(filter(markers.list4.5, p_val_adj < 0.05, avg_log2FC < -0.5)),
        rownames(filter(markers.list4.6, p_val_adj < 0.05, avg_log2FC < -0.5))
    )
    save.dir <- fs::path(set.dir, "集合1IV密集vsIV瘤旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(set1.up, fs::path(save.dir, "集合1IV密集vsIV瘤旁上调.csv"))
    write.csv(set1.down, fs::path(save.dir, "集合1IV密集vsIV瘤旁下调.csv"))
    enrichmentGenelist(set1.up, fs::path(set.dir, "集合1IV密集vsIV瘤旁上调"))
    enrichmentGenelist(set1.down, fs::path(set.dir, "集合1IV密集vsIV瘤旁下调"))

    # set 2: GBM 密集 vs IV 癌旁
    set2.up <- intersect(
        rownames(filter(markers.list4.2, p_val_adj < 0.05, avg_log2FC > 0.5)),
        rownames(filter(markers.list4.3, p_val_adj < 0.05, avg_log2FC > 0.5))
    )
    set2.down <- intersect(
        rownames(filter(markers.list4.2, p_val_adj < 0.05, avg_log2FC < -0.5)),
        rownames(filter(markers.list4.3, p_val_adj < 0.05, avg_log2FC < -0.5))
    )
    save.dir <- fs::path(set.dir, "集合2GBM密集vsIV瘤旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(set2.up, fs::path(save.dir, "集合2GBM密集vsIV瘤旁上调.csv"))
    write.csv(set2.down, fs::path(save.dir, "集合2GBM密集vsIV瘤旁下调.csv"))
    enrichmentGenelist(set2.up, fs::path(set.dir, "集合2GBM密集vsIV瘤旁上调"))
    enrichmentGenelist(set2.down, fs::path(set.dir, "集合2GBM密集vsIV瘤旁下调"))

    # set 3: set1 && set2
    set3.up <- intersect(set1.up, set2.up)
    set3.down <- intersect(set1.down, set2.down)
    save.dir <- fs::path(set.dir, "集合3集合1+集合2")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(set3.up, fs::path(save.dir, "集合3集合1+集合2上调.csv"))
    write.csv(set3.down, fs::path(save.dir, "集合3集合1+集合2下调.csv"))
    enrichmentGenelist(set3.up, fs::path(set.dir, "集合3集合1+集合2上调"))
    enrichmentGenelist(set3.down, fs::path(set.dir, "集合3集合1+集合2下调"))

    return.list <- list(
        "set1" = list("up" = set1.up, "down" = set1.down),
        "set2" = list("up" = set2.up, "down" = set2.down),
        "set3" = list("up" = set3.up, "down" = set3.down),
    )
    return(return.list)
}
set.list <- acquireSetList()
