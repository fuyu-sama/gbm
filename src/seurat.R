#! /usr/bin/env Rscript

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
