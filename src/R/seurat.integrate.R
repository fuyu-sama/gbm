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
source(fs::path(WORKDIR, "src", "R", "utils.R"))

save.dirs <- lapply(
    idx.list,
    function(idx) save.dir <- fs::path(WORKDIR, "results", idx)
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

# %% integrate data
loadSeuratList()
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

    return(integrate.obj)
}
integrate.obj <- integrateData(seurat.list)

# %% annotation
regionAnnotation <- function(seurat.obj) {
    region.annotations <- list(
        "21B-603-5.4" = "Blood vessel rich area",
        "22F-10823-3.2" = "Blood vessel rich area",
        "22F-21576-1.0" = "Blood vessel rich area",
        "22F-23738-2.2" = "Blood vessel rich area",
        "21B-603-5.2" = "Tumor cell densely populated area",
        "22F-10823-3.0" = "Tumor cell densely populated area",
        "22F-21576-1.2" = "Tumor cell densely populated area",
        "22F-23738-2.0" = "Tumor cell densely populated area",
        "21B-603-5.0" = "Tumor area 1",
        "22F-10823-3.5" = "Tumor area 2",
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
    p <- SpatialDimPlot(seurat.obj, label.size = 3, alpha = 0.6, label = TRUE)
    save.path <- fs::path(save.dirs[["int"]], "分区.pdf")
    ggsave(save.path, p, width = 28, height = 7)

    p1 <- DimPlot(
        seurat.obj,
        reduction = "tsne",
        cols = DiscretePalette(length(region.annotations))
    )
    p2 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "sample")
    ggsave(fs::path(save.dirs[["int"]], "integrate.tsne.pdf"), p1 + p2, width = 20)
    return(seurat.obj)
}
integrate.obj <- regionAnnotation(integrate.obj)

# %% integrate de 1
differentialExpressionInt1 <- function(seurat.obj) {
    markers <- FindAllMarkers(
        seurat.obj, logfc.threshold = 0.1, verbose = FALSE
    )
    save.path <- fs::path(save.dirs[["int"]], "整合分区域差异表达.csv")
    write.csv(markers, save.path)

    tops <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    p <- DoHeatmap(seurat.obj, features = tops$gene) + NoLegend()
    save.path <- fs::path(save.dirs[["int"]], "整合分区域差异表达热图.pdf")
    ggsave(save.path, p, height = 30, width = 15)

    return(markers)
}
markers.list1 <- differentialExpressionInt1(integrate.obj)

# %% enrichment 1
enrichmentInt1 <- function() {
    clusters <- unique(markers.list1$cluster)
    for (cluster in clusters) {
        upscale.table <- markers.list1 %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC > 0)
        downscale.table <- markers.list1 %>%
            group_by(cluster) %>%
            filter(p_val_adj < 0.05, avg_log2FC < 0)
        upscale.genes <- upscale.table[upscale.table$cluster == cluster, ]$gene
        upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
        downscale.genes <- downscale.table[downscale.table$cluster == cluster, ]$gene
        downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]

        save.dir <- fs::path(save.dirs[["int"]], "GO")
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
            fs::path(save.dir, paste(cluster, "GO.pdf", sep = ".")),
            p1 + p2,
            width = 14,
            height = 15
        )
        write.xlsx2(
            as.data.frame(upscale.ego),
            fs::path(save.dir, paste(cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC > 0"
        )
        write.xlsx2(
            as.data.frame(downscale.ego),
            fs::path(save.dir, paste(cluster, "GO.xlsx", sep = ".")),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )

        save.dir <- fs::path(save.dirs[["int"]], "KEGG")
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
                fs::path(save.dir, paste(cluster, "KEGG.xlsx", sep = ".")),
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
                fs::path(save.dir, paste(cluster, "KEGG.xlsx", sep = ".")),
                sheetName = "avg_log2FC < 0",
                append = TRUE
            )
        } else {
            p2 <- NULL
        }
        ggsave(
            fs::path(save.dir, paste(cluster, "KEGG.pdf", sep = ".")),
            p1 + p2,
            width = 14
        )
    }
}
enrichmentInt1()

# %% integrate de 2 tumor vs. para
differentialExpressionInt2 <- function(seurat.obj) {
    ident.1 <- c(
        "Blood vessel rich area",
        "Tumor cell densely populated area",
        "Tumor area 1",
        "Tumor area 2",
        "Tumor area 3",
        "Tumor area 4",
        "Tumor area 5",
        "Tumor area 6",
        "Tumor area 7",
        "Tumor area 8"
    )
    markers <- FindMarkers(
        seurat.obj,
        ident.1 = ident.1,
        ident.2 = c("Junction area", "Junction area 2", "Parancerous area"),
        logfc.threshold = 0.1,
        verbose = FALSE
    )
    save.dir <- fs::path(save.dirs[["int"]], "肿瘤vs瘤旁")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    save.path <- fs::path(save.dir, "差异表达.csv")
    write.csv(markers, save.path)

    tops <- markers %>% top_n(n = 50, wt = avg_log2FC) %>% rownames()
    bottoms <- markers %>% top_n(n = -50, wt = avg_log2FC) %>% rownames()
    p <- DoHeatmap(seurat.obj, features = c(tops, bottoms)) + NoLegend()
    save.path <- fs::path(save.dir, "差异表达热图.pdf")
    ggsave(save.path, p, height = 20, width = 15)

    return(markers)
}
markers.list2 <- differentialExpressionInt2(integrate.obj)

# %% entichment 2
enrichmentInt2 <- function() {
    upscale.table <- markers.list1 %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- markers.list1 %>%
        filter(p_val_adj < 0.05, avg_log2FC < 0)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]

    save.dir <- fs::path(save.dirs[["int"]], "肿瘤vs瘤旁", "GO")
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
        fs::path(save.dir, "GO.pdf"),
        p1 + p2,
        width = 14,
        height = 15
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir, "GO.xlsx"),
        sheetName = "avg_log2FC > 0"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir, "GO.xlsx"),
        sheetName = "avg_log2FC < 0",
        append = TRUE
    )

    save.dir <- fs::path(save.dirs[["int"]], "肿瘤vs瘤旁", "KEGG")
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
            fs::path(save.dir, "KEGG.xlsx"),
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
            fs::path(save.dir, "KEGG.xlsx"),
            sheetName = "avg_log2FC < 0",
            append = TRUE
        )
    } else {
        p2 <- NULL
    }
    ggsave(
        fs::path(save.dir, "KEGG.pdf"),
        p1 + p2,
        width = 14
    )
}
enrichmentInt2()

# %% integrate de 3 IV vs GBM
differentialExpressionInt3 <- function(compare.ident) {
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
    file.name <- paste(compare.ident, "四级.vs.GBM", sep = ".")
    save.dir <- fs::path(save.dirs[["int"]], "四级-vs-gbm")
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
integrate.markers3 <- mclapply(named.regions, differentialExpressionInt3)

# %% enrichment for integrate 3
enrichmentInt3 <- function(region) {
    upscale.table <- integrate.markers3[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers3[[region]] %>%
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
    save.dir <- fs::path(save.dirs[["int"]], "四级-vs-GBM", "GO")
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

    save.dir <- fs::path(save.dirs[["int"]], "四级-vs-GBM", "KEGG")
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
sapply(names(integrate.markers3), enrichmentInt3)

# %% integrate de 4
differentialExpressionInt4 <- function(compare.ident) {
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
    save.dir <- fs::path(save.dirs[["int"]], "其他区域-vs-瘤旁")
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
integrate.markers4 <- mclapply(named.regions, differentialExpressionInt4)

# %% enrichment for integrate 4
enrichmentInt4 <- function(region) {
    upscale.table <- integrate.markers4[[region]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 0)
    downscale.table <- integrate.markers4[[region]] %>%
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
    save.dir <- fs::path(save.dirs[["int"]], "其他区域-vs-瘤旁")
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
sapply(names(integrate.markers4), enrichmentInt4)
