#! /usr/bin/env Rscript

# %% environment config
library <- function(...) suppressMessages(base::library(...))
library('parallel')
library('dplyr')
library('xlsx')
library('ggplot2')
library('ggvenn')
library('Seurat')
library('clusterProfiler')
library('ComplexHeatmap')
library('org.Hs.eg.db')
library('ggsci')

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
source(fs::path(WORKDIR, "src", "utils.R"))

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
        seurat.obj <- RenameCells(seurat.obj, add.cell.id = idx)
        names(seurat.obj@images) <- idx
        seurat.obj@meta.data <- seurat.obj@meta.data[, 1:5]
        seurat.obj@meta.data["sample"] <- rep(idx, dim(seurat.obj)[2])
        seurat.obj@meta.data["level"] <- rep(cons[[idx]], dim(seurat.obj)[2])
        Idents(seurat.obj) <- paste(idx, Idents(seurat.obj), sep = ".")
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
        reduction = "rpca",
        k.anchor = 4,
        dims = 1:30,
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
write.csv(
    t(GetAssayData(integrate.obj, assay = "SCT", slot = "count")),
    fs::path(WORKDIR, "Data", "counts", "full-sct.csv")
)
saveIntegrate()

# %% annotation
integrate.obj@images[[3]]@spot.radius <- 0.00958
integrate.obj@images[[4]]@spot.radius <- 0.00958
regionAnnotation <- function(seurat.obj) {
    region.annotations <- list(
        "21B-603-5.4" = "Blood vessel rich area",
        "22F-21576-1.3" = "Blood vessel rich area",
        "22F-23738-2.3" = "Blood vessel rich area",
        "21B-603-5.2" = "IV Tumor cell densely populated area",
        "22F-10823-3.5" = "IV Tumor cell densely populated area",
        "22F-21576-1.2" = "GBM Tumor cell densely populated area",
        "22F-23738-2.0" = "GBM Tumor cell densely populated area",
        "22F-23738-2.4" = "GBM Tumor cell densely populated area",
        "21B-603-5.0" = "IV Tumor area 1",
        "22F-10823-3.0" = "IV Tumor area 2",
        "22F-10823-3.2" = "IV Tumor area 3",
        "22F-21576-1.0" = "GBM Tumor area 1",
        "22F-21576-1.1" = "GBM Tumor area 2",
        "22F-21576-1.4" = "GBM Tumor area 3",
        "22F-23738-2.1" = "GBM Tumor area 4",
        "22F-23738-2.2" = "GBM Tumor area 5",
        "22F-23738-2.5" = "GBM Tumor area 6",
        "22F-23738-2.6" = "GBM Tumor area 7",
        "22F-23738-2.7" = "GBM Tumor area 8",
        "21B-603-5.1" = "Junction area",
        "22F-10823-3.1" = "Junction area",
        "22F-10823-3.3" = "Junction area",
        "21B-603-5.3" = "Normal tissue adjacent to tumor area",
        "22F-10823-3.4" = "Normal tissue adjacent to tumor area"
    )
    seurat.obj@meta.data$orig.idents <- Idents(seurat.obj)
    seurat.obj <- RenameIdents(seurat.obj, region.annotations)

    p <- SpatialDimPlot(
        seurat.obj, label.size = 3, alpha = 1, cols = cols, label = FALSE
    )
    save.path <- fs::path(save.dirs[["int"]], "分区.pdf")
    ggsave(save.path, p, width = 28, height = 7)

    p1 <- DimPlot(seurat.obj, reduction = "tsne", cols = cols)
    p2 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "sample")
    p3 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "level")
    ggsave(
        fs::path(save.dirs[["int"]], "integrate.tsne.pdf"),
        p1 + p3 + p2,
        width = 30
    )

    p <- DimPlot(
        seurat.obj, reduction = "tsne", cols = cols, split.by = "sample")
    ggsave(
        fs::path(save.dirs[["int"]], "integrate.tsne.split.pdf"),
        p,
        width = 40
    )

    write.csv(Idents(seurat.obj), fs::path(WORKDIR, "results", "idents.csv"))

    return(seurat.obj)
}
integrate.obj <- regionAnnotation(integrate.obj)

# %% SCENIC
runSCENIC <- function() {
    system(paste(fs::path(WORKDIR, "src", "scenic.sh"), "full-sct", sep = " "))
    auc.df <- read.csv(
        fs::path(WORKDIR, "results", "scenic", "full-sct.auc.csv"),
        header = TRUE,
        row.names = 1,
        check.names = FALSE
        ) %>% t()
    auc.matrix <- as.matrix(auc.df)
    colnames(auc.matrix) <- colnames(auc.df)
    rownames(auc.matrix) <- rownames(auc.df)
    auc.assay <- CreateAssayObject(data = auc.matrix)
    integrate.obj[["SCENIC"]] <- auc.assay

    bin.df <- read.csv(
        fs::path(WORKDIR, "results", "scenic", "full-sct.bin.csv"),
        header = TRUE,
        row.names = 1,
        check.names = FALSE
        ) %>% t()
    bin.matrix <- as.matrix(bin.df)
    colnames(bin.matrix) <- colnames(bin.df)
    rownames(bin.matrix) <- rownames(bin.df)
    bin.assay <- CreateAssayObject(data = bin.matrix)
    integrate.obj[["SCENIC.bin"]] <- bin.assay

    rss.level.df <- read.csv(
        fs::path(WORKDIR, "results", "scenic", "full-sct.rss_level.csv"),
        header = TRUE,
        row.names = 1,
        check.names = FALSE
        ) %>% t() %>% as.data.frame()
    rss.level.df[is.na(rss.level.df)] <- 0

    top.gbm <- rss.level.df %>%
        arrange(desc(GBM)) %>%
        top_n(n = 50, wt = GBM) %>%
        rownames() %>%
        gsub("\\(.*\\)", "", .)
    top.gbm.plus <- paste0(top.gbm, "(+)")

    top.iv <- rss.level.df %>%
        arrange(desc(IV)) %>%
        top_n(n = 50, wt = IV) %>%
        rownames() %>%
        gsub("\\(.*\\)", "", .)
    top.iv.plus <- paste0(top.iv, "(+)")

    auc.obj <- integrate.obj
    DefaultAssay(auc.obj) <- "SCENIC"

    p <- FeaturePlot(
        auc.obj,
        features = top.gbm.plus,
        reduction = "tsne",
        min.cutoff = 0,
        max.cutoff = "q90",
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "auc-gbm-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    p <- FeaturePlot(
        auc.obj,
        features = top.iv.plus,
        reduction = "tsne",
        min.cutoff = 0,
        max.cutoff = "q90",
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "auc-iv-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    auc.obj <- integrate.obj
    DefaultAssay(auc.obj) <- "SCENIC.bin"
    p <- FeaturePlot(
        auc.obj,
        features = top.gbm.plus,
        reduction = "tsne",
        cols = c("lightgrey", "red"),
        min.cutoff = 0,
        max.cutoff = 1,
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "bin-gbm-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    p <- FeaturePlot(
        auc.obj,
        features = top.iv.plus,
        reduction = "tsne",
        cols = c("lightgrey", "red"),
        min.cutoff = 0,
        max.cutoff = 1,
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "bin-iv-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    p <- FeaturePlot(
        integrate.obj,
        features = top.gbm,
        reduction = "tsne",
        min.cutoff = 0,
        max.cutoff = "q90",
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "expression-gbm-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    p <- FeaturePlot(
        integrate.obj,
        features = top.iv,
        reduction = "tsne",
        min.cutoff = 0,
        max.cutoff = "q90",
        ncol = 5
    )
    ggsave(
        fs::path(WORKDIR, "results", "scenic-plot", "expression-iv-tsne.pdf"),
        p,
        width = 25,
        height = 50,
        limitsize = FALSE
    )

    threshold.df <- read.csv(
        fs::path(WORKDIR, "results", "scenic", "full-sct.threshold.csv"),
        row.names = 1
    )

    save.dir <- fs::path(WORKDIR, "results", "scenic-plot", "iv")
    if (!fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    for (gene in top.iv.plus) {
        p <- histPlot(integrate.obj, gene, threshold.df)
        ggsave(fs::path(save.dir, paste0(gene, ".pdf")))
    }

    save.dir <- fs::path(WORKDIR, "results", "scenic-plot", "gbm")
    if (!fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    for (gene in top.gbm.plus) {
        p <- histPlot(integrate.obj, gene, threshold.df)
        ggsave(fs::path(save.dir, paste0(gene, ".pdf")))
    }

    top100.gbm <- rss.level.df %>%
        arrange(desc(GBM)) %>%
        top_n(n = 100, wt = GBM) %>%
        rownames()

    top100.iv <- rss.level.df %>%
        arrange(desc(IV)) %>%
        top_n(n = 100, wt = IV) %>%
        rownames()

    p <- ggvenn(
        list("GBM Top 100 Regulon" = top100.gbm, "IV Top 100 Regulon" = top100.iv),
        c("GBM Top 100 Regulon", "IV Top 100 Regulon")
    )
    ggsave(fs::path(WORKDIR, "results", "scenic-plot", "venn.pdf"), p)

    except.spots <- c(
        colnames(subset(auc.obj, idents = "Normal tissue adjacent to tumor area")),
        colnames(subset(auc.obj, idents = "Junction area")),
        colnames(subset(auc.obj, idents = "Blood vessel rich area"))
    )
    all.spots <- colnames(auc.obj)
    tumor.spots <- all.spots[! all.spots %in% except.spots]
    level <- as.vector(auc.obj[, tumor.spots]$level)
    names(level) <- tumor.spots

    tops <- c(
        intersect(top100.gbm, top100.iv),
        setdiff(top100.gbm, top100.iv),
        setdiff(top100.iv, top100.gbm)
    )

    pdf(fs::path(WORKDIR, "results", "scenic-plot", "heatmap.pdf"), width = 21)
    ha <- rowAnnotation(
        level = level, col = list(level = c("IV" = "#ffff7f", "GBM" = "#7f7fff")))
    ht <- Heatmap(
        scale(t(auc.df[tops, tumor.spots])),
        left_annotation = ha,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE
    )
    draw(ht)
    dev.off()
}
runSCENIC()

# %% draw RCTD
drawRCTD <- function() {
    FLAG <- 1
    for (idx in idx.full) {
        if (FLAG) {
            rctd.results <- read.csv(
                fs::path(save.dirs[[idx]], paste0(idx, ".rctd.csv")),
                row.names = 1,
                header = TRUE
            )
            rownames(rctd.results) <- paste(idx, rownames(rctd.results), sep = "_")
            FLAG <- 0
        } else {
            read.results <- read.csv(
                fs::path(save.dirs[[idx]], paste0(idx, ".rctd.csv")),
                row.names = 1,
                header = TRUE
            )
            rownames(read.results) <- paste(idx, rownames(read.results), sep = "_")
            rctd.results <- rbind(rctd.results, read.results)
        }
    }

    rctd.results <- rctd.results[colnames(integrate.obj), ]
    rctd.results[is.na(rctd.results)] <- 0
    integrate.obj$Immune <- rctd.results$Immune
    integrate.obj$NormalBrain <- rctd.results$NormalBrain
    integrate.obj$Tumor <- rctd.results$Tumour
    integrate.obj$first_type <- rctd.results$first_type
    integrate.obj$region <- Idents(integrate.obj)

    p <- FeaturePlot(
        integrate.obj,
        features = c("NormalBrain", "Tumor", "Immune"),
        max.cutoff = "q90",
        ncol = 3,
        reduction = "tsne"
    )
    ggsave(fs::path(WORKDIR, "results", "rctd.tsne.pdf"), width = 21, height = 7)

    p <- RidgePlot(
        integrate.obj,
        features = c("NormalBrain", "Tumor", "Immune"),
        cols = cols,
        ncol = 3
    )
    ggsave(fs::path(WORKDIR, "results", "rctd.ridge.pdf"), width = 21, height = 10)

    stack.df <- integrate.obj@meta.data[, c("first_type", "region")]
    stack.df <- stack.df[stack.df$first_type != 0, ]
    stack.df$number <- 1
    stack.df <- plyr::ddply(
        stack.df, "region", transform, percent = 1 / sum(number) * 100)
    p <- ggplot(stack.df, aes(region, percent, fill = first_type)) +
        geom_bar(stat = "identity", position = "stack") +
        coord_flip() +
        theme_bw()
    ggsave(fs::path(WORKDIR, "results", "rctd.stack.pdf"), p)
}
drawRCTD()

# %% gene list
runModuleScore <- function() {
    gene.list <- list(
        "ubiquitin" = loadUbiquitin(),
        "rbp" = loadRBP(),
        "kinase" = loadKinase()
    )
    integrate.obj <- AddModuleScore(
        integrate.obj,
        gene.list,
        assay = "integrated",
        slot = "data",
        name = names(gene.list),
        search = TRUE,
        verbose = FALSE,
        timeout = 30
    )
    p <- FeaturePlot(
        integrate.obj,
        features = c("ubiquitin1", "rbp2", "kinase3"),
        reduction = "tsne",
        max.cutoff = "q90",
        min.cutoff = "q5",
        cols = c("lightgrey", "darkgreen")
    )
    ggsave(
        fs::path(WORKDIR, "results", "module-score.pdf"),
        p, width = 14, height = 14
    )

    p <- SpatialFeaturePlot(
        integrate.obj,
        features = c("ubiquitin1", "rbp2", "kinase3")
    )
    ggsave(
        fs::path(WORKDIR, "results", "module-score.2.pdf"),
        p, width = 14, height = 14
    )

    p <- DotPlot(
        integrate.obj,
        features = c("ubiquitin1", "rbp2", "kinase3"),
        scale = FALSE
    )
    p <- ggplot(p$data, aes(x = id, y = features.plot, fill = avg.exp.scaled)) +
        geom_tile() +
        coord_flip() +
        theme(panel.grid = element_blank(), text = element_text(size = 20)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_fill_gsea()
    ggsave(fs::path(WORKDIR, "results", "module-score.hm.pdf"), p, width = 8)
}

# %% tumor vs nat
differentialExpressionNAT <- function() {
    markers.tumor.GBM <- FindMarkers(
        integrate.obj,
        ident.1 = "GBM Tumor cell densely populated area",
        ident.2 = "Normal tissue adjacent to tumor area"
    )
    markers.tumor.IV <- FindMarkers(
        integrate.obj,
        ident.1 = "IV Tumor cell densely populated area",
        ident.2 = "Normal tissue adjacent to tumor area"
    )
    markers.tumor.IV$gene <- rownames(markers.tumor.IV)
    markers.tumor.GBM$gene <- rownames(markers.tumor.GBM)

    save.dir <- fs::path(WORKDIR, "results", "tumor.vs.nat")
    if (!fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers.tumor.IV, fs::path(save.dir, "iv.de.csv"))
    write.csv(markers.tumor.GBM, fs::path(save.dir, "gbm.de.csv"))

    return(list(IV = markers.tumor.IV, GBM = markers.tumor.GBM))
}
markers.tumor <- differentialExpressionNAT()

# %% tumor vs nat enrichment
enrichmentFindMarkers(
    markers.tumor[["IV"]],
    fs::path(WORKDIR, "results", "tumor.vs.nat", "IV")
)
enrichmentFindMarkers(
    markers.tumor[["GBM"]],
    fs::path(WORKDIR, "results", "tumor.vs.nat", "GBM")
)

# %% draw tumor vs nat enrichment
drawTumorNATKEGG <- function() {
    read.path <- fs::path(
        WORKDIR, "results", "tumor.vs.nat", "IV", "KEGG", "KEGG.xlsx"
    )
    iv.enrich <- read.xlsx2(read.path, 1)
    rownames(iv.enrich) <- iv.enrich$ID
    iv.enrich <- iv.enrich[, 3:dim(iv.enrich)[2]]
    iv.enrich$Count <- as.numeric(iv.enrich$Count)
    read.path <- fs::path(
        WORKDIR, "results", "tumor.vs.nat", "GBM", "KEGG", "KEGG.xlsx"
    )
    gbm.enrich <- read.xlsx2(read.path, 1)
    rownames(gbm.enrich) <- gbm.enrich$ID
    gbm.enrich <- gbm.enrich[, 3:dim(gbm.enrich)[2]]
    gbm.enrich$Count <- as.numeric(gbm.enrich$Count)

    enrich.list <- list(
        "IV" = rownames(iv.enrich),
        "GBM" = rownames(gbm.enrich)
    )
    overlaps <- gplots::venn(enrich.list, show.plot = FALSE)
    overlaps <- attributes(overlaps)$intersection

    draw.list <- list(
        "IDH mutant" = iv.enrich[overlaps[["IV:GBM"]], ] %>% top_n(10, wt = Count),
        "IDH wildtype" = gbm.enrich[overlaps[["IV:GBM"]], ] %>% top_n(10, wt = Count)
    )
    p <- ggenrich2(draw.list, overlap = TRUE) + ggtitle("KEGG Enrichment")
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "IV:GBM-common-KEGG.pdf"), p)

    pathways <- c("hsa05166", "hsa04810", "hsa05205", "hsa04015", "hsa05203")
    draw.list <- list(
        "IDH mutant" = iv.enrich[pathways, ],
        "IDH wlidtype" = gbm.enrich[overlaps[["GBM"]], ]
    )
    p <- ggenrich2(draw.list, overlap = FALSE) + ggtitle("KEGG Enrichment")
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "IV:GBM-iden-KEGG.pdf"), p)
}

drawTumorNATGO <- function() {
    read.path <- fs::path(
        WORKDIR, "results", "tumor.vs.nat", "IV", "GO", "GO.xlsx"
    )
    iv.enrich <- read.xlsx2(read.path, 1) %>% filter(ONTOLOGY == "BP")
    rownames(iv.enrich) <- iv.enrich$ID
    iv.enrich <- iv.enrich[, 3:dim(iv.enrich)[2]]
    iv.enrich$Count <- as.numeric(iv.enrich$Count)
    read.path <- fs::path(
        WORKDIR, "results", "tumor.vs.nat", "GBM", "GO", "GO.xlsx"
    )
    gbm.enrich <- read.xlsx2(read.path, 1) %>% filter(ONTOLOGY == "BP")
    rownames(gbm.enrich) <- gbm.enrich$ID
    gbm.enrich <- gbm.enrich[, 3:dim(gbm.enrich)[2]]
    gbm.enrich$Count <- as.numeric(gbm.enrich$Count)

    enrich.list <- list(
        "IV" = rownames(iv.enrich),
        "GBM" = rownames(gbm.enrich)
    )
    overlaps <- gplots::venn(enrich.list, show.plot = FALSE)
    overlaps <- attributes(overlaps)$intersection

    draw.list <- list(
        "IDH mutant" = iv.enrich[overlaps[["IV:GBM"]], ] %>% top_n(10, wt = Count),
        "IDH wildtype" = gbm.enrich[overlaps[["IV:GBM"]], ] %>% top_n(10, wt = Count)
    )
    p <- ggenrich2(draw.list, overlap = TRUE) + ggtitle("GO Enrichment")
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "IV:GBM-common-GO.pdf"), p)

    draw.list <- list(
        "IDH mutant" = iv.enrich[overlaps[["IV"]], ] %>% top_n(5, wt = Count),
        "IDH wlidtype" = gbm.enrich[overlaps[["GBM"]], ] %>% top_n(5, wt = Count)
    )
    p <- ggenrich2(draw.list, overlap = FALSE) + ggtitle("GO Enrichment")
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "IV:GBM-iden-GO.pdf"), p)
}

drawTumorNATKEGG()
drawTumorNATGO()

# %% draw tumor vs nat
drawTumorNAT <- function() {
    th <- theme(text = element_text(size = 24))
    cells <- c(
        subset(
            integrate.obj, idents = "GBM Tumor cell densely populated area"
            ) %>% colnames(),
        subset(
            integrate.obj, idents = "IV Tumor cell densely populated area"
            ) %>% colnames(),
        subset(
            integrate.obj, idents = "Normal tissue adjacent to tumor area"
            ) %>% colnames()
    )

    genes <- c(
        "SOD2", "UBA52", "S100A6", "CTSB", "PPP1CB",
        "SMOC1", "APOE", "HIPK2",
        "SPP1", "IGFBP2", "CALD1", "TMSB4X"
    )
    p <- FeaturePlot(
        integrate.obj,
        features = genes,
        cells = cells,
        reduction = "tsne",
        min.cutoff = 0,
        max.cutoff = "q90",
        ncol = 4
        ) + NoLegend()
    save.path <- fs::path(WORKDIR, "results", "tumor.vs.nat", "tsne.pdf")
    ggsave(save.path, p, width = 28, height = 21)

    p <- ggvolcano(
        markers.tumor[["IV"]], genes = genes, label_size = 5
        ) + th #+ coord_flip()
    save.path <- fs::path(WORKDIR, "results", "tumor.vs.nat", "iv-volcano.pdf")
    ggsave(save.path, p, width = 7)

    p <- ggvolcano(
        markers.tumor[["GBM"]], genes = genes, label_size = 5
        ) + th #+ coord_flip()
    save.path <- fs::path(WORKDIR, "results", "tumor.vs.nat", "gbm-volcano.pdf")
    ggsave(save.path, p, width = 7)
}
drawTumorNAT()

# %%
vennGBMIV <- function() {
    tops.gbm <- markers.tumor[["GBM"]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        filter(pct.1 > pct.2) %>%
        arrange(desc(avg_log2FC)) %>%
        rownames()
    tops.iv <- markers.tumor[["IV"]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        filter(pct.1 > pct.2) %>%
        arrange(desc(avg_log2FC)) %>%
        rownames()
    tops.list <- list(
        "IDH wildtype" = tops.gbm,
        "IDH mutant" = tops.iv
    )
    overlaps <- gplots::venn(tops.list, show.plot = FALSE)
    overlaps <- attributes(overlaps)$intersection
    p <- ggvenn(tops.list, columns = names(tops.list))
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "venn.pdf"), p)
}
vennGBMIV()

# %%
differentialExpressionGBMvsIV <- function() {
    seurat.obj <- integrate.obj
    names(seurat.obj@images) <- NULL
    except.spots <- c(
        colnames(subset(seurat.obj, idents = "Normal tissue adjacent to tumor area")),
        colnames(subset(seurat.obj, idents = "Junction area")),
        colnames(subset(seurat.obj, idents = "Blood vessel rich area"))
    )
    all.spots <- colnames(seurat.obj)
    tumor.spots <- all.spots[! all.spots %in% except.spots]
    seurat.obj <- seurat.obj[, tumor.spots]
    markers <- FindMarkers(
        seurat.obj,
        min.pct = 0.3,
        ident.1 = "GBM",
        ident.2 = "IV",
        group.by = "level"
    )
    save.dir <- fs::path(WORKDIR, "results", "gbm.vs.iv")
    if (!fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    write.csv(markers, fs::path(save.dir, "markers.csv"))
    return(markers)
}
markers.GBMvsIV <- differentialExpressionGBMvsIV() %>%
    filter(p_val_adj < 0.05) %>%
    filter(pct.1 > 0.3, pct.2 > 0.3)

# %%
enrichmentFindMarkers(
    markers.GBMvsIV,
    fs::path(WORKDIR, "results", "gbm.vs.iv")
)

# %%
drawGBMvsIV <- function() {
    tops <- top_n(markers.GBMvsIV, n = 100, wt = avg_log2FC) %>% rownames()
    bottoms <- top_n(markers.GBMvsIV, n = -100, wt = avg_log2FC) %>% rownames()
    genes <- c(tops, bottoms)

    seurat.obj <- integrate.obj
    names(seurat.obj@images) <- NULL
    except.spots <- c(
        colnames(subset(seurat.obj, idents = "Normal tissue adjacent to tumor area")),
        colnames(subset(seurat.obj, idents = "Junction area")),
        colnames(subset(seurat.obj, idents = "Blood vessel rich area"))
    )
    all.spots <- colnames(seurat.obj)
    tumor.spots <- all.spots[! all.spots %in% except.spots]
    seurat.obj <- seurat.obj[, tumor.spots]

    p <- DoHeatmap(
        seurat.obj,
        features = genes,
        cells = tumor.spots,
        group.by = "level"
    )
    ggsave(fs::path(WORKDIR, "results", "gbm.vs.iv", "heatmap.pdf"), p)

}
drawGBMvsIV()
