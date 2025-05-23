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
        "kinase" = loadKinase(),
        "phosphatase" = loadDEPOD()
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
        features = c("ubiquitin1", "rbp2", "kinase3", "phosphatase4"),
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
        features = c("ubiquitin1", "rbp2", "kinase3", "phosphatase4")
    )
    ggsave(
        fs::path(WORKDIR, "results", "module-score.2.pdf"),
        p, width = 14, height = 14
    )

    p <- DotPlot(
        integrate.obj,
        features = c("ubiquitin1", "rbp2", "kinase3", "phosphatase4"),
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
runModuleScore()

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
    fs::path(WORKDIR, "results", "tumor.vs.nat", "IV"),
    logfc.threshold = 1
)
enrichmentFindMarkers(
    markers.tumor[["GBM"]],
    fs::path(WORKDIR, "results", "tumor.vs.nat", "GBM"),
    logfc.threshold = 1
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

    #pathways <- c("hsa05166", "hsa04810", "hsa05205", "hsa04015", "hsa05203")
    pathways1 <- c("hsa04010", "hsa04530", "hsa04621", "hsa04141", "hsa04630")
    pathways2 <- c("hsa04926", "hsa05150", "hsa01522", "hsa04979", "hsa05218")
    draw.list <- list(
        "IDH mutant" = iv.enrich[pathways1, ],
        #"IDH wlidtype" = gbm.enrich[overlaps[["GBM"]], ]
        "IDH wlidtype" = gbm.enrich[pathways2, ]
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
    ggsave(fs::path(WORKDIR, "results", "tumor.vs.nat", "IV:GBM-iden-GO.pdf"), p, width = 8)
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
        markers.tumor[["IV"]], genes = genes, label_size = 5, cut_off_logFC = 1
        ) + th #+ coord_flip()
    save.path <- fs::path(WORKDIR, "results", "tumor.vs.nat", "iv-volcano.pdf")
    ggsave(save.path, p, width = 7)

    p <- ggvolcano(
        markers.tumor[["GBM"]], genes = genes, label_size = 5, cut_off_logFC = 1
        ) + th #+ coord_flip()
    save.path <- fs::path(WORKDIR, "results", "tumor.vs.nat", "gbm-volcano.pdf")
    ggsave(save.path, p, width = 7)
}
drawTumorNAT()

# %% draw venn diagram
vennGBMIV <- function() {
    tops.gbm <- markers.tumor[["GBM"]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
        filter(pct.1 > pct.2) %>%
        arrange(desc(avg_log2FC)) %>%
        rownames()
    tops.iv <- markers.tumor[["IV"]] %>%
        filter(p_val_adj < 0.05, avg_log2FC > 1) %>%
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

# %% gbm vs iv
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

# %% enrichment gbm vs iv
enrichmentFindMarkers(
    markers.GBMvsIV,
    fs::path(WORKDIR, "results", "gbm.vs.iv")
)

# %% rctd 2
rctd.list <- list()
for (idx in names(integrate.obj@images)) {
    results.df <- read.csv(
        fs::path(save.dirs[[idx]], paste0(idx, ".rctd.2.csv")),
        header = TRUE,
        row.names = 1
    )
    rownames(results.df) <- paste(idx, results.df$Row.names, sep = "_")
    rctd.list[[idx]] <- results.df
}
rctd.df <- bind_rows(rctd.list)
rctd.df <- rctd.df[colnames(integrate.obj), ]

integrate.obj$Glioma.1 <- rctd.df$Glioma.1
integrate.obj$Glioma.2 <- rctd.df$Glioma.2
integrate.obj$Glioma.3 <- rctd.df$Glioma.3
integrate.obj$Glioma.4 <- rctd.df$Glioma.4
integrate.obj$Glioma.5 <- rctd.df$Glioma.5
integrate.obj$Glioma.6 <- rctd.df$Glioma.6
integrate.obj$Macrophage <- rctd.df$Macrophage
integrate.obj$Oligodendrocyte <- rctd.df$Oligodendrocyte
integrate.obj$T.cells <- rctd.df$T.cells
celltypes <- c(
    "Glioma.1", "Glioma.2", "Glioma.3", "Glioma.4", "Glioma.5",
    "Glioma.6", "Macrophage", "Oligodendrocyte", "T.cells"
)
for (celltype in celltypes) {
    p <- SpatialFeaturePlot(integrate.obj, celltype)
    ggsave(paste0(celltype, ".rctd.pdf"), p, width = 28, height = 7)
}

# %% iv sub 1
integrate.sub <- integrate.obj[, integrate.obj$level == "IV"]
region.annotations <- list(
    "Blood vessel rich area" = "Tumor area",
    "IV Tumor cell densely populated area" = "Tumor area",
    "IV Tumor area 1" = "Tumor area",
    "IV Tumor area 2" = "Tumor area",
    "IV Tumor area 3" = "Tumor area"
)
integrate.sub <- RenameIdents(integrate.sub, region.annotations)
markers.iv <- FindAllMarkers(integrate.sub, min.pct = 0.3)

# %%
cluster.value <- "Tumor area"
enrich.genes <- markers.iv %>%
    filter(cluster == cluster.value) %>%
    filter(avg_log2FC > 1, p_val_adj < 0.05)
enrich.list <- enrichmentGenelist(
    enrich.genes$gene, fs::path(WORKDIR, "results", "250218", cluster.value)
)
draw.pathways <- enrich.list$ego@result %>%
    arrange(desc(Count)) %>%
    top_n(10)
p <- dotplot(enrich.list$ego, showCategory = draw.pathways$Description)
ggsave(fs::path(WORKDIR, "results", "250220", paste0(cluster.value, ".pdf")), p)

cluster.value <- "Normal tissue adjacent to tumor area"
enrich.genes <- markers.iv %>%
    filter(cluster == cluster.value) %>%
    filter(avg_log2FC > 1, p_val_adj < 0.05)
enrich.list <- enrichmentGenelist(
    enrich.genes$gene, fs::path(WORKDIR, "results", "250218", cluster.value)
)
draw.pathways <- enrich.list$ego@result %>%
    arrange(desc(Count)) %>%
    top_n(10)
p <- dotplot(enrich.list$ego, showCategory = draw.pathways$Description)
ggsave(fs::path(WORKDIR, "results", "250220", paste0(cluster.value, ".pdf")), p)

cluster.value <- "Junction area"
enrich.genes <- markers.iv %>%
    filter(cluster == cluster.value) %>%
    filter(avg_log2FC > 1, p_val_adj < 0.05)
enrich.list <- enrichmentGenelist(
    enrich.genes$gene, fs::path(WORKDIR, "results", "250218", cluster.value)
)
draw.pathways <- enrich.list$ego@result %>%
    arrange(desc(Count)) %>%
    top_n(10)
p <- dotplot(enrich.list$ego, showCategory = draw.pathways$Description)
ggsave(fs::path(WORKDIR, "results", "250220", paste0(cluster.value, ".pdf")), p)

# %% iv sub 2
integrate.sub.iv <- integrate.obj[, integrate.obj$level == "IV"]
region.annotations <- list(
    "IV Tumor cell densely populated area" = "Tumor area",
    "IV Tumor area 1" = "Tumor area",
    "IV Tumor area 2" = "Tumor area",
    "IV Tumor area 3" = "Tumor area"
)
integrate.sub.iv <- RenameIdents(integrate.sub.iv, region.annotations)
integrate.sub.iv.levels <- c(
    "Blood vessel rich area",
    "Tumor area",
    "Junction area",
    "Normal tissue adjacent to tumor area"
)
Idents(integrate.sub.iv) <- factor(
    Idents(integrate.sub.iv), levels = integrate.sub.iv.levels
)
#markers.iv.2 <- FindAllMarkers(integrate.sub.iv, min.pct = 0.3)
#write.csv(markers.iv.2, "markers.iv.2.csv")
markers.iv.2 <- read.csv("markers.iv.2.csv", row.names = 1, header = TRUE)

# %%
markers.iv.2$cluster <- factor(
    markers.iv.2$cluster, levels = integrate.sub.iv.levels
)
tops <- markers.iv.2 %>%
    filter(!grepl("^MT-", gene)) %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup()
p <- DoHeatmap(integrate.sub.iv, features = tops$gene) +
    ggtitle("IDH mutant")
write.csv(tops, "results/250221/IV热图表格.csv")
ggsave("results/250221/IV热图.pdf", p, width = 10, height = 10)

# %%
top30 <- markers.iv.2 %>%
    filter(cluster == "Junction area") %>%
    filter(!grepl("^MT-", gene)) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 30)
p <- DotPlot(integrate.sub.iv, features = top30$gene) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("IDH mutant")
ggsave('dotplot.pdf', p, width = 14)

# %%
markers.iv.2.enrichment <- function(cluster.value) {
    enrich.genes <- markers.iv.2 %>%
        filter(cluster == cluster.value) %>%
        filter(avg_log2FC > 1, p_val_adj < 0.05)
    enrich.list <- enrichmentGenelist(
        enrich.genes$gene, fs::path(WORKDIR, "results", "250218", cluster.value)
    )
    draw.pathways <- enrich.list$ego@result %>%
        arrange(desc(Count)) %>%
        top_n(10)
    p <- dotplot(enrich.list$ego, showCategory = draw.pathways$Description)
    ggsave(fs::path(WORKDIR, "results", "250221", paste0(cluster.value, ".pdf")), p)
}
markers.iv.2.enrichment("Normal tissue adjacent to tumor area")
markers.iv.2.enrichment("Tumor area")
markers.iv.2.enrichment("Junction area")

# %% gbm sub
integrate.sub.gbm <- integrate.obj[, integrate.obj$level == "GBM"]
integrate.sub.gbm.levels <- c(
    "Blood vessel rich area",
    "GBM Tumor area 2",
    "GBM Tumor area 3",
    "GBM Tumor area 6",
    "GBM Tumor area 7",
    "GBM Tumor area 8",
    "GBM Tumor area 1",
    "GBM Tumor area 4",
    "GBM Tumor area 5",
    "GBM Tumor cell densely populated area"
)
Idents(integrate.sub.gbm) <- factor(
    Idents(integrate.sub.gbm), levels = integrate.sub.gbm.levels
)
#markers.gbm <- FindAllMarkers(integrate.sub.gbm, min.pct = 0.3)
#write.csv(markers.gbm, "GBM热图.csv")
markers.gbm <- read.csv("GBM热图.csv", row.names = 1, header = TRUE)

# %%
markers.gbm$cluster <- factor(
    markers.gbm$cluster, levels = integrate.sub.gbm.levels
)
tops <- markers.gbm %>%
    filter(!grepl("^MT-", gene)) %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup()
p <- DoHeatmap(integrate.sub.gbm, features = tops$gene) +
    ggtitle("IDH wildtype")
write.csv(tops, "results/250221/GBM热图表格.csv")
ggsave("results/250221/GBM热图.pdf", p, width = 14, height = 10)

# %%
compare.gene <- function(seurat.obj, gene) {
    draw.data <- GetAssayData(seurat.obj, slot = "data", assay = "integrated")
    raw.df <- GetAssayData(seurat.obj, slot = "counts", assay = "Spatial")
    express <- draw.data[gene, raw.df[gene, ] > 0]
    draw.df <- data.frame(
        "expression" = express,
        "region" = Idents(seurat.obj)[names(express)],
        row.names = names(express)
    )
    mean.df <- draw.df %>%
        group_by(region) %>%
        summarize(across(everything(), .fns = mean, na.rm = TRUE)) %>%
        as.data.frame()
    if (mean.df[1, 2] > mean.df[2, 2]) {
        if (mean.df[2, 2] > mean.df[3, 2]) {
            return(gene)
        }
    }
}
up.genes <- mclapply(
    unique(markers.iv$gene),
    function(x) compare.gene(integrate.sub, x),
    mc.cores = 3
)
up.genes <- unlist(up.genes)
enrich.list <- enrichmentGenelist(
    up.genes, fs::path(WORKDIR, "results", "250218", "小提琴图基因")
)
p <- dotplot(enrich.list$ego, showCategory = draw.pathways$Description)
ggsave(fs::path(WORKDIR, "results", "250220", "小提琴图基因.pdf"), p)

# %%
tops <- markers.iv %>%
    filter(!grepl("^MT-", gene)) %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 30) %>%
    ungroup()
p <- DoHeatmap(integrate.sub, features = tops$gene)
ggsave("热图.pdf", p, width = 10, height = 10)

# %%
ave.exp <- AverageExpression(integrate.sub)$integrated
select.genes <- ave.exp[ave.exp[, 1] > ave.exp[, 2], ]
select.genes <- as.data.frame(
    select.genes[select.genes[, 2] > select.genes[, 3], ]
)

select.genes <- filter(select.genes, rownames(select.genes) %in% markers.iv$gene)

# %%
violin.plot <- function(seurat.obj, gene) {
    draw.data <- GetAssayData(seurat.obj, slot = "data", assay = "integrated")
    raw.df <- GetAssayData(seurat.obj, slot = "counts", assay = "Spatial")
    express <- draw.data[gene, raw.df[gene, ] > 0]
    draw.df <- data.frame(
        "expression" = express,
        "region" = Idents(seurat.obj)[names(express)],
        row.names = names(express)
    )
    mean.df <- draw.df %>%
        group_by(region) %>%
        summarize(across(everything(), .fns = mean, na.rm = TRUE)) %>%
        as.data.frame()
    if (mean.df[1, 2] > mean.df[2, 2]) {
        if (mean.df[2, 2] > mean.df[3, 2]) {
            p <- ggplot(mapping = aes(x = region, y = expression)) +
                geom_violin(data = draw.df) +
                geom_line(data = mean.df, group = 1) +
                geom_point(data = mean.df) +
                labs(title = gene) +
                theme_classic() +
                theme(plot.title = element_text(face = "bold.italic", hjust = 0.5))
            ggsave(
                fs::path(WORKDIR, "gene-expression", "violin", paste0(gene, ".pdf")),
                p
            )
        }
    }
}
mclapply(
    unique(markers.iv$gene),
    function(x) violin.plot(integrate.sub, x),
    mc.cores = 5
)

# %%
markers <- FindAllMarkers(integrate.obj, min.pct = 0.3)
write.csv(markers, "markers.csv")

# %%
gene.list <- list(
    "ubiquitin" = loadUbiquitin(),
    "rbp" = loadRBP(),
    "kinase" = loadKinase(),
    "phosphatase" = loadDEPOD()
)
max_length <- max(sapply(gene.list, length))
gene.list.filled <- lapply(gene.list, function(x) {
      c(x, rep(NA, max_length - length(x)))
})
gene.df <- as.data.frame(gene.list.filled)
write.csv(gene.df, "gene_list.csv", row.names = FALSE, na = "", quote = FALSE)

draw.obj.all <- integrate.obj
region.annotations <- list(
    "IV Tumor cell densely populated area" = "IDH mutant Tumor cell densely\npopulated area",
    "GBM Tumor cell densely populated area" = "IDH wildtype Tumor cell densely\npopulated area",
    "IV Tumor area 1" = "IDH mutant Tumor area 1",
    "IV Tumor area 2" = "IDH mutant Tumor area 2",
    "IV Tumor area 3" = "IDH mutant Tumor area 3",
    "GBM Tumor area 1" = "IDH wildtype Tumor area 1",
    "GBM Tumor area 2" = "IDH wildtype Tumor area 2",
    "GBM Tumor area 3" = "IDH wildtype Tumor area 3",
    "GBM Tumor area 4" = "IDH wildtype Tumor area 4",
    "GBM Tumor area 5" = "IDH wildtype Tumor area 5",
    "GBM Tumor area 6" = "IDH wildtype Tumor area 6",
    "GBM Tumor area 7" = "IDH wildtype Tumor area 7",
    "GBM Tumor area 8" = "IDH wildtype Tumor area 8",
    "Normal tissue adjacent to tumor area" = "Normal tissue adjacent\nto tumor area"
)
draw.obj.all <- RenameIdents(draw.obj.all, region.annotations)
markers$cluster_1 <- ifelse(
    markers$cluster %in% names(region.annotations),
    region.annotations[markers$cluster],
    markers$cluster
)

# %%
draw.gene.list <- function(genes, level, draw.levels = NULL) {
    if (level == "IV") {
        height <- 14
        width <- 8
    } else {
        height <- 16
        width <- 12
    }
    save.path <- fs::path(
        WORKDIR,
        "results",
        "250307",
        paste(level, genes, "pdf", sep = ".")
    )
    draw.obj <- draw.obj.all[, draw.obj.all$level == level]
    if (is.null(draw.levels)) {
        draw.levels <- levels(draw.obj)
    } else {
        draw.levels <- c(
            draw.levels,
            levels(draw.obj)[!levels(draw.obj) %in% draw.levels]
        )
    }
    Idents(draw.obj) <- factor(Idents(draw.obj), levels = draw.levels)
    genes <- markers %>%
        mutate(cluster_1 = factor(cluster_1, levels = draw.levels)) %>%
        filter(p_val_adj < 0.05, gene %in% gene.list[[genes]]) %>%
        group_by(cluster_1) %>%
        arrange(desc(avg_log2FC)) %>%
        slice_head(n = 10)
    p <- DoHeatmap(draw.obj, genes$gene)
    ggsave(save.path, p, width = width, height = height)
}

draw.levels <- c(
    "Normal tissue adjacent\nto tumor area",
    "Blood vessel rich area"
)
draw.gene.list(genes = "kinase", level = "IV", draw.levels = draw.levels)

draw.levels <- c(
    "IDH wildtype Tumor area 2",
    "IDH wildtype Tumor area 3",
    "IDH wildtype Tumor area 4",
    "IDH wildtype Tumor area 6",
    "IDH wildtype Tumor area 7",
    "IDH wildtype Tumor area 8",
    "Blood vessel rich area"
)
draw.gene.list(genes = "kinase", level = "GBM", draw.levels = draw.levels)

draw.levels <- c(
    "Normal tissue adjacent\nto tumor area",
    "Blood vessel rich area"
)
draw.gene.list(genes = "phosphatase", level = "IV", draw.levels = draw.levels)

draw.levels <- c(
    "IDH wildtype Tumor area 2",
    "IDH wildtype Tumor area 3",
    "IDH wildtype Tumor area 4",
    "IDH wildtype Tumor area 6",
    "IDH wildtype Tumor area 7",
    "IDH wildtype Tumor area 8",
    "Blood vessel rich area"
)
draw.gene.list(genes = "phosphatase", level = "GBM", draw.levels = draw.levels)

draw.levels <- c(
    "Normal tissue adjacent\nto tumor area",
    "Blood vessel rich area",
    "IDH mutant Tumor cell densely\npopulated area"
)
draw.gene.list(genes = "rbp", level = "IV", draw.levels = draw.levels)

draw.levels <- c(
    "IDH wildtype Tumor area 2",
    "IDH wildtype Tumor area 3",
    "IDH wildtype Tumor area 6",
    "IDH wildtype Tumor area 7",
    "IDH wildtype Tumor area 8",
    "Blood vessel rich area"
)
draw.gene.list(genes = "rbp", level = "GBM", draw.levels = draw.levels)

draw.levels <- c(
    "Normal tissue adjacent\nto tumor area",
    "Blood vessel rich area",
    "IDH mutant Tumor cell densely\npopulated area"
)
draw.gene.list(genes = "ubiquitin", level = "IV", draw.levels = draw.levels)

draw.levels <- c(
    "IDH wildtype Tumor area 6",
    "IDH wildtype Tumor area 7",
    "IDH wildtype Tumor area 8",
    "Blood vessel rich area"
)
draw.gene.list(genes = "ubiquitin", level = "GBM", draw.levels = draw.levels)

# %%
gene <- "MSTN"
p <- FeaturePlot(
    integrate.obj,
    features = gene,
    reduction = "tsne",
    min.cutoff = 0,
    max.cutoff = "q90",
    ) + NoLegend()
ggsave(paste0(gene, "-tsne.pdf"), p)

draw.obj <- integrate.obj
DefaultAssay(draw.obj) <- "Spatial"
p <- SpatialFeaturePlot(
    draw.obj,
    features = gene,
    slot = "count",
    min.cutoff = 0,
    max.cutoff = "q95"
)
ggsave(paste0(gene, "-spatial.pdf"), p, width = 28)

# %%
p <- DotPlot(
    integrate.obj,
    features = c("kinase3", "phosphatase4", "ubiquitin1", "rbp2"),
    scale = FALSE
)
p <- ggplot(p$data, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
    geom_tile() +
    coord_flip() +
    theme(panel.grid = element_blank(), text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gsea()
ggsave(fs::path(WORKDIR, "results", "module-score.hm.pdf"), p, width = 9, height = 5)
