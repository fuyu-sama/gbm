library('forcats')
library('org.Hs.eg.db')
# savedir
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
colors <- c(
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F",
    "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4",
    "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
    "#CAB2D6"
)
cols <- c(
    "Blood vessel rich area" = "#666666",
    "IV Tumor cell densely populated area" = "#FF7F00",
    "GBM Tumor cell densely populated area" = "#E31A1C",
    "IV Tumor area 1" = "#FB9A99",
    "IV Tumor area 2" = "#386CB0",
    "IV Tumor area 3" = "#E7298A",
    "GBM Tumor area 1" = "#D95F02",
    "GBM Tumor area 2" = "#FDBF6F",
    "GBM Tumor area 3" = "#FFFF99",
    "GBM Tumor area 4" = "#F0027F",
    "GBM Tumor area 5" = "#BF5B17",
    "GBM Tumor area 6" = "#E6AB02",
    "GBM Tumor area 7" = "#FB9A99",
    "GBM Tumor area 8" = "#A6CEE3",
    "Junction area" = "#7FC97F",
    "Normal tissue adjacent to tumor area" = "#33A02C"
)

# options
future::plan("multicore", workers = 4)
options(mc.cores = 4)
options(future.globals.maxSize = 500 * 1024 ^ 3)

OrgDb <- org.Hs.eg.db

ggsave <- function(...) suppressMessages(ggplot2::ggsave(...))

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

loadKinase <- function() {
    kinase.df <- jsonlite::fromJSON(fs::path(WORKDIR, "Data", "klifs.net.json"))
    return(kinase.df$name)
}

loadUbiquitin <- function() {
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
    return(ubiquitin.df$V6)
}

loadTF <- function() {
    tf.df <- read.csv(fs::path(WORKDIR, "Data", "DatabaseExtract_v_1.01.csv"))
    tf.df <- filter(tf.df, TF.assessment == "Known motif")
    return(tf.df$HGNC.symbol)
}

loadRBP <- function() {
    rbp.df <- read.csv(fs::path(WORKDIR, "Data", "cancers-751631-suppl-final.csv"))
    return(rbp.df$ID)
}

drawGenelist <- function(seurat.obj, gene.list, save.name) {
    idx <- names(seurat.obj@images)

    marker.genes <- markers.list1[[idx]] %>%
        group_by(cluster) %>%
        filter(p_val_adj < 0.05, avg_log2FC > .5) %>%
        arrange(desc(avg_log2FC), .by_group = TRUE) %>%
        filter(gene %in% gene.list)

    write.csv(
        marker.genes,
        fs::path(save.dirs[[idx]], paste(idx, save.name, "csv", sep = "."))
    )

    p <- DoHeatmap(seurat.obj, features = marker.genes$gene)
    save.path <- fs::path(
        save.dirs[[idx]], paste(idx, save.name, "pdf", sep = ".")
    )
    ggsave(save.path, p, height = 20, width = 15)

}

# enrichment functions
kkMapIds <- function(kk) {
    ids <- strsplit(kk$geneID, "/")
    gene.names <- c()
    for (line in ids) {
        genes <- mapIds(
            OrgDb, keys = line, column = "SYMBOL", keytype = "ENTREZID"
        )
        genes <- paste(genes, collapse = "/")
        gene.names <- c(gene.names, genes)
    }
    kk$geneName <- gene.names
    return(kk)
}

enrichmentFindAllMarkers <- function(
    cluster, markers, save.dir, logfc.threshold = 0.5) {
    upscale.table <- markers %>%
        filter(p_val_adj < 0.05, avg_log2FC > logfc.threshold)
    downscale.table <- markers %>%
        filter(p_val_adj < 0.05, avg_log2FC < -logfc.threshold)
    upscale.genes <- upscale.table[upscale.table$cluster == cluster, ]$gene
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- downscale.table[downscale.table$cluster == cluster, ]$gene
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]

    save.dir.go <- fs::path(save.dir, "GO")
    if (! fs::dir_exists(save.dir.go)) fs::dir_create(save.dir.go)
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
    p1 <- dotplot(upscale.ego) +
        ggtitle("GO Enrichment for upscale genes")
    p2 <- dotplot(downscale.ego, split = "ONTOLOGY") +
        ggtitle("GO Enrichment for downscale genes")
    cluster <- stringr::str_replace_all(cluster, " ", "_")
    ggsave(
        fs::path(save.dir.go, paste(cluster, "GO.pdf", sep = ".")),
        p1 + p2,
        width = 14
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir.go, paste(cluster, "GO.xlsx", sep = ".")),
        sheetName = "upscale"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir.go, paste(cluster, "GO.xlsx", sep = ".")),
        sheetName = "downscale",
        append = TRUE
    )

    save.dir.kegg <- fs::path(save.dir, "KEGG")
    if (! fs::dir_exists(save.dir.kegg)) fs::dir_create(save.dir.kegg)
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
            ggtitle("KEGG Enrichment for upscale genes")
        write.xlsx2(
            kkMapIds(as.data.frame(upscale.kk)),
            fs::path(save.dir.kegg, paste(cluster, "KEGG.xlsx", sep = ".")),
            sheetName = "upscale"
        )
    } else {
        p1 <- NULL
    }
    if (dim(downscale.kk)[1] > 0) {
        p2 <- dotplot(downscale.kk) +
            ggtitle("KEGG Enrichment for downscale genes")
        write.xlsx2(
            kkMapIds(as.data.frame(downscale.kk)),
            fs::path(save.dir.kegg, paste(cluster, "KEGG.xlsx", sep = ".")),
            sheetName = "downscale",
            append = TRUE
        )
    } else {
        p2 <- NULL
    }
    ggsave(
        fs::path(save.dir.kegg, paste(cluster, "KEGG.pdf", sep = ".")),
        p1 + p2,
        width = 14
    )
}

enrichmentFindMarkers <- function(markers, save.dir, logfc.threshold = 0.5) {
    upscale.table <- markers %>%
        filter(p_val_adj < 0.05, avg_log2FC > logfc.threshold)
    downscale.table <- markers %>%
        filter(p_val_adj < 0.05, avg_log2FC < -logfc.threshold)
    upscale.genes <- rownames(upscale.table)
    upscale.genes <- upscale.genes[!grepl("^DEPRECATED-", upscale.genes)]
    downscale.genes <- rownames(downscale.table)
    downscale.genes <- downscale.genes[!grepl("^DEPRECATED-", downscale.genes)]

    save.dir.go <- fs::path(save.dir, "GO")
    if (! fs::dir_exists(save.dir.go)) fs::dir_create(save.dir.go)

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
    p1 <- dotplot(upscale.ego) +
        ggtitle("GO Enrichment for upscale genes")
    p2 <- dotplot(downscale.ego) +
        ggtitle("GO Enrichment for downscale genes")
    ggsave(
        fs::path(save.dir.go, "GO.pdf"),
        p1 + p2,
        width = 14
    )
    write.xlsx2(
        as.data.frame(upscale.ego),
        fs::path(save.dir.go, "GO.xlsx"),
        sheetName = "upscale"
    )
    write.xlsx2(
        as.data.frame(downscale.ego),
        fs::path(save.dir.go, "GO.xlsx"),
        sheetName = "downscale",
        append = TRUE
    )

    save.dir.kegg <- fs::path(save.dir, "KEGG")
    if (! fs::dir_exists(save.dir.kegg)) fs::dir_create(save.dir.kegg)
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
            ggtitle("KEGG Enrichment for upscale genes")
        write.xlsx2(
            kkMapIds(as.data.frame(upscale.kk)),
            fs::path(save.dir.kegg, "KEGG.xlsx"),
            sheetName = "upscale"
        )
    } else {
        p1 <- NULL
    }
    if (dim(downscale.kk)[1] > 0) {
        p2 <- dotplot(downscale.kk) +
            ggtitle("KEGG Enrichment for downscale genes")
        write.xlsx2(
            kkMapIds(as.data.frame(downscale.kk)),
            fs::path(save.dir.kegg, "KEGG.xlsx"),
            sheetName = "downscale",
            append = TRUE
        )
    } else {
        p2 <- NULL
    }
    ggsave(fs::path(save.dir.kegg, "KEGG.pdf"), p1 + p2, width = 14)
}

enrichmentGenelist <- function(markers, save.dir) {
    markers <- markers[!grepl("^DEPRECATED-", markers)]

    save.dir.go <- fs::path(save.dir, "GO")
    if (! fs::dir_exists(save.dir.go)) fs::dir_create(save.dir.go)

    ego <- enrichGO(
        markers,
        OrgDb = OrgDb,
        keyType = "SYMBOL",
        ont = "ALL",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
    )

    write.csv(as.data.frame(ego), fs::path(save.dir.go, "GO.csv"))
    if (dim(ego)[1] > 1) {
        if (dim(ego)[1] < 20) {
            showCategory <- dim(ego)[1]
        } else {
            showCategory <- 20
        }
        p <- dotplot(ego, showCategory = showCategory) + ggtitle("GO Enrichment")
        ggsave(fs::path(save.dir.go, "GO.pdf"), p, height = 10)
    }

    save.dir.kegg <- fs::path(save.dir, "KEGG")
    if (! fs::dir_exists(save.dir.kegg)) fs::dir_create(save.dir.kegg)
    ids <- mapIds(
        OrgDb, keys = markers, column = "ENTREZID", keytype = "SYMBOL")
    kk <- enrichKEGG(
        ids,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2
    )
    write.csv(kkMapIds(as.data.frame(kk)), fs::path(save.dir.kegg, "KEGG.csv"))
    if (dim(kk)[1] > 1) {
        if (dim(kk)[1] < 20) {
            showCategory <- dim(kk)[1]
        } else {
            showCategory <- 20
        }
        p <- dotplot(kk, showCategory = showCategory) + ggtitle("KEGG Enrichment")
        ggsave(fs::path(save.dir.kegg, "KEGG.pdf"), p, height = 10)
    }
}

DoComplexHeatmap <- function(
    object,
    save.path,
    features = NULL,
    cells = NULL,
    group.by = "ident",
    slot = "scale.data",
    assay = NULL
    ) {
    if (is.null(features)) features <- rownames(object)
    if (is.null(cells)) cells <- colnames(object)
    if (!is.null(assay)) DefaultAssay(object) <- assay
    names(object@images) <- NULL
    object <- object[features, cells]
    if (group.by == "ident") {
        o <- sort(Idents(object))
        object <- object[features, names(o)]
    }

    expression.df <- GetAssayData(object, slot = slot)[features, names(o)]
    ha <- HeatmapAnnotation(region = as.vector(o))
    ht <- Heatmap(
        expression.df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        top_annotation = ha
    )

    pdf(save.path)
    draw(ht)
    dev.off()
}

ggvolcano <- function(
    dataset, cut_off_pvalue = 0.01, cut_off_logFC = 0.5, genes = NULL) {

    dataset$change <- ifelse(
        dataset$p_val_adj < cut_off_pvalue & abs(dataset$avg_log2FC) >= cut_off_logFC,
        ifelse(dataset$avg_log2FC > cut_off_logFC, 'Up', 'Down'),
        'Stable'
    )

    if (!is.null(genes)) {
        dataset$label <- ifelse(dataset$Symbol %in% genes, dataset$Symbol, "")
    }

    p <- ggplot(dataset, aes(x = avg_log2FC, y = -log10(p_val_adj), color = change)) +
        geom_point(alpha = 0.4, size = .5) +
        theme_bw(base_size = 12) +
        xlab("Averange log2FC") +
        ylab("-Log10(Adjusted p value)") +
        theme(plot.title = element_text(size = 15, hjust = 0.5)) +
        scale_colour_manual(values = c('steelblue', 'gray', 'brown')) +
        geom_hline(yintercept = -log10(0.05), lty = 4) +
        geom_vline(xintercept = c(-cut_off_logFC, cut_off_logFC), lty = 4)

    if (!is.null(genes)) {
        p <- p + geom_label_repel(
            data = dataset,
            aes(label = label),
            size = 3,
            box.padding = unit(0.5, "lines"),
            point.padding = unit(0.8, "lines"),
            segment.color = "black",
            show.legend = FALSE,
            max.overlaps = 50000)
    }

    return(p)
}

pDotPlot <- function(...) {
    p <- DotPlot(...)
    exp <- p$data
    exp$features.plot <- as.factor(exp$features.plot)
    exp$features.plot <- fct_inorder(exp$features.plot)
    p <- ggplot(exp, aes(x = id, y = features.plot)) +
        geom_point(aes(size = `pct.exp`, color = `avg.exp.scaled`)) +
        theme_classic() +
        coord_flip() +
        theme(panel.grid = element_blank(), text = element_text(size = 20)) +
        scale_color_gradient2(low = "blue", mid = "lightgrey", high =  "red") +
        labs(x = NULL, y = NULL) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(p)
}
