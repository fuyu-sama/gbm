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
library('dplyr')
library('ggplot2')
library('clusterProfiler')
library('org.Hs.eg.db')

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
loadSeuratList <- function() {
    seurat.list <<- readRDS(fs::path(WORKDIR, "results", "seurat.list.rds"))
}

loadIntegrate <- function() {
    integrate.obj <<- readRDS(fs::path(WORKDIR, "results", "integrate.rds"))
}

# %% enrichment functions
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
        qvalueCutoff = 0.1
    )

    write.csv(as.data.frame(ego), fs::path(save.dir.go, "GO.csv"))

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
        qvalueCutoff = 0.1
    )
    write.csv(as.data.frame(kk), fs::path(save.dir.kegg, "KEGG.csv"))
}

# %% DE 1
markers.list1 <- list()
for (idx in names(seurat.list)) {
    save.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区域差异表达.csv"))
    markers.list1[[idx]] <- read.csv(save.path, header = TRUE, row.names = 1)
}

# %%
idx <- idx.full[1]
save.dir <- fs::path(save.dirs[[idx]], "GO")
cluster <- "Tumor_cell_densely_populated_area"
pathway.list <- list()
read.df <- xlsx::read.xlsx2(
    fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
    row.names = 1,
    sheetName = "avg_log2FC > 0"
)
pathway.list[[idx.full[1]]] <- read.df$Description

idx <- idx.full[2]
save.dir <- fs::path(save.dirs[[idx]], "GO")
cluster <- "Tumor_cell_densely_populated_area"
read.df <- xlsx::read.xlsx2(
    fs::path(save.dir, paste(idx, cluster, "GO.xlsx", sep = ".")),
    row.names = 1,
    sheetName = "avg_log2FC > 0"
)
pathway.list[[idx.full[2]]] <- read.df$Description

overlaps <- gplots::venn(pathway.list, show.plot = FALSE)
overlaps <- attributes(overlaps)$intersections
capture.output(
    overlaps,
    file = fs::path(save.path, "GO.txt")
)

VennDiagram::venn.diagram(
    pathway.list,
    fs::path(save.path, "GO-venn.pdf"),
    imagetype = "pdf",
    disable.logging = TRUE
)

# %%
save.path <- fs::path(WORKDIR, "results", "密集区交集")
if (! fs::dir_exists(save.path)) fs::dir_create(save.path)

densly.genes <- list()
genes.list1 <- markers.list1[[idx.full[1]]] %>%
    filter(cluster == "Tumor cell densely populated area") %>%
    filter(p_val_adj < 0.05, avg_log2FC > .5)
densly.genes[[idx.full[1]]] <- genes.list1$gene
genes.list2 <- markers.list1[[idx.full[2]]] %>%
    filter(cluster == "Tumor cell densely populated area") %>%
    filter(p_val_adj < 0.05, avg_log2FC > .5)
densly.genes[[idx.full[2]]] <- genes.list2$gene

overlaps <- gplots::venn(densly.genes, show.plot = FALSE)
overlaps <- attributes(overlaps)$intersections

VennDiagram::venn.diagram(
    densly.genes,
    fs::path(save.path, "基因venn.png"),
    imagetype = "png",
    disable.logging = TRUE
)

#enrichmentGenelist(overlaps[["21B-603-5:22F-10823-3"]], save.path)
