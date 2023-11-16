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
library('xlsx')
library('parallel')
library('dplyr')
library('ggplot2')
library('Seurat')
library('org.Hs.eg.db')

# options
options(mc.cores = 5)
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
loadSeuratList <- function() {
    seurat.list <<- readRDS(fs::path(WORKDIR, "results", "seurat.list.rds"))
}

loadIntegrate <- function() {
    integrate.obj <<- readRDS(fs::path(WORKDIR, "results", "integrate.rds"))
}

# %% load seurat data
loadIntegrate()
loadSeuratList()
for (sample in names(seurat.list)) {
    names(seurat.list[[sample]]@images) <- NULL
    seurat.list[[sample]] <- RenameCells(
        seurat.list[[sample]], add.cell.id = sample)
}
seurat.obj <- merge(seurat.list[[idx.full[1]]], seurat.list[[idx.full[2]]])
seurat.obj <- merge(seurat.obj, seurat.list[[idx.full[3]]])
seurat.obj <- merge(seurat.obj, seurat.list[[idx.full[4]]])
seurat.obj@images <- integrate.obj@images

# %%
path2desc <- list(
    "GO:0042063" = "gliogenesis",
    "GO:0010001" = "glial cell differentiation",
    "GO:0021782" = "glial cell development",
    "GO:0014003" = "oligodendrocyte development",
    "GO:0048709" = "oligodendrocyte differentiation",
    "GO:0007015" = "actin filament organization"
)

path2ext <- clusterProfiler:::get_GO_data(
    "org.Hs.eg.db", "ALL", "SYMBOL")$PATHID2EXTID

de.list <- list()
for (idx in idx.full) {
    de.path <- fs::path(save.dirs[[idx]], paste0(idx, ".分区域差异表达.csv"))
    de.list[[idx]] <- read.csv(de.path, row.names = 1)
}

# %%
part.dir <- fs::path(WORKDIR, "results", "GO.part")
if (! fs::dir_exists(part.dir)) fs::dir_create(part.dir)
draw <- function(path) {
    genes <- path2ext[[path]]
    save.dir <- fs::path(part.dir, path2desc[[path]])
    gene.dir <- fs::path(save.dir, "gene-expression")
    if (! fs::dir_exists(save.dir)) {
        fs::dir_create(save.dir)
        fs::dir_create(gene.dir)
    }
    for (idx in idx.full) {
        de.sub <- filter(de.list[[idx]], gene %in% genes)
        write.csv(de.sub, fs::path(save.dir, paste0(idx, ".分区域差异表达.csv")))
    }
    for (gene in genes) {
        if (! gene %in% rownames(seurat.obj)) next
        p <- SpatialFeaturePlot(seurat.obj, gene, max.cutoff = "q95")
        ggsave(fs::path(gene.dir, paste0(gene, ".pdf")), p)
    }
}
mclapply(names(path2desc), draw)
