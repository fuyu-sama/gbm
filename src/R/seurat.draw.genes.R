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
library('ggplot2')
library('Seurat')

# options
options(mc.cores = 5)
options(future.globals.maxSize = 500 * 1024 ^ 3)

# public variables
ggsave <- function(...) suppressMessages(ggplot2::ggsave(...))
WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
idx.full <- c("21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2")
args <- commandArgs(trailingOnly = TRUE)

# public functions
loadSeuratList <- function() {
    seurat.list <<- readRDS(fs::path(WORKDIR, "results", "seurat.list.rds"))
}

loadIntegrate <- function() {
    integrate.obj <<- readRDS(fs::path(WORKDIR, "results", "integrate.rds"))
}

# %% draw functions
draw.raw <- function(gene) {
    if (grepl("^DEPRECATED", gene)) return()
    DefaultAssay(seurat.obj) <- "Spatial"
    save.dir <- fs::path(WORKDIR, "gene-expression", "raw")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    p <- SpatialFeaturePlot(seurat.obj, gene, slot = "count",  max.cutoff = "q95")
    ggsave(fs::path(save.dir, paste0(gene, ".pdf")), p, width = 28)
}

draw.int <- function(gene) {
    if (grepl("^DEPRECATED", gene)) return()
    DefaultAssay(integrate.obj) <- "SCT"
    save.dir <- fs::path(WORKDIR, "gene-expression", "int")
    if (! fs::dir_exists(save.dir)) fs::dir_create(save.dir)
    p <- SpatialFeaturePlot(integrate.obj, gene, max.cutoff = "q95")
    ggsave(fs::path(save.dir, paste0(gene, ".pdf")), p, width = 28)
}

# %% run
loadIntegrate()
if (args[1] == "raw") {
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
    mclapply(rownames(seurat.obj), draw.raw)
} else if (args[1] == "int") {
    mclapply(rownames(integrate.obj), draw.int)
}

