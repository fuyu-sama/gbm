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
library('parallel')
library('xlsx')
library('ggplot2')
library('Seurat')

options(mc.cores = 4)

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

# %% read data
loadSeuratList()
regionAnnotation <- function(seurat.obj) {
    idx <- names(seurat.obj@images)
    region.annotations <- list()
    region.annotations[["21B-603-5"]] <- c(
        "4" = "Blood vessel rich area",
        "2" = "Tumor cell densely populated area",
        "0" = "Tumor area 1",
        "1" = "Junction area",
        "3" = "Parancerous area"
    )
    region.annotations[["22F-10823-3"]] <- c(
        "2" = "Blood vessel rich area",
        "0" = "Tumor cell densely populated area",
        "5" = "Tumor area 2",
        "3" = "Junction area 2",
        "1" = "Junction area",
        "4" = "Parancerous area"
    )
    region.annotations[["22F-21576-1"]] <- c(
        "0" = "Blood vessel rich area",
        "2" = "Tumor cell densely populated area",
        "3" = "Tumor area 3",
        "1" = "Tumor area 4"
    )
    region.annotations[["22F-23738-2"]] <- c(
        "2" = "Blood vessel rich area",
        "0" = "Tumor cell densely populated area",
        "1" = "Tumor area 5",
        "4" = "Tumor area 6",
        "3" = "Tumor area 7",
        "5" = "Tumor area 8"
    )
    seurat.obj <- RenameIdents(seurat.obj, region.annotations[[idx]])
    return(seurat.obj)
}
seurat.list <- mclapply(seurat.list, regionAnnotation)

# %%
gene.df <- read.xlsx2(
    fs::path(WORKDIR, "Data", "Death", "铁死亡及自噬基因集.xlsx"), 1)

# %%
for (idx in idx.full) {
    seurat.obj <- seurat.list[[idx]]
    DefaultAssay(seurat.obj) <- "SCT"

    features <- c()
    for (feature in gene.df[["Autophagy.657."]]) {
        if (feature %in% rownames(seurat.obj)) features <- c(features, feature)
    }
    p <- DoHeatmap(seurat.obj, features = features, slot = "scale.data") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    ggsave(
        fs::path(WORKDIR, "gene-expression", paste0(idx, "-autophagy.pdf")), p)

    features <- c()
    for (feature in gene.df[["Ferroptosis.564."]]) {
        if (feature %in% rownames(seurat.obj)) features <- c(features, feature)
    }
    p <- DoHeatmap(seurat.obj, features = features, slot = "scale.data") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    ggsave(
        fs::path(WORKDIR, "gene-expression", paste0(idx, "-ferroptosis.pdf")), p)
}
