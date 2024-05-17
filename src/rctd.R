#! /usr/bin/env Rscript

# %%
library <- function(...) suppressMessages(base::library(...))
library('parallel')
library('dplyr')
library('Seurat')
library('spacexr')
library('ggplot2')

WORKDIR <- fs::path(Sys.getenv("HOME"), "workspace", "gbm")
source(fs::path(WORKDIR, "src", "R", "utils.R"))

# %% read st data
loadSeuratList()
createPuck <- function(idx) {
    coor.df <- GetTissueCoordinates(seurat.list[[idx]])
    colnames(coor.df) <- c("y", "x")
    coor.df <- coor.df[, c("x", "y")]
    st.df <- GetAssayData(
        seurat.list[[idx]], assay = "Spatial", slot = "count"
    )
    numi.df <- seurat.list$nCount_Spatial
    puck <- SpatialRNA(coor.df, st.df, numi.df)
    return(puck)
}
puck.list <- mclapply(idx.list, createPuck)

# %% read sc data
sc.path <- fs::path(WORKDIR, "Data", "singlecell", "SCP503", "other")
meta.df <- read.csv(
    fs::path(sc.path, "Richards_NatureCancer_GBM_scRNAseq_meta.csv"),
    row.names = 1,
    header = TRUE
)
meta.df <- meta.df[-1, ]
rownames(meta.df) <- gsub("-", ".", rownames(meta.df))

numi.sc <- meta.df$nUMI
names(numi.sc) <- rownames(meta.df)
cluster.sc <- as.factor(meta.df$CellType)
names(cluster.sc) <- rownames(meta.df)

sc.df <- read.table(
    fs::path(sc.path, "Richards_NatureCancer_GBM_scRNAseq_counts.csv"),
    row.names = 1,
    sep = ",",
    check.names = FALSE,
    header = TRUE
)
sc.df <- sc.df[, rownames(meta.df)]

ref.sc <- Reference(sc.df, cluster.sc, numi.sc)

# %%
runRCTD <- function(idx) {
    rctd.obj <- create.RCTD(puck.list[[idx]], ref.sc, max_cores = 1)
    rctd.obj <- run.RCTD(rctd.obj, doublet_mode = 'doublet')
    return(rctd.obj)
}
rctd.results <- mclapply(idx.list, runRCTD)

# %%
for (idx in names(rctd.results)) {
    results.df <- rctd.results[[idx]]@results$results_df
    weights.df <- rctd.results[[idx]]@results$weights %>% as.data.frame()
    results.df$Immune <- weights.df$Immune
    results.df$NormalBrain <- weights.df$NormalBrain
    results.df$Tumour <- weights.df$Tumour
    write.csv(
        results.df,
        fs::path(save.dirs[[idx]], paste0(idx, ".rctd.csv"))
    )
}
