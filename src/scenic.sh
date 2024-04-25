#!/usr/bin/env bash

WORKDIR=$HOME/workspace/gbm
cd $WORKDIR

source venv/bin/activate

ranking_10k=Data/scenic/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
ranking_500bp=Data/scenic/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
ann_fname=Data/scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

pyscenic grn \
   Data/counts/$1.csv \
   Data/scenic/allTFs_hg38.txt \
   -o results/scenic/$1.adj.csv \
   --num_workers 10

pyscenic ctx \
    results/scenic/$1.adj.csv \
    $ranking_10k $ranking_500bp \
    --annotations_fname $ann_fname \
    --expression_mtx_fname Data/counts/$1.csv \
    --output results/scenic/$1.motifs.csv \
    --num_workers 10

pyscenic aucell \
    Data/counts/$1.csv \
    results/scenic/$1.motifs.csv \
    --output results/scenic/$1.auc.csv \
    --num_workers 10

python src/scenic-downstream.py $1
