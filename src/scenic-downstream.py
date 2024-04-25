#! /usr/bin/env python3
# -*- encoding: utf-8 -*-

# %%
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")

idx = sys.argv[1]

# %% read data
count_df = pd.read_csv(
    Path.joinpath(WORKDIR, "Data", "counts", f"{idx}.csv"),
    index_col=0,
    header=0,
)
auc_df = pd.read_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.auc.csv"),
    index_col=0,
    header=0,
)

adj_df = pd.read_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.adj.csv"),
    index_col=False,
)

idents = pd.read_csv(
    Path.joinpath(WORKDIR, "results", "idents.csv"),
    index_col=0,
).loc[auc_df.index, :]

# %% binarize auc
bin_df, thresholds = binarize(auc_df)
bin_df.to_csv(Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.bin.csv"))
thresholds.to_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.threshold.csv"),
    header=False,
)

# %% rss for GBM vs. IV
idents_sub = idents[idents["x"] != "Normal tissue adjacent to tumor area"]
idents_sub = idents_sub[idents_sub["x"] != "Junction area"]
idents_sub["idx"] = [i.split("_")[0] for i in idents_sub.index]
idents_sub["level"] = np.where(
    idents_sub["idx"].isin(["21B-603-5", "22F-10823-3"]),
    "IV",
    "GBM",
)
rss_level = regulon_specificity_scores(
    auc_df.loc[idents_sub.index, :],
    idents_sub["level"],
)
rss_level.to_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.rss_level.csv"))

cats = sorted(list(set(idents_sub["level"])))
fig = plt.figure(figsize=(6, 4))
for c, num in zip(cats, range(1, len(cats) + 1)):
    x = rss_level.T[c]
    ax = fig.add_subplot(1, 2, num)
    plot_rss(rss_level, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05,
                x.max() + (x.max() - x.min()) * 0.05)
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(
        ax.texts,
        autoalign='xy',
        ha='right',
        va='bottom',
        arrowprops=dict(arrowstyle='-', color='lightgrey'),
        precision=0.001,
    )
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(
    0.00,
    0.5,
    'Regulon specificity score (RSS)',
    ha='center',
    va='center',
    rotation='vertical',
    size='x-large',
)
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})
plt.savefig(
    Path.joinpath(WORKDIR, "results", "scenic-plot", f"{idx}.rss.level.pdf"),
    dpi=150,
    bbox_inches="tight",
)

# %% calculate rss for all region
rss_celltype = regulon_specificity_scores(auc_df, idents["x"])
rss_celltype.to_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.rss_celltype.csv"))

cats = sorted(list(set(idents["x"])))
fig = plt.figure(figsize=(12, 16))
for c, num in zip(cats, range(1, len(cats) + 1)):
    x = rss_celltype.T[c]
    ax = fig.add_subplot(4, 4, num)
    plot_rss(rss_celltype, c, top_n=5, max_n=None, ax=ax)
    ax.set_ylim(x.min() - (x.max() - x.min()) * 0.05,
                x.max() + (x.max() - x.min()) * 0.05)
    for t in ax.texts:
        t.set_fontsize(12)
    ax.set_ylabel('')
    ax.set_xlabel('')
    adjust_text(
        ax.texts,
        autoalign='xy',
        ha='right',
        va='bottom',
        arrowprops=dict(arrowstyle='-', color='lightgrey'),
        precision=0.001,
    )
fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
fig.text(
    0.00,
    0.5,
    'Regulon specificity score (RSS)',
    ha='center',
    va='center',
    rotation='vertical',
    size='x-large',
)
plt.tight_layout()
plt.rcParams.update({
    'figure.autolayout': True,
    'figure.titlesize': 'large',
    'axes.labelsize': 'medium',
    'axes.titlesize': 'large',
    'xtick.labelsize': 'medium',
    'ytick.labelsize': 'medium'
})
plt.savefig(
    Path.joinpath(WORKDIR, "results", "scenic-plot", f"{idx}.rss.pdf"),
    dpi=150,
    bbox_inches="tight",
)
