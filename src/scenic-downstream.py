#! /usr/bin/env python3
# -*- encoding: utf-8 -*-

# %%
import sys
from math import floor, ceil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")

idx = sys.argv[1]


def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None, regulon_name=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    if regulon_name is not None:
        regulon_name_input = regulon_name
    else:
        regulon_name_input = None
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, ".")
    ax.set_ylim([floor(data.min() * 100.0) / 100.0, ceil(data.max() * 100.0) / 100.0])
    ax.set_ylabel("RSS")
    ax.set_xlabel("Regulon")
    ax.set_title(cell_type)
    ax.set_xticklabels([])

    font = {
        "color": "red",
        "weight": "normal",
        "size": 4,
    }

    if top_n > 0:
        for idx, (regulon_name, rss_val) in enumerate(
                zip(data[0:top_n].index, data[0:top_n].values)
        ):
            ax.plot([idx, idx], [rss_val, rss_val], "r.")
            ax.text(
                idx + (max_n / 25),
                rss_val,
                regulon_name,
                fontdict=font,
                horizontalalignment="left",
                verticalalignment="center",
            )

    if regulon_name_input is not None:
        for regulon_name in regulon_name_input:
            if "(+)" not in regulon_name:
                regulon_name = regulon_name + "(+)"
            idx = data.index.get_loc(regulon_name)
            rss_val = data[regulon_name]
            ax.plot([idx, idx], [rss_val, rss_val], "r.")
            ax.text(
                idx + (max_n / 25),
                rss_val,
                regulon_name,
                fontdict=font,
                horizontalalignment="left",
                verticalalignment="center",
            )


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
idents_sub = idents_sub[idents_sub["x"] != "Blood vessel rich area"]
idents_sub["idx"] = [i.split("_")[0] for i in idents_sub.index]
idents_sub["level"] = np.where(
    idents_sub["idx"].isin(["21B-603-5", "22F-10823-3"]),
    "IV",
    "GBM",
)

# idents_sub_2 = [i for i in idents.index if i not in idents_sub.index]
# idents_sub_2 = idents.loc[idents_sub_2, ]
# idents_sub_3 = idents.copy()
# idents_sub_3["idx"] = [i.split("_")[0] for i in idents_sub_3.index]
# idents_sub_3["level"] = np.where(
# idents_sub_3["idx"].isin(["21B-603-5", "22F-10823-3"]),
# "IV",
# "GBM",
# )
# for i in idents_sub_3.index:
# if i in idents_sub_2.index:
# idents_sub_3.loc[i, "level"] = "para"
tops = {
    "GBM": ["FOS", "SOX8", "JUNB", "JUN", "TBX2",
            "FOXA3", "FOSB", "HOXC5", "MAFB", "FOXF2"],
    "IV": ["JUN", "FOS", "JUNB", "CEBPD", "ETS1",
           "MAFB", "TBX2", "LTF", "NR5A2", "TEAD4"]
}

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
    plot_rss(
        rss_level,
        c,
        top_n=0,
        max_n=None,
        ax=ax,
        regulon_name=tops[c]
    )
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
    Path.joinpath(WORKDIR, "results", "scenic-plot", f"{idx}.rss.level.1.pdf"),
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
