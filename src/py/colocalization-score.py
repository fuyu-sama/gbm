#! /usr/bin/env python3
# -*- encoding: utf-8 -*-

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

# %%
import pickle
from pathlib import Path
from PIL import Image

import numpy as np
import matplotlib.pyplot as plt

import svgbit as sb

plt.rcParams.update({'font.size': 18})
WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")
idx_full = ["21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2"]
Image.MAX_IMAGE_PIXELS = 933120000

he_path = {
    "21B-603-5":
    Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                  "V43J31-122_D1_21B-603-5.tif"),
    "22F-10823-3":
    Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                  "V43J31-122_A1_22F-10823-3.tif"),
    "22F-21576-1":
    Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                  "V43J31-095_A1_22F-21576-1.tif"),
    "22F-23738-2":
    Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                  "V43J31-095_D1_22F-23738-2.tif")
}

# %%
datasets = {}
for idx in idx_full:
    with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
              "rb") as f:
        datasets[idx] = pickle.load(f)

# %%
gene_1 = "CD58"
gene_2 = "CDC25C"
gene_2 = "CDK1"
n_spots = {
    f"{gene_1} only": [],
    f"{gene_2} only": [],
    f"{gene_1} & {gene_2}": [],
}
for idx in datasets:
    d = datasets[idx]
    hotspot_df = d.hotspot_df.sparse.to_dense()
    gene_1_series = hotspot_df[gene_1]
    gene_2_series = hotspot_df[gene_2]
    co_series = gene_1_series + gene_2_series
    total_co = len(co_series[co_series == 2])
    gene_1_only = len(gene_1_series[gene_1_series == 1]) - total_co
    gene_2_only = len(gene_2_series[gene_2_series == 1]) - total_co

    n_spots[f"{gene_1} only"].append(gene_1_only)
    n_spots[f"{gene_2} only"].append(gene_2_only)
    n_spots[f"{gene_1} & {gene_2}"].append(total_co)

fig, ax = plt.subplots(figsize=(10, 10))

width = 0.25  # the width of the bars
multiplier = 0
colors = [
    "tab:green",
    "tab:blue",
    "tab:orange",
    "tab:pink",
    "tab:cyan",
]
x = np.arange(len(idx_full))
ymax = 0
for attribute, measurement in n_spots.items():
    ymax = max(ymax, max(measurement))
    offset = width * multiplier
    rects = ax.bar(
        x + offset,
        measurement,
        width,
        label=attribute,
        color=colors[multiplier],
    )
    ax.bar_label(rects, padding=3)
    multiplier += 1

ax.set_ylabel('Number of hotspots')
ax.set_xticks(x + width, idx_full)
ax.legend(loc='upper center', ncols=1)
ax.set_ylim(0, ymax + 50)

fig.savefig(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{gene_1}&{gene_2}.svg"))

for idx in datasets:
    d = datasets[idx]
    save_path = Path.joinpath(WORKDIR, "results", "svgbit",
                              f"{idx}-{gene_1}&{gene_2}.svg")
    sb.core.plot._hotspot_colocalization_map(
        d.hotspot_df,
        d.coordinate_df[["Y", "X"]],
        (gene_1, gene_2),
        save_path,
        he_path[idx],
        colors=colors,
    )

    save_path = Path.joinpath(WORKDIR, "results", "svgbit",
                              f"{idx}-{gene_1}.svg")
    sb.core.plot._hotspot_expression(
        d.hotspot_df,
        d.coordinate_df[["Y", "X"]],
        gene_1,
        save_path,
        he_path[idx],
    )

    save_path = Path.joinpath(WORKDIR, "results", "svgbit",
                              f"{idx}-{gene_2}.svg")
    sb.core.plot._hotspot_expression(
        d.hotspot_df,
        d.coordinate_df[["Y", "X"]],
        gene_2,
        save_path,
        he_path[idx],
    )
