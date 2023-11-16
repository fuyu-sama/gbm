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

# %% environment config
import pickle
from collections import Counter
from pathlib import Path
from PIL import Image
from multiprocessing import Pool

import pandas as pd
import matplotlib.pyplot as plt

import svgbit as sb

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")
idx_full = ["21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2"]
he_path = Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                        "V43J31-122_D1_21B-603-5.tif")
Image.MAX_IMAGE_PIXELS = 933120000
he_dir = Path.joinpath(WORKDIR, "Data", "ST", "brightfield")
he_path = {
    "21B-603-5": Path.joinpath(he_dir, "V43J31-122_D1_21B-603-5.tif"),
    "22F-10823-3": Path.joinpath(he_dir, "V43J31-122_A1_22F-10823-3.tif"),
    "22F-21576-1": Path.joinpath(he_dir, "V43J31-095_A1_22F-21576-1.tif"),
    "22F-23738-2": Path.joinpath(he_dir, "V43J31-095_D1_22F-23738-2.tif"),
}

# %%
datasets = {}
for idx in idx_full:
    with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
              "rb") as f:
        datasets[idx] = pickle.load(f)

# %%
genes = []
for idx in idx_full:
    gene_pairs = pd.read_csv(
        Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-gene_pairs.csv"),
        index_col=0,
        header=0,
    )
    gene_pairs = gene_pairs[gene_pairs["gene_1"] == "VIM"]
    gene_pairs = gene_pairs[gene_pairs["colocalization_degree"] != "low"]
    [genes.append(i) for i in gene_pairs["gene_2"]]
genes = Counter(genes)
genes_df = pd.DataFrame.from_dict(genes, orient="index")
genes_df.to_csv("count.csv")

save_dir = Path.joinpath(WORKDIR, "results", "svgbit", "colocalization_map")
if not save_dir.exists():
    save_dir.mkdir()


def draw(idx):
    for gene in genes:
        if genes[gene] < 2:
            continue
        fig = sb.plot.hotspot_colocalization_map(
            datasets[idx],
            ("VIM", gene),
            Path.joinpath(save_dir, f"VIM-{gene}-{idx}.jpg"),
            he_path[idx],
            s=8,
        )
        plt.close(fig)


pool = Pool(4)
pool.map(draw, idx_full)
pool.close()
pool.join()
