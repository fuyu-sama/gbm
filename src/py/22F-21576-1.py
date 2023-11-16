#! venv/bin/python
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
from pathlib import Path
from PIL import Image

import pandas as pd

import svgbit as sb

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")
idx_full = ["21B-603-5", "22F-10823-3", "22F-21576-1", "22F-23738-2"]
he_path = Path.joinpath(WORKDIR, "Data", "ST", "brightfield",
                        "V43J31-095_A1_22F-21576-1.tif")
Image.MAX_IMAGE_PIXELS = 933120000


# %% helper funcs
def read_exp(idx):
    count_df = pd.read_csv(
        Path.joinpath(WORKDIR, "Data", "scale_df", "raw", f"{idx}.csv"),
        index_col=0,
        header=0,
    ).T

    coor_df = pd.read_csv(
        Path.joinpath(WORKDIR, "Data", "coor_df", f"{idx}.csv"),
        index_col=0,
        header=0,
    )

    coor_df = coor_df.reindex(index=count_df.index)
    coor_df = coor_df[["pxl_col_in_fullres", "pxl_row_in_fullres"]]

    d = sb.STDataset(count_df, coor_df)
    d = sb.filters.low_variance_filter(d)
    # d = sb.filters.high_expression_filter(d)
    d = sb.normalizers.logcpm_normalizer(d)
    return d


# %%
idx = "22F-21576-1"
d = read_exp(idx)
sb.run(d, cores=5)

d.hotspot_df.T.to_csv(
    Path.joinpath(WORKDIR, "Data", "scale_df", "hotspot", f"{idx}.csv"))
d.AI.sort_values(ascending=False).to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-AI.csv"))
d.svg_cluster.to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-svg-cluster.csv"))
sb.plot.svg_heatmap(
    d, Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-heatmap.svg"),
    he_path)
sb.plot.spot_type_map(
    d, Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-typemap.svg"),
    he_path)

with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
          "wb") as f:
    pickle.dump(d, f)

# %%
idx = "22F-21576-1"
with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
          "rb") as f:
    d = pickle.load(f)

# %%
selected_genes = [
    "CD34", "CHKA", "EGFR", "ETFA", "GFAP", "IDH1", "MKI67", "MLH1", "MSH2",
    "MSH6", "MUC1", "NF1", "NF2", "NTRK2", "NTRK3", "OLIG2", "PDGFRA", "PMS2",
    "S100A1", "S100A4", "S100A6", "S100A8", "S100A9", "S100A10", "S100A11",
    "S100A13", "S100A16", "S100B", "S100PBP", "SMARCB1", "TP53", "VIM"
]
selected_genes = set(selected_genes + list(d.AI.sort_values()[-500:].index))

gene_pairs = sb.find_combinations(
    d,
    center_spots=d.spots,
    selected_genes=selected_genes,
    use_neighbor=False,
)
gene_pairs.to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-gene_pairs.csv"))
