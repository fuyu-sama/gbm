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
                        "V43J31-095_D1_22F-23738-2.tif")
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
    coor_df = coor_df[["pxl_row_in_fullres", "pxl_col_in_fullres"]]

    d = sb.STDataset(count_df, coor_df)
    d = sb.filters.low_variance_filter(d)
    d = sb.filters.high_expression_filter(d)
    # d = sb.normalizers.logcpm_normalizer(d)
    return d


# %%
idx = "22F-23738-2"
d = read_exp(idx)
sb.run(d, cores=5, n_svgs=500, n_svg_clusters=2)

d.hotspot_df.T.to_csv(
    Path.joinpath(WORKDIR, "Data", "scale_df", "hotspot", f"{idx}.csv"))
d.AI.sort_values(ascending=False).to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-AI.csv"))
d.svg_cluster.to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-svg-cluster.csv"))
sb.plot.svg_heatmap(
    d, Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-heatmap.svg"))
sb.plot.spot_type_map(
    d, Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-typemap.svg"))

with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
          "wb") as f:
    pickle.dump(d, f)

# %%
if True:
    idx = "22F-23738-2"
    with open(Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}.pickle"),
              "rb") as f:
        d = pickle.load(f)

# %%
selected_genes = ["CDC27", "DBF4", "ANAPC4", "STAG2", "RBL2", "TTK", "ORC2",
                  "ORC4", "CCNH", "ORC3", "RB1", "CCNA2", "BUB1B", "RAD21",
                  "ORC5", "BUB1", "CDK1", "HDAC2"]
selected_genes = set(selected_genes + list(d.AI.sort_values()[-500:].index))

gene_pairs = sb.find_combinations(
    d,
    center_spots=d.spots,
    selected_genes=selected_genes,
    use_neighbor=False,
)
gene_pairs.to_csv(
    Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-gene_pairs.csv"))

# %%
colors = [
    "tab:green", "tab:blue", "tab:orange", "tab:pink", "tab:cyan",
]
save_path = Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-CD58&CDK1.svg")
sb.core.plot._hotspot_colocalization_map(
    d.hotspot_df,
    d.coordinate_df[["Y", "X"]],
    ("CD58", "CDK1"),
    save_path,
    he_path,
    colors=colors,
)

save_path = Path.joinpath(WORKDIR, "results", "svgbit", f"{idx}-CD58&CDC25C.svg")
sb.core.plot._hotspot_colocalization_map(
    d.hotspot_df,
    d.coordinate_df[["Y", "X"]],
    ("CD58", "CDC25C"),
    save_path,
    he_path,
    colors=colors,
)
