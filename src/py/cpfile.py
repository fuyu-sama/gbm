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

# %%
import shutil
from pathlib import Path

import pandas as pd

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")

# %%
gtf_names = {}
gtf_path = Path.joinpath(Path.home(), "refgenome", "spaceranger",
                         "refdata-gex-GRCh38-2020-A", "genes", "genes.gtf")
with open(gtf_path) as f:
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip()
        line = line.split("\t")[8].split("; ")
        for item in line:
            key, value = item.split(" ")
            if key == "gene_id":
                gene_id = value.strip('"')
            elif key == "gene_name":
                symbol = value.strip('"')
        gtf_names[gene_id] = symbol

probe_path = Path.joinpath(
    Path.home(),
    "refgenome",
    "spaceranger",
    "Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv",
)
with open(probe_path) as f:
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip()
        if line[0:7] == "gene_id":
            continue
        line = line.split(",")[2].split("|")
        gene_id = line[0]
        symbol = line[1]
        gtf_names[gene_id] = symbol

# %%
idx = "21B-603-5"
outs_dir = Path.joinpath(WORKDIR, "spaceranger", idx, "outs")
count_path = Path.joinpath(outs_dir, "raw_feature_bc_matrix", f"{idx}.raw.csv")
coor_path = Path.joinpath(outs_dir, "spatial", "tissue_positions.csv")

count_dst = Path.joinpath(WORKDIR, "Data", "scale_df", "raw", f"{idx}.csv")
coor_dst = Path.joinpath(WORKDIR, "Data", "coor_df", f"{idx}.csv")

count_df = pd.read_csv(count_path, index_col=0, header=0)
count_df.index = [gtf_names[i] for i in count_df.index]
count_df.to_csv(count_dst)
shutil.copy(coor_path, coor_dst)

# %%
idx = "22F-10823-3"
outs_dir = Path.joinpath(WORKDIR, "spaceranger", idx, "outs")
count_path = Path.joinpath(outs_dir, "filtered_feature_bc_matrix",
                           f"{idx}.filtered.csv")
coor_path = Path.joinpath(outs_dir, "spatial", "tissue_positions.csv")

count_dst = Path.joinpath(WORKDIR, "Data", "scale_df", "raw", f"{idx}.csv")
coor_dst = Path.joinpath(WORKDIR, "Data", "coor_df", f"{idx}.csv")

count_df = pd.read_csv(count_path, index_col=0, header=0)
count_df.index = [gtf_names[i] for i in count_df.index]
count_df.to_csv(count_dst)
shutil.copy(coor_path, coor_dst)

# %%
idx = "22F-21576-1"
outs_dir = Path.joinpath(WORKDIR, "spaceranger", idx, "outs")
count_path = Path.joinpath(outs_dir, "filtered_feature_bc_matrix",
                           f"{idx}.filtered.csv")
coor_path = Path.joinpath(outs_dir, "spatial", "tissue_positions.csv")

count_dst = Path.joinpath(WORKDIR, "Data", "scale_df", "raw", f"{idx}.csv")
coor_dst = Path.joinpath(WORKDIR, "Data", "coor_df", f"{idx}.csv")

count_df = pd.read_csv(count_path, index_col=0, header=0)
count_df.index = [gtf_names[i] for i in count_df.index]
count_df.to_csv(count_dst)
shutil.copy(coor_path, coor_dst)

# %%
idx = "22F-23738-2"
outs_dir = Path.joinpath(WORKDIR, "spaceranger", idx, "outs")
count_path = Path.joinpath(outs_dir, "filtered_feature_bc_matrix",
                           f"{idx}.filtered.csv")
coor_path = Path.joinpath(outs_dir, "spatial", "tissue_positions.csv")

count_dst = Path.joinpath(WORKDIR, "Data", "scale_df", "raw", f"{idx}.csv")
coor_dst = Path.joinpath(WORKDIR, "Data", "coor_df", f"{idx}.csv")

count_df = pd.read_csv(count_path, index_col=0, header=0)
count_df.index = [gtf_names[i] for i in count_df.index]
count_df.to_csv(count_dst)
shutil.copy(coor_path, coor_dst)
