#! /usr/bin/env python3
# -*- encoding: utf-8 -*-

import sys
from pathlib import Path

import pandas as pd
from pyscenic.binarization import binarize

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")

idx = sys.argv[1]
auc_df = pd.read_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.auc.csv"),
    index_col=0,
    header=0,
)
bin_df, thresholds = binarize(auc_df)
bin_df.to_csv(Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.bin.csv"))
thresholds.to_csv(
    Path.joinpath(WORKDIR, "results", "scenic", f"{idx}.threshold.csv"),
    header=False,
)
