#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# %%
import sys
from pathlib import Path

import gepia

WORKDIR = Path.joinpath(Path.home(), "workspace", "gbm")

bp = gepia.boxplot()
bp.setParam("dataset", ["GBM", "LGG"])
bp.setOutDir(str(Path.joinpath(WORKDIR, "results", "GEPIA")) + "/")

bp.setParam("signature", sys.argv[1])
bp.query()
