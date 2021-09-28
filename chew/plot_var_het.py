"""Plotting of var(het)."""

import os
import json

from logzero import logger
import numpy as np
import plotly.express as px
import plotly.io as pio
import pandas as pd
import vcfpy

from .fingerprint import load_fingerprints

def run(args):
    fps = load_fingerprints(args.fps)
    plot_data = {"sample": [], "metric": [], "value": []}
    for sample in fps.keys():
        is_alt = fps[sample][0][1]
        is_alt_hom = fps[sample][0][2]
        het_vars = fps[sample][1][is_alt & ~is_alt_hom]
        plot_data["sample"].append(sample)
        plot_data["metric"].append("var(het)")
        plot_data["value"].append(np.var(het_vars))
    logger.info("Plotting to %s", args.out_html)
    plot_df = pd.DataFrame.from_dict(plot_data)
    fig = px.box(plot_df, x="metric", y="value", points="all")
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)
    logger.info("All done. Have a nice day!")
