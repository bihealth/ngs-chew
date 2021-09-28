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

# old function (to compare and correct any errors):
"""
def load_het_aabs(args):
    het_aabs = {}
    if args.fps:
        fps = load_fingerprints(args.fps)

    if args.vcf:
        for path_vcf in args.vcf:
            logger.info("- %s", path_vcf)
            with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
                for sample in vcf_reader.header.samples.names:
                    het_aabs.setdefault(sample, {})
                for record in vcf_reader:
                    key = "%s:%s" % (record.CHROM, record.POS)
                    for call in record.calls:
                        if call.data["GT"].replace("|", "/") in ("0/1", "1/0"):
                            if sum(call.data.get("AD", [0])):
                                het_aabs[call.sample][key] = call.data["AD"][0] / sum(call.data["AD"])
                            else:
                                het_aabs[call.sample][key] = 0
        return het_aabs
"""

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
    # debug
    print(plot_df)
    fig = px.box(plot_df, x="metric", y="value", points="all")
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)

    logger.info("All done. Have a nice day!")
