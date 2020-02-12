"""Plotting of var(het)."""

import os
import json

from logzero import logger
import numpy as np
import plotly.express as px
import plotly.io as pio
import pandas as pd
import vcfpy


def load_het_aabs(args):
    het_aabs = {}

    for path_vcf in args.vcf:
        logger.info("- %s", path_vcf)
        with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
            for sample in vcf_reader.header.samples.names:
                het_aabs.setdefault(sample, {})
            for record in vcf_reader:
                key = "%s:%s" % (record.CHROM, record.POS)
                for call in record.calls:
                    if call.data['GT'].replace('|', '/') in ('0/1', '1/0'):
                        if sum(call.data.get("AD", [0])):
                            het_aabs[call.sample][key] = call.data["AD"][0] / sum(call.data["AD"])
                        else:
                            het_aabs[call.sample][key] = 0

    return het_aabs


def run(args):
    logger.info("Loading HET AAB values.")

    if args.var_het_cache and os.path.exists(args.var_het_cache):
        logger.info("Loading AABs from %s", args.var_het_cache)
        with open(args.var_het_cache, "rt") as inputf:
            het_aabs = json.load(inputf)
    else:
        logger.info("Collecting AABs")
        het_aabs = load_het_aabs(args)
        logger.info("Writing AABs to %s", args.var_het_cache)
        if args.var_het_cache:
            with open(args.var_het_cache, "wt") as outputf:
                json.dump(het_aabs, outputf)

    logger.info("Preparing plot...")
    df = pd.DataFrame.from_dict(het_aabs)

    plot_data = {"sample": [], "metric": [], "value": []}
    for c in df.columns:
        plot_data["sample"].append(c)
        aabs = df[c][~np.isnan(df[c])]
        plot_data["metric"].append("var(het)")
        plot_data["value"].append(np.var(aabs))

    logger.info("Plotting to %s", args.out_html)
    plot_df = pd.DataFrame.from_dict(plot_data)
    print(plot_df)
    fig = px.box(plot_df, x="metric", y="value", points="all")
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)

    logger.info("All done. Have a nice day!")
