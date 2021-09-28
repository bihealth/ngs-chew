"""Plotting of alternate allele balance."""

import os
import json

from logzero import logger
import numpy as np
import plotly.express as px
import plotly.io as pio
import pandas as pd
import vcfpy

from .fingerprint import load_fingerprints

"""def load_aabs(args):
    aabs = {}
    for path_vcf in args.vcf:
        logger.info("- %s", path_vcf)
        with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
            for sample in vcf_reader.header.samples.names:
                aabs.setdefault(sample, {})
            for record in vcf_reader:
                key = "%s:%s" % (record.CHROM, record.POS)
                for call in record.calls:
                    if sum(call.data.get("AD", [0])):
                        aabs[call.sample][key] = call.data["AD"][0] / sum(call.data["AD"])
                    else:
                        aabs[call.sample][key] = 0

    keys = None
    for mapping in aabs.values():
        if keys is None:
            keys = set(mapping.keys())
        else:
            keys &= set(mapping.keys())

    for mapping in aabs.values():
        for key in list(mapping.keys()):
            if key not in keys:
                del mapping[key]

    return aabs


def run(args):
    logger.info("Loading AAB values.")

    if args.aab_cache and os.path.exists(args.aab_cache):
        logger.info("Loading AABs from %s", args.aab_cache)
        with open(args.aab_cache, "rt") as inputf:
            aabs = json.load(inputf)
    else:
        logger.info("Collecting AABs")
        aabs = load_aabs(args)
        logger.info("Writing AABs to %s", args.aab_cache)
        if args.aab_cache:
            with open(args.aab_cache, "wt") as outputf:
                json.dump(aabs, outputf)

    logger.info("Preparing plot...")
    df = pd.DataFrame.from_dict(aabs)

    plot_data = {"sample": [], "bin": [], "count": []}
    for c in df.columns:
        hist = np.histogram(df[c], 50, (0.0, 1.0))
        plot_data["count"] += hist[0].tolist()
        plot_data["sample"] += [c] * hist[0].shape[0]
        plot_data["bin"] += hist[1].tolist()[:-1]

    logger.info("Plotting to %s", args.out_html)
    plot_df = pd.DataFrame.from_dict(plot_data)
    fig = px.line(plot_df, x="bin", y="count", color="sample")
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)

    logger.info("All done. Have a nice day!")"""

def run(args):
    fps = load_fingerprints(args.fps)
    df = pd.DataFrame.from_dict({key: fps[key][1] for key in fps.keys()})
    df.drop(df[df.sum(axis=1) == 0].index,inplace=True)

    # masking - remove targets not in a given set of keys K
    # - load sites
    # - rangeindex sites
    # - filter rangeindex for sites given in K -> s
    # - df.index = intersection(df.index,s)

    plot_data = {"sample": [], "bin": [], "count": []}
    for c in df.columns:
        hist = np.histogram(df[c], 50, (0.0, 1.0))
        plot_data["count"] += hist[0].tolist()
        plot_data["sample"] += [c] * hist[0].shape[0]
        plot_data["bin"] += hist[1].tolist()[:-1]

    logger.info("Plotting to %s", args.out_html)
    plot_df = pd.DataFrame.from_dict(plot_data)
    fig = px.line(plot_df, x="bin", y="count", color="sample")
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)

    logger.info("All done. Have a nice day!")
