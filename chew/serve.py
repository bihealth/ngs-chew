"""Implementation of the serve command"""


import typing

import attrs
import cattrs
from dash import Dash, dcc, html
from logzero import logger
import pandas as pd
import plotly.express as px
from tqdm import tqdm

from chew.stats import (
    SampleStats,
    compute_chrx_het_hom,
    compute_sample_stats,
    extract_header,
    load_fingerprint_all,
)


@attrs.frozen
class Config:
    verbosity: int
    #: Path to cohort-wide PED file
    cohort_ped: str
    #: List of paths to fingerprint ``.npz`` files
    fingerprints: typing.List[str]
    #: Suffix to strip from sample names
    strip_suffix: str


def load_all_stats(config: Config):
    logger.info("Loading all fingerprint files")
    all_stats = []
    for container in map(load_fingerprint_all, tqdm(config.fingerprints)):
        all_stats.append(compute_sample_stats(container))
    return all_stats


def run(config: Config):
    logger.info("Running report server...")
    logger.info("")
    all_stats = load_all_stats(config)
    df = pd.DataFrame.from_records(
        list(map(cattrs.unstructure, all_stats)),
    )
    df["sample_name"] = df["sample_name"].str.replace(config.strip_suffix, "")
    print(df)

    app = Dash(__name__)
