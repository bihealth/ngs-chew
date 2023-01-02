"""Plotting of var(het)."""

import attrs
from logzero import logger
import pandas as pd
import plotly.express as px
import plotly.io as pio


@attrs.frozen
class Config:
    verbosity: int
    title: str
    out_html: str
    stats_out: str


def run(config: Config):
    logger.info("Loading statistics...")
    df = pd.read_csv(config.stats_out, sep="\t", header=0, index_col=0)

    logger.info("Plotting to %s", config.out_html)
    print(df)
    fig = px.box(df, x="var_het", points="all")
    fig.update_layout(title=config.title)
    pio.write_html(fig, file=config.out_html)

    logger.info("All done. Have a nice day!")
