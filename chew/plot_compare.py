"""Plotting for comparison."""

import attrs
from logzero import logger
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio


@attrs.frozen
class Config:
    verbosity: int
    compare_out: str
    out_html: str
    title: str


def run(config: Config):
    logger.info("Loading matrix.")

    df = pd.read_csv(config.compare_out, sep="\t", header=0, index_col=0)

    fig = go.Figure(data=go.Heatmap(z=df.to_numpy(), x=df.columns, y=df.index))
    fig.update_layout(title=config.title)
    pio.write_html(fig, file=config.out_html)

    logger.info("All done. Have a nice day!")
