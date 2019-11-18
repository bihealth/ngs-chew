"""Plotting for comparison."""

from logzero import logger
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio


def run(args):
    logger.info("Loading matrix.")

    df = pd.read_csv(args.stats, sep="\t", header=0, index_col=0)

    fig = go.Figure(data=go.Heatmap(z=df.to_numpy(), x=df.columns, y=df.index))
    fig.update_layout(title=args.title)
    pio.write_html(fig, file=args.out_html)

    logger.info("All done. Have a nice day!")
