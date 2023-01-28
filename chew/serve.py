"""Implementation of the serve command"""


import csv
import typing

import attrs
import cattrs
from dash import Dash, Input, Output, dash_table, dcc, html
from dash.dash_table.Format import Format, Scheme
import dash_bootstrap_components as dbc
from logzero import logger
import numpy as np
import pandas as pd
import plotly.express as px
from tqdm import tqdm

from chew.common import (
    CHROM_LENS_GRCH37,
    CHROM_LENS_GRCH38,
    PedigreeMember,
    pedigree_member_from_tsv,
)
from chew.stats import compute_sample_stats, extract_header, load_fingerprint_all


@attrs.frozen
class ChromDosage:
    sample_name: str
    chr_1: float
    chr_2: float
    chr_3: float
    chr_4: float
    chr_5: float
    chr_6: float
    chr_7: float
    chr_8: float
    chr_9: float
    chr_10: float
    chr_11: float
    chr_12: float
    chr_13: float
    chr_14: float
    chr_15: float
    chr_16: float
    chr_17: float
    chr_18: float
    chr_19: float
    chr_20: float
    chr_21: float
    chr_22: float
    chr_x: float
    chr_y: float


def compute_chrom_dosages(container) -> ChromDosage:
    samtools_idxstats = container["samtools_idxstats"].tolist()
    counts = {}
    for line in samtools_idxstats.splitlines():
        arr = line.split("\t")
        counts[arr[0]] = int(arr[2])
    total = sum(counts.values())

    def to_key(s: str) -> str:
        if s.startswith("chr"):
            return s.replace("chr", "chr_").lower()
        else:
            return f"chr_{s}".lower()

    kwargs = {
        to_key(k): v / total
        for k, v in counts.items()
        if k in CHROM_LENS_GRCH37 or k in CHROM_LENS_GRCH38
    }
    header = extract_header(container)
    return ChromDosage(sample_name=header.sample, **kwargs)


@attrs.frozen
class Config:
    verbosity: int
    #: Path to cohort-wide PED file
    cohort_ped: str
    #: Optional path to annotation TSV file
    annos_tsv: typing.Optional[str]
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


def load_all_chrom_dosages(config: Config):
    logger.info("Loading all chrom dosage data")
    chrom_dosages = []
    for container in map(load_fingerprint_all, tqdm(config.fingerprints)):
        chrom_dosages.append(compute_chrom_dosages(container))
    return chrom_dosages


def load_ped_file(config: Config) -> typing.List[PedigreeMember]:
    logger.info("Load pedigree file")
    result = []
    with open(config.cohort_ped, "rt") as inputf:
        reader = csv.reader(inputf, delimiter="\t")
        for row in reader:
            result.append(pedigree_member_from_tsv(row))
    return result


COLUMN_LABELS = {
    "sample_name": "Sample",
    "sample": "Sample",
    "father": "Father",
    "mother": "Mother",
    "affected": "Affected",
    "sex": "Sex",
    "release": "Genome",
    "hom_refs": "# of 0/0 SNPs",
    "hets": "# of 0/1 SNPs",
    "hom_alts": "# of 1/1 SNPs",
    "mask_ones": "covered SNPs",
    "var_het": "var(het)",
    "chrx_het_hom": "het/hom calls on chrX",
    "chrx_frac": "fraction of chrX reads",
    "chry_frac": "fraction of chrY reads",
}

COLUM_FORMATS = {
    "var_het": {"type": "numeric", "format": Format(precision=4, scheme=Scheme.fixed)},
    "chrx_het_hom": {"type": "numeric", "format": Format(precision=3, scheme=Scheme.fixed)},
    "chrx_frac": {"type": "numeric", "format": Format(precision=3, scheme=Scheme.fixed)},
    "chry_frac": {"type": "numeric", "format": Format(precision=3, scheme=Scheme.fixed)},
}

SCATTER_DIMS = [
    "var_het",
    "chrx_het_hom",
    "chrx_frac",
    "chry_frac",
    "hom_refs",
    "sex",
    "hets",
    "hom_alts",
    "mask_ones",
]

SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}


def build_select_div(
    field_id: str,
    label: str,
    initial_value: str,
    anno_dims,
) -> html.Div:
    return html.Div(
        [
            dbc.Select(
                options=[
                    {
                        "label": COLUMN_LABELS.get(dim, dim),
                        "value": dim,
                    }
                    for dim in SCATTER_DIMS + anno_dims
                ],
                value=initial_value,
                id=field_id,
            ),
            html.Label(
                [label],
                htmlFor=field_id,
            ),
        ],
        className="form-floating col-12",
    )


def unpivot(df, key="key", variable="variable", value="value"):
    n, k = df.shape
    data = {
        value: df.to_numpy().ravel("F"),
        variable: np.asarray(df.columns).repeat(n),
        key: np.tile(np.asarray(df.index), k),
    }
    return pd.DataFrame(data, columns=[key, variable, value])


def run(config: Config):
    logger.info("Running report server...")
    logger.info("")
    all_stats = load_all_stats(config)
    df_stats = pd.DataFrame.from_records(
        list(map(cattrs.unstructure, all_stats)),
    )
    all_chrom_dosages = load_all_chrom_dosages(config)
    df_chrom_dosages = pd.DataFrame.from_records(list(map(cattrs.unstructure, all_chrom_dosages)))
    df_stats["sample_name"] = df_stats["sample_name"].str.replace(config.strip_suffix, "")

    # Load pedigree
    logger.info("Loading pedigree file...")
    pedigree = load_ped_file(config)
    df_ped = pd.DataFrame.from_records(list(map(cattrs.unstructure, pedigree)))
    df_ped["name"] = df_ped["name"].str.replace(config.strip_suffix, "")
    df_ped["father"] = df_ped["father"].str.replace(config.strip_suffix, "")
    df_ped["mother"] = df_ped["mother"].str.replace(config.strip_suffix, "")

    df = df_ped.set_index("name").join(df_stats.set_index("sample_name")).reset_index()

    if config.annos_tsv:
        logger.info("Loading annotations TSV %s", config.annos_tsv)
        df_anno = pd.read_csv(config.annos_tsv, sep="\t")
        df_anno_0 = df_anno.columns[0]
        df_anno[df_anno_0] = df_anno[df_anno_0].str.replace(config.strip_suffix, "")
        df = df.set_index("name").join(df_anno.set_index(df_anno_0)).reset_index()
        anno_dims = df_anno.columns.to_list()[1:]
    else:
        anno_dims = []

    logger.info("Data frame looks like %s", df)

    app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)
    sidebar = html.Div(
        [
            html.H2("ngs-chew", className="display-6"),
            html.Hr(),
            dbc.Nav(
                [
                    dbc.NavLink("Stats Table", href="/", active="exact"),
                    dbc.NavLink("Stats Plots: X vs. Y", href="/stats-plots", active="exact"),
                    dbc.NavLink("Dosage per Chrom.", href="/chrom-dosage-plots", active="exact"),
                ],
                vertical=True,
                pills=True,
            ),
        ],
        style=SIDEBAR_STYLE,
    )

    content = html.Div(id="page-content", style=CONTENT_STYLE)

    app.layout = html.Div([dcc.Location(id="url"), sidebar, content])

    @app.callback(
        Output("stats-scatter-plot", "children"),
        [
            Input("stats-scatter-select-dim-h", "value"),
            Input("stats-scatter-select-dim-v", "value"),
            Input("stats-scatter-select-dim-color", "value"),
            Input("stats-scatter-select-dim-symbol", "value"),
            Input("stats-scatter-select-marker-size", "value"),
        ],
    )
    def render_scatter_plot(dim_h, dim_v, color, symbol, marker_size):
        if np.issubdtype(df[dim_h].dtype, np.number):
            px_func = px.scatter
            px_kwargs = {"symbol": symbol}
        else:
            px_func = px.strip
            px_kwargs = {}
        fig = px_func(
            df,
            x=dim_h,
            y=dim_v,
            color=color,
            hover_data=["name", "sex"],
            labels={
                dim_h: COLUMN_LABELS.get(dim_h, dim_h),
                dim_v: COLUMN_LABELS.get(dim_v, dim_v),
            },
            **px_kwargs,
        )
        fig.update_traces(marker={"size": int(marker_size)})
        return dcc.Graph(figure=fig)

    @app.callback(
        Output("chrom-dosage-plots", "children"),
        [
            Input("chrom-dosage-select-dim-color", "value"),
            Input("chrom-dosage-select-dim-symbol", "value"),
            Input("chrom-dosage-select-marker-size", "value"),
        ],
    )
    def render_chrom_dosage_plot(color, symbol, marker_size):
        df_chrom_dosage_plot_tmp = unpivot(
            df_chrom_dosages.set_index("sample_name"),
            "sample_name",
            "chrom",
            "dosage",
        )
        df_chrom_dosage_plot_tmp["chrom"] = df_chrom_dosage_plot_tmp["chrom"].str.upper()
        df_chrom_dosage_plot_tmp["chrom"] = df_chrom_dosage_plot_tmp["chrom"].str.replace(
            "CHR_", ""
        )
        df_chrom_dosage_plot = (
            df_ped.set_index("name")
            .join(df_chrom_dosage_plot_tmp.set_index("sample_name"))
            .reset_index()
            .rename(columns={"index": "name"})
        )
        if config.annos_tsv:
            logger.info("Loading annotations TSV %s", config.annos_tsv)
            df_anno = pd.read_csv(config.annos_tsv, sep="\t")
            df_anno_0 = df_anno.columns[0]
            df_anno[df_anno_0] = df_anno[df_anno_0].str.replace(config.strip_suffix, "")
            df_chrom_dosage_plot = (
                df_chrom_dosage_plot.set_index("name")
                .join(df_anno.set_index(df_anno_0))
                .reset_index()
                .rename(columns={"index": "name"})
            )
        fig = px.strip(
            df_chrom_dosage_plot,
            x="chrom",
            y="dosage",
            color=color,
            hover_data=["name", "dosage"],
        )
        fig.update_traces(marker={"size": int(marker_size)})
        return dcc.Graph(figure=fig)

    @app.callback(Output("page-content", "children"), [Input("url", "pathname")])
    def render_page_content(pathname):
        if pathname == "/":
            display_df = df.rename(columns=COLUMN_LABELS)
            columns = [
                {
                    "name": COLUMN_LABELS.get(col, col),
                    "id": COLUMN_LABELS.get(col, col),
                    **COLUM_FORMATS.get(col, {}),
                }
                for col in df.columns
            ]
            return dash_table.DataTable(
                display_df.to_dict("records"),
                columns,
                sort_action="native",
                filter_action="native",
                page_action="native",
                page_current=0,
                page_size=20,
            )
        elif pathname == "/stats-plots":
            row_select = html.Form(
                [
                    build_select_div(
                        "stats-scatter-select-dim-h", "horizontal axis", "sex", anno_dims
                    ),
                    build_select_div(
                        "stats-scatter-select-dim-v", "vertical axis", "chrx_het_hom", anno_dims
                    ),
                    build_select_div("stats-scatter-select-dim-color", "color", "sex", anno_dims),
                    build_select_div("stats-scatter-select-dim-symbol", "shape", "sex", anno_dims),
                    html.Div(
                        [
                            dbc.Select(
                                options=[
                                    {"label": str(2 * i), "value": str(2 * i)} for i in range(1, 11)
                                ],
                                value="6",
                                id="stats-scatter-select-marker-size",
                            ),
                            html.Label(
                                ["size"],
                                htmlFor="stats-scatter-select-marker-size",
                            ),
                        ],
                        className="form-floating col-12",
                    ),
                ],
                className="row row-cols-lg-auto g-3 align-items-center",
            )
            row_plot = html.Div(
                [
                    html.Div(
                        [
                            dbc.Placeholder(xs=6),
                            html.Br(),
                            dbc.Placeholder(className="w-75"),
                            html.Br(),
                            dbc.Placeholder(style={"width": "25%"}),
                        ],
                        id="stats-scatter-plot",
                        className="col-12",
                    )
                ],
                className="row",
            )
            return html.Div(
                [
                    row_select,
                    row_plot,
                ]
            )
        elif pathname == "/chrom-dosage-plots":
            row_select = html.Form(
                [
                    build_select_div("chrom-dosage-select-dim-color", "color", "sex", anno_dims),
                    build_select_div("chrom-dosage-select-dim-symbol", "shape", "sex", anno_dims),
                    html.Div(
                        [
                            dbc.Select(
                                options=[
                                    {"label": str(2 * i), "value": str(2 * i)} for i in range(1, 11)
                                ],
                                value="6",
                                id="chrom-dosage-select-marker-size",
                            ),
                            html.Label(
                                ["size"],
                                htmlFor="chrom-dosage-select-marker-size",
                            ),
                        ],
                        className="form-floating col-12",
                    ),
                ],
                className="row row-cols-lg-auto g-3 align-items-center",
            )
            row_plot = html.Div(
                [
                    html.Div(
                        [
                            dbc.Placeholder(xs=6),
                            html.Br(),
                            dbc.Placeholder(className="w-75"),
                            html.Br(),
                            dbc.Placeholder(style={"width": "25%"}),
                        ],
                        id="chrom-dosage-plots",
                        className="col-12",
                    )
                ],
                className="row",
            )
            return html.Div(
                [
                    row_select,
                    row_plot,
                ]
            )
        # If the user tries to reach a different page, return a 404 message
        return html.Div(
            [
                html.H1("404: Not found", className="text-danger"),
                html.Hr(),
                html.P(f"The pathname {pathname} was not recognised..."),
            ],
            className="p-3 bg-light rounded-3",
        )

    app.run_server(
        dev_tools_hot_reload=True,
        debug=True,
    )
