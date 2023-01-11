"""Implementation of the serve command"""


import csv
import typing

import attrs
import cattrs
from dash import Dash, Input, Output, State, dcc, html
import dash_bootstrap_components as dbc
from dash import dash_table
from dash.dash_table.Format import Format, Scheme
from logzero import logger
import pandas as pd
import plotly.express as px
from tqdm import tqdm
from chew.common import PedigreeMember, pedigree_member_from_tsv

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


def run(config: Config):
    logger.info("Running report server...")
    logger.info("")
    all_stats = load_all_stats(config)
    df_stats = pd.DataFrame.from_records(
        list(map(cattrs.unstructure, all_stats)),
    )
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
                    dbc.NavLink("Stats Plots", href="/stats-plots", active="exact"),
                    dbc.NavLink("Page 2", href="/page-2", active="exact"),
                ],
                vertical=True,
                pills=True,
            ),
        ],
        style=SIDEBAR_STYLE,
    )

    nav = dbc.Nav(
        [
            dbc.NavItem(dbc.NavLink("Active", active=True, href="#")),
            dbc.NavItem(dbc.NavLink("A link", href="#")),
            dbc.NavItem(dbc.NavLink("Another link", href="#")),
            dbc.NavItem(dbc.NavLink("Disabled", disabled=True, href="#")),
            dbc.DropdownMenu(
                [dbc.DropdownMenuItem("Item 1"), dbc.DropdownMenuItem("Item 2")],
                label="Dropdown",
                nav=True,
            ),
        ]
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
        fig = px.scatter(
            df,
            x=dim_h,
            y=dim_v,
            color=color,
            symbol=symbol,
            hover_data=["name", "sex"],
            labels={
                dim_h: COLUMN_LABELS.get(dim_h, dim_h),
                dim_v: COLUMN_LABELS.get(dim_v, dim_v),
            },
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

            def build_select_div(
                field_id: str,
                label: str,
                initial_value: str,
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

            row_select = html.Form(
                [
                    build_select_div("stats-scatter-select-dim-h", "horizontal axis", "sex"),
                    build_select_div("stats-scatter-select-dim-v", "vertical axis", "chrx_het_hom"),
                    build_select_div("stats-scatter-select-dim-color", "color", "sex"),
                    build_select_div("stats-scatter-select-dim-symbol", "shape", "sex"),
                    html.Div(
                        [
                            dbc.Select(
                                options=[{"label": str(i), "value": str(i)} for i in range(1, 11)],
                                value="2",
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
        elif pathname == "/page-2":
            return html.P("Oh cool, this is page 2!")
        # If the user tries to reach a different page, return a 404 message
        return html.Div(
            [
                html.H1("404: Not found", className="text-danger"),
                html.Hr(),
                html.P(f"The pathname {pathname} was not recognised..."),
            ],
            className="p-3 bg-light rounded-3",
        )

    app.run_server()
