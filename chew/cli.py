import typing

import click
from logzero import logger

from chew import __version__, compare, fingerprint, plot_compare, plot_var_het, stats

try:
    import dash  # noqa

    have_dash_installed = True
except ImportError:
    have_dash_installed = False


@click.group()
@click.version_option(__version__)
@click.option("--verbose/--no-verbose", help="Enable verbose mode.", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command("fingerprint", help="Compute fingerprint to numpy .npz files.")  # type: ignore[attr-defined]
@click.option("--min-coverage", type=int, default=5, help="Minimal required coverage.")
@click.option("--reference", required=True, help="Path to reference FASTA file.")
@click.option(
    "--output-fingerprint",
    required=True,
    help="Path to output .npz file (extension will be added automatically if necessary)",
)
@click.option(
    "--output-aafs/--no-output-aafs",
    default=True,
    help="Write alternate allele fractions to .npz file.",
)
@click.option("--input-bam", required=True, help="Path to input BAM file.")
@click.option("--genome-release", required=False, default=None, help="Genome release used.")
@click.option("--max-sites", type=int, default=0, help="Maximal number of sites to consider.")
@click.option("--write-vcf/--no-write-vcf", default=False, help="Enable writing of call VCF.")
@click.option(
    "--step-autosomal-snps/--no-step-autosomal-snps",
    default=True,
    help="Enable autosomal SNP step (default: yes)",
)
@click.option(
    "--step-chrx-snps/--no-step-chrx-snps",
    default=True,
    help="Enable chrX SNP step (default: yes)",
)
@click.option(
    "--step-samtools-idxstats/--no-step-samtools-idxstats",
    default=True,
    help="Enable samtools idxstats step (default: yes)",
)
@click.option(
    "--step-bcftools-roh/--no-step-bcftools-roh",
    default=True,
    help="Enable bcftools roh step (default: yes)",
)
@click.option("--write-vcf/--no-write-vcf", default=False, help="Enable writing of call VCF.")
@click.pass_context
def cli_fingerprint(
    ctx: click.Context,
    min_coverage: int,
    reference: str,
    output_fingerprint: str,
    output_aafs: bool,
    input_bam: str,
    genome_release: typing.Optional[str],
    max_sites: int,
    step_autosomal_snps: bool,
    step_chrx_snps: bool,
    step_samtools_idxstats: bool,
    step_bcftools_roh: bool,
    write_vcf: bool,
):
    config = fingerprint.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        min_coverage=min_coverage,
        reference=reference,
        output_fingerprint=output_fingerprint,
        output_aafs=output_aafs,
        input_bam=input_bam,
        genome_release=genome_release,
        max_sites=max_sites,
        step_autosomal_snps=step_autosomal_snps,
        step_chrx_snps=step_chrx_snps,
        step_samtools_idxstats=step_samtools_idxstats,
        step_bcftools_roh=step_bcftools_roh,
        write_vcf=write_vcf,
    )
    fingerprint.run(config)


@cli.command("compare", help="Perform fingeprint comparison.")  # type: ignore[attr-defined]
@click.option(
    "--output-prefix", type=str, default="chew-comparison", help="Path to comparison file."
)
@click.option("--min-mask-ones", type=int, help="Minimal number of ones in mask.")
@click.option("--max-mask-ones", type=int, help="Maximal number of ones in mask.")
@click.option("--by-path/--no-by-path", default=False, help="Use path as fingerprint name.")
@click.argument("fingerprints", nargs=-1)
@click.pass_context
def cli_compare(
    ctx: click.Context,
    output_prefix: str,
    min_mask_ones: typing.Optional[int],
    max_mask_ones: typing.Optional[int],
    fingerprints: typing.List[str],
    by_path: bool,
):
    config = compare.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        output_prefix=output_prefix,
        min_mask_ones=min_mask_ones,
        max_mask_ones=max_mask_ones,
        fingerprints=fingerprints,
        by_path=by_path,
    )
    compare.run(config)


@cli.command("stats", help="Compute statistics from fingerprint .npz files.")  # type: ignore[attr-defined]
@click.option("--output", default="chew-stats.txt", help="Path to output file.")
@click.argument("fingerprints", nargs=-1)
@click.pass_context
def cli_stats(
    ctx: click.Context,
    output: str,
    fingerprints: typing.List[str],
):
    config = stats.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        output=output,
        fingerprints=fingerprints,
    )
    stats.run(config)


@cli.command("plot-compare", help="Plot result of 'ngs-chew compare'.")  # type: ignore[attr-defined]
@click.option(
    "--title", default="NGS Chew Comparison Plot", help="title to use for the output HTML file."
)
@click.argument("compare-out")
@click.argument("out_html")
@click.pass_context
def cli_plot_compare(
    ctx: click.Context,
    title: str,
    compare_out: str,
    out_html: str,
):
    config = plot_compare.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        compare_out=compare_out,
        out_html=out_html,
        title=title,
    )
    plot_compare.run(config)


@cli.command("plot-var-het", help="Plot var(het) metric from .npz files.")  # type: ignore[attr-defined]
@click.option(
    "--title", default="NGS Chew var(het) Plot", help="title to use for the output HTML file."
)
@click.argument("stats_out")
@click.argument("out_html")
@click.pass_context
def cli_plot_var_het(
    ctx: click.Context,
    title: str,
    stats_out: str,
    out_html: str,
):
    config = plot_var_het.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        title=title,
        out_html=out_html,
        stats_out=stats_out,
    )
    plot_var_het.run(config)


if have_dash_installed:

    @cli.command("serve", help="Run report server")  # type: ignore[attr-defined]
    @click.option(
        "--annos-tsv",
        default=None,
        required=False,
        help="Optional TSV file with further annotations",
    )
    @click.option("--strip-suffix", default="", help="Suffix to strip from sample names")
    @click.argument("cohort_ped")
    @click.argument("fingerprints", nargs=-1)
    @click.pass_context
    def cli_serve(
        ctx: click.Context,
        annos_tsv: typing.Optional[str],
        strip_suffix: str,
        cohort_ped: str,
        fingerprints: typing.List[str],
    ):
        if not fingerprints:
            logger.warn("No fingerprints given!")
            return

        from chew import serve

        config = serve.Config(
            verbosity=2 if ctx.obj["verbose"] else 1,
            strip_suffix=strip_suffix,
            cohort_ped=cohort_ped,
            fingerprints=fingerprints,
            annos_tsv=annos_tsv,
        )
        serve.run(config)
