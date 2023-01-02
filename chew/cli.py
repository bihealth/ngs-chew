import sys
import typing

import click
from logzero import logger

from chew import __version__, compare, fingerprint, plot_compare, plot_var_het, stats


@click.group()
@click.version_option(__version__)
@click.option("--verbose/--no-verbose", help="Enable verbose mode.", default=False)
@click.pass_context
def cli(ctx: click.Context, verbose: bool):
    """Main entry point for CLI via click."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose


@cli.command("fingerprint", help="Compute fingerprint to numpy .npz files.")
@click.option("--min-coverage", type=int, default=5, help="Minimal required coverage.")
@click.option("--reference", required=True, help="Path to reference FASTA file.")
@click.option(
    "--output-fingerprint",
    required=True,
    help="Path to output .npz file (extension will be added automatically if necessary)",
)
@click.option(
    "--output-aafs/--no-output-aafs",
    default=False,
    help="Write alternate allele fractions to .npz file.",
)
@click.option("--input-bam", required=True, help="Path to input BAM file.")
@click.option("--genome-release", required=False, default="GRCh37", help="Genome release used.")
@click.option("--max-sites", type=int, default=0, help="Maximal number of sites to consider.")
@click.option("--write-vcf/--no-write-vcf", default=False, help="Enable writing of call VCF.")
@click.pass_context
def cli_fingerprint(
    ctx: click.Context,
    min_coverage: int,
    reference: str,
    output_fingerprint: bool,
    output_aafs: bool,
    input_bam: str,
    genome_release: str,
    max_sites: int,
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
        write_vcf=write_vcf,
    )
    fingerprint.run(config)


@cli.command("compare", help="Perform fingeprint comparison.")
@click.option(
    "--output-prefix", type=str, default="chew-comparison", help="Path to comparison file."
)
@click.option("--min-mask-ones", type=int, help="Minimal number of ones in mask.")
@click.option("--max-mask-ones", type=int, help="Maximal number of ones in mask.")
@click.argument("fingerprints", nargs=-1)
@click.pass_context
def cli_compare(
    ctx: click.Context,
    output_prefix: str,
    min_mask_ones: typing.Optional[int],
    max_mask_ones: typing.Optional[int],
    fingerprints: typing.List[str],
):
    config = compare.Config(
        verbosity=2 if ctx.obj["verbose"] else 1,
        output_prefix=output_prefix,
        min_mask_ones=min_mask_ones,
        max_mask_ones=max_mask_ones,
        fingerprints=fingerprints,
    )
    compare.run(config)


@cli.command("stats", help="Compute statistics from fingerprint .npz files.")
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


@cli.command(help="Plot result of 'ngs-chew compare'.")
@click.argument("stats")
@click.option(
    "--title", default="NGS Chew Comparison Plot", help="title to use for the output HTML file."
)
def plot_compare(
    ctx: click.Context,
    out_html: str,
):
    pass


@cli.command(help="Plot var(het) metric from .npz files.")
@click.argument("fingerprints", nargs=-1)
@click.option(
    "--title", default="NGS Chew var(het) Plot", help="title to use for the output HTML file."
)
def plot_var_het(
    ctx: click.Context,
    out_html: str,
    fingerprints: typing.List[str],
):
    pass


if __name__ == "__main__":
    sys.exit(main())
