import argparse
import sys

from logzero import logger

from . import compare, fingerprint, stats, plot_compare, plot_aab, plot_var_het
from chew import __version__


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))

    subparser = parser.add_subparsers(dest="command")

    # -------------------------------------------------------------------------
    #  parser_fingerprint
    parser_fingerprint = subparser.add_parser("fingerprint")
    parser_fingerprint.add_argument("--min-coverage", default=5, help="Minimal required coverage")
    parser_fingerprint.add_argument(
        "--reference", required=True, help="Path to reference FASTA file."
    )
    parser_fingerprint.add_argument(
        "--output-fingerprint",
        required=True,
        help="Path to output .npz file (extension .npz is added automatically if necessary).",
    )
    parser_fingerprint.add_argument(
        "--max-sites", default=0, type=int, help="Optional, maximal number of sites to consider"
    )
    parser_fingerprint.add_argument("--genome-release", required=False, help="Force genome release")

    subparser_fingerprints = parser_fingerprint.add_subparsers(dest="specific")

    # --- subparser fastq --- #
    parser_fastq = subparser_fingerprints.add_parser("fastq")
    parser_fastq.add_argument('--input', nargs='+', help="Path to input files containing reads in .fastq or .fastq.gz format.")
    parser_fastq.add_argument('--sample', help="name of the sample.")
    parser_fastq.add_argument('--k', default=15,type=int, help="h-mers are counted, where h = 2*k+1. Should range between 10 for 21-mers and 30 for 61-mers.")
    parser_fastq.add_argument('--bf_size' , default="12G",type=str, help="Size of bloom filter used by Jellyfish. Should be equivalent to the total number of distinct k-mers found in the data set. Use 100G for full homo sapiens WGS data set.")
    parser_fastq.add_argument('--depth' , default=300,type=int, help="sequencing depth. Used to find the expected number of hmers that occur more than once.")
    parser_fastq.add_argument('--error_rate' , default=0.01,type=float, help="error rate in the given ngs-reads")

        # --- subparser bam --- #
    parser_bam = subparser_fingerprints.add_parser("bam")
    parser_bam.add_argument('--input', help="Path to input alignment file in .bam format.")
    parser_bam.add_argument(
        "--write-vcf",
        default=False,
        action="store_true",
        help="Write out the VCF file created during fingerprinting if alignment file is given.",
    )
    # --- subparser vcf --- #
    parser_vcf = subparser_fingerprints.add_parser("vcf")
    parser_vcf.add_argument('--input', help="Path to input variant calling file in .vcf format.")

    # -------------------------------------------------------------------------
    #  parser_compare
    parser_compare = subparser.add_parser("compare")
    parser_compare.add_argument(
        "--output-prefix", default="chew-comparison", help="Path to comparison file."
    )
    parser_compare.add_argument(
        "fingerprints", nargs="+", help="Path(s) to .fingerprint.npz files to compare."
    )
    parser_compare.add_argument("--min-mask-ones", type=int, help="Minimal number of ones in mask.")
    parser_compare.add_argument("--max-mask-ones", type=int, help="Maximal number of ones in mask.")

    # -------------------------------------------------------------------------
    #  parser_stats
    parser_stats = subparser.add_parser("stats")
    parser_stats.add_argument(
        "fingerprints", nargs="+", help="Path(s) to .fingerprint.npz to show stats for."
    )
    parser_stats.add_argument("--output", default="chew-stats.txt", help="Path to stats file.")

    # -------------------------------------------------------------------------
    #  parser_plot_compare
    parser_plot_compare = subparser.add_parser("plot_compare")
    parser_plot_compare.add_argument("stats", help="Path to `-stats.txt` file.")
    parser_plot_compare.add_argument("out_html", help="Path to output HTML file.")
    parser_plot_compare.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )

    # -------------------------------------------------------------------------
    #  parser_plot_aab
    parser_plot_aab = subparser.add_parser("plot_aab")
    parser_plot_aab.add_argument("out_html", help="Path to output HTML file.")
    #parser_plot_aab.add_argument("vcf", nargs="+", help="Path(s) to input VCF files.")
    parser_plot_aab.add_argument("fps", nargs="+", help="Path(s) to input fingerprint (.npz) files.")
    parser_plot_aab.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )
    parser_plot_aab.add_argument("--aab-cache", help="Path to AAB cache JSON file.")

    # -------------------------------------------------------------------------
    #  parser_plot_var_het
    parser_plot_var_het = subparser.add_parser("plot_var_het")
    parser_plot_var_het.add_argument("out_html", help="Path to output HTML file.")
    parser_plot_var_het.add_argument("fps", nargs="+", help="Path(s) to input fingerprint (.npz) files.")
    parser_plot_var_het.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )
    parser_plot_var_het.add_argument("--var_het_cache", help="Path to AAB cache JSON file.")
    # -------------------------------------------------------------------------

    args = parser.parse_args(argv)
    logger.info("Options: %s" % vars(args))
    if args.command == "fingerprint":
        return fingerprint.run(args)
    elif args.command == "compare":
        return compare.run(args)
    elif args.command == "stats":
        return stats.run(args)
    elif args.command == "plot_compare":
        return plot_compare.run(args)
    elif args.command == "plot_aab":
        return plot_aab.run(args)
    elif args.command == "plot_var_het":
        return plot_var_het.run(args)
    else:
        parser.error("No command given!")


if __name__ == "__main__":  # pragma: nocover
    sys.exit(main())
