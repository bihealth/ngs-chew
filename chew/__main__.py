import argparse
import sys

from logzero import logger

from chew import __version__

from . import compare, fingerprint, plot_aab, plot_compare, plot_var_het, stats


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action="version", version="%(prog)s {}".format(__version__))

    subparser = parser.add_subparsers(dest="command")

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
    parser_fingerprint.add_argument("--input-bam", required=True, help="Path to input .bam file")
    parser_fingerprint.add_argument("--genome-release", required=False, help="Force genome release")
    parser_fingerprint.add_argument(
        "--max-sites", default=0, type=int, help="Optional, maximal number of sites to consider"
    )
    parser_fingerprint.add_argument(
        "--write-vcf",
        default=False,
        action="store_true",
        help="Write out the VCF file created during fingerprinting",
    )

    parser_compare = subparser.add_parser("compare")
    parser_compare.add_argument(
        "--output-prefix", default="chew-comparison", help="Path to comparison file."
    )
    parser_compare.add_argument(
        "fingerprints", nargs="+", help="Path(s) to .fingerprint.npz files to compare."
    )
    parser_compare.add_argument("--min-mask-ones", type=int, help="Minimal number of ones in mask.")
    parser_compare.add_argument("--max-mask-ones", type=int, help="Maximal number of ones in mask.")

    parser_stats = subparser.add_parser("stats")
    parser_stats.add_argument(
        "fingerprints", nargs="+", help="Path(s) to .fingerprint.npz to show stats for."
    )
    parser_stats.add_argument("--output", default="chew-stats.txt", help="Path to stats file.")

    parser_plot_compare = subparser.add_parser("plot_compare")
    parser_plot_compare.add_argument("stats", help="Path to `-stats.txt` file.")
    parser_plot_compare.add_argument("out_html", help="Path to output HTML file.")
    parser_plot_compare.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )

    parser_plot_aab = subparser.add_parser("plot_aab")
    parser_plot_aab.add_argument("out_html", help="Path to output HTML file.")
    parser_plot_aab.add_argument("vcf", nargs="+", help="Path(s) to input VCF files.")
    parser_plot_aab.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )
    parser_plot_aab.add_argument("--aab-cache", help="Path to AAB cache JSON file.")

    parser_plot_var_het = subparser.add_parser("plot_var_het")
    parser_plot_var_het.add_argument("out_html", help="Path to output HTML file.")
    parser_plot_var_het.add_argument("vcf", nargs="+", help="Path(s) to input VCF files.")
    parser_plot_var_het.add_argument(
        "--title", default="NGS Chew Comparison Plot", help="Title to use for the output HTML file."
    )
    parser_plot_var_het.add_argument("--var_het_cache", help="Path to AAB cache JSON file.")

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


if __name__ == "__main__":
    sys.exit(main())
