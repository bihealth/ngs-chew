import argparse
import sys

from logzero import logger

from . import compare, fingerprint, stats
from chew import __version__


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
    parser_fingerprint.add_argument("--input-bam", required=True, help="Path to input .bam file.")
    parser_fingerprint.add_argument("--genome-release", required=False, help="Force genome release")
    parser_fingerprint.add_argument(
        "--max-sites", default=0, type=int, help="Optional, maximal number of sites to consider"
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

    args = parser.parse_args(argv)
    logger.info("Options: %s" % vars(args))
    if args.command == "fingerprint":
        return fingerprint.run(args)
    elif args.command == "compare":
        return compare.run(args)
    elif args.command == "stats":
        return stats.run(args)
    else:
        parser.error("No command given!")


if __name__ == "__main__":
    sys.exit(main())
