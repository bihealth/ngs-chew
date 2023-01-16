import contextlib
import gzip
import os
import shlex
import subprocess
import tempfile
import typing
import warnings

import attrs
from chew.stats import extract_header, load_fingerprint_all
from logzero import logger
import numpy as np
import numpy.typing
import pysam
import vcfpy

import chew
from chew.common import CHROM_LENS_GRCH37, CHROM_LENS_GRCH38, CHROM_LENS_HG19

#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Template for creating ``bcftools mpileup`` call.
BCFTOOLS_PILEUP = (
    r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"
)

#: Template for creating ``bcftools call`` call.
BCFTOOLS_CALL = r"bcftools call -c -Oz -o %(calls)s"

#: Template for creating ``samtools idxstats`` call.
SAMTOOLS_IDXSTATS = r"samtools idxstats %(input_bam)s"


class InvalidInputWarning(UserWarning):
    """Raised when using input warning."""


@attrs.frozen
class Config:
    #: Verbosity level
    verbosity: int
    #: Input file
    input_npz: str
    #: Output file
    output: str


def create_vcf_header(sample: str, release: str) -> vcfpy.Header:
    if release == "GRCh37":
        chrom_lens = CHROM_LENS_GRCH38
    elif release == "GRCh38":
        chrom_lens = CHROM_LENS_GRCH38
    else:
        raise RuntimeError(f"Invalid release {release}")
    header = vcfpy.Header(samples=vcfpy.SamplesInfos([sample]))
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    for name, length in chrom_lens.items():
        header.add_contig_line({"ID": name, "length": length})
    header.add_filter_line({"ID": "PASS", "Description": "All filters passed"})
    header.add_format_line({"ID": "GT", "Number": 1, "Type": "Integer", "Description": "Genotype"})
    header.add_format_line({"ID": "AD", "Number": "R", "Type": "Integer", "Description": "Allelic depths for the ref and alt alleles in the order listed"})
    header.add_format_line({"ID": "DP", "Number": 1, "Type": "Integer", "Description": "pproximate read depth"})
    return header


def run(config: Config):
    logger.info("Reading fingerprint file...")
    container = load_fingerprint_all(config.input_npz)
    fp_header = extract_header(container)
    vcf_header = create_vcf_header(fp_header.sample, fp_header.release)
    with vcfpy.Writer.from_path(config.output, vcf_header) as writer:
        pass

    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Prepare VCF files...")

        logger.info("Run 'bcftools roh' ...")

    logger.info("All done. Have a nice day!")
