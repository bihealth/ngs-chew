import os
import shlex
import subprocess

import attrs
from logzero import logger
import numpy.typing as npt
import vcfpy

from chew.common import CHROM_LENS_GRCH37, CHROM_LENS_GRCH38, load_sites

#: Use this value for genotype likelyhood (FORMAT/PL)
BCFTOOLS_ROH_PL = 30

#: Use this value for assumed allele frequency
BCFTOOLS_ROH_AF = 0.05

#: Template for creating ``bcftools roh`` call.
BCFTOOLS_ROH = r"bcftools roh -G %(pl)f --AF-dflt %(af)f %(input_vcf)s -O r -o %(output_txt)s"


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
        chrom_lens = CHROM_LENS_GRCH37
    elif release == "GRCh38":
        chrom_lens = CHROM_LENS_GRCH38
    else:
        raise RuntimeError(f"Invalid release {release}")
    header = vcfpy.Header(samples=vcfpy.SamplesInfos([sample]))
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    header.add_line(vcfpy.HeaderLine("reference", release))
    for name, length in chrom_lens.items():
        header.add_contig_line({"ID": name, "length": length})
    header.add_filter_line({"ID": "PASS", "Description": "All filters passed"})
    header.add_format_line({"ID": "GT", "Number": 1, "Type": "String", "Description": "Genotype"})
    return header


def write_vcf(tmpdir: str, sample: str, release: str, autosomal_fingerprint) -> str:
    logger.info("Constructing VCF header...")
    vcf_header = create_vcf_header(sample, release)
    sites = load_sites(release)
    autosomal_mask = autosomal_fingerprint[0]
    autosomal_is_alt = autosomal_fingerprint[1]
    autosomal_hom_alt = autosomal_fingerprint[2]
    path_vcf = os.path.join(tmpdir, "gts.vcf")
    if release == "GRCh37":
        prefix = ""
    else:
        prefix = "chr"
    with vcfpy.Writer.from_path(path_vcf, vcf_header) as writer:
        for i, site in enumerate(sorted(sites)):
            if not autosomal_mask[i]:
                gt = "./."
            elif autosomal_hom_alt[i]:
                gt = "1/1"
            elif autosomal_is_alt[i]:
                gt = "0/1"
            else:
                gt = "0/0"
            record = vcfpy.Record(
                CHROM=f"{prefix}{site.chrom}",
                POS=site.pos,
                ID=[],
                REF="N",
                ALT=[vcfpy.Substitution(vcfpy.SNV, "A")],
                QUAL=None,
                FILTER=[],
                INFO={},
                FORMAT=["GT"],
                calls=[
                    vcfpy.Call(
                        sample=sample,
                        data={
                            "GT": gt,
                        },
                    )
                ],
            )
            writer.write_record(record)
    return path_vcf


def run_roh(tmpdir: str, path_vcf: str) -> str:
    logger.info("Run 'bcftools roh' ...")
    path_txt = os.path.join(tmpdir, "roh.txt")
    cmd_roh = BCFTOOLS_ROH % {
        "pl": BCFTOOLS_ROH_PL,
        "af": BCFTOOLS_ROH_AF,
        "input_vcf": path_vcf,
        "output_txt": path_txt,
    }
    logger.info("  roh: %s", " ".join(shlex.split(cmd_roh)))
    p = subprocess.Popen(shlex.split(cmd_roh))
    p.wait()
    return path_txt


@attrs.frozen
class RohTxtContents:
    chroms: npt.ArrayLike
    starts: npt.ArrayLike
    ends: npt.ArrayLike
    quals: npt.ArrayLike


def parse_roh_txt(path_txt: str) -> RohTxtContents:
    chroms = []
    starts = []
    ends = []
    quals = []
    with open(path_txt, "rt") as inputf:
        for line in inputf:
            if not line.startswith("RG"):
                continue
            arr = line.split("\t")
            chroms.append(arr[2])
            starts.append(int(arr[3]))
            ends.append(int(arr[4]))
            quals.append(float(arr[7]))
    return RohTxtContents(chroms=chroms, starts=starts, ends=ends, quals=quals)
