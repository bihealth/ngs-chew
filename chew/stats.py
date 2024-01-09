import typing

import attrs
from logzero import logger
import numpy as np
from tqdm import tqdm

from chew.common import CHROM_LENS_GRCH37, CHROM_LENS_GRCH38


@attrs.frozen
class Config:
    verbosity: int
    fingerprints: typing.List[str]
    output: str


@attrs.frozen
class Header:
    magic_string: str
    format_version: str
    chew_version: str
    release: str
    sample: str
    fields: typing.Tuple[str, ...]


def load_fingerprint_all(path):
    return np.load(path)


def extract_header(container) -> Header:
    magic_string = container["header"][0]
    version = container["header"][1]
    if magic_string != "ngs_chew_fingerprint":
        raise Exception(f"Wrong file format {magic_string}")
    if version != "4":
        raise Exception(f"Wrong version {version}, must be '4'")
    if len(container["header"]) < 5:
        raise Exception("Expected at least 5 fields in header")

    return Header(
        magic_string=container["header"][0],
        format_version=container["header"][1],
        chew_version=container["header"][2],
        release=container["header"][3],
        sample=container["header"][4],
        fields=tuple(container["header"][5].split(",")),
    )


def parse_samtools_idxstats(samtools_idxstats: str) -> typing.Tuple[float, float]:
    chrom_reads = {}
    for line in samtools_idxstats.splitlines():
        arr = line.split("\t")
        chrom = arr[0]
        reads = int(arr[2])
        if chrom in CHROM_LENS_GRCH37 or chrom in CHROM_LENS_GRCH38:
            chrom_reads[chrom] = reads
    total_count = sum(chrom_reads.values())
    chrx_count = chrom_reads.get("X", chrom_reads.get("chrX", 0))
    chry_count = chrom_reads.get("Y", chrom_reads.get("chrY", 0))
    if not total_count:
        return 0.0, 0.0
    else:
        return chrx_count / total_count, chry_count / total_count


def compute_chrx_het_hom(container):
    chrx_fingerprint = container["chrx_fingerprint"]
    chrx_mask = chrx_fingerprint[0]
    chrx_is_alt = chrx_fingerprint[1]
    chrx_hom_alt = chrx_fingerprint[2]
    num_homs = np.count_nonzero(chrx_hom_alt & chrx_mask) + np.count_nonzero(
        ~chrx_is_alt & chrx_mask
    )
    if num_homs > 0:
        return np.count_nonzero(chrx_is_alt & chrx_mask) / num_homs
    else:
        return None


def compute_autosomal_aafs(container):
    autosomal_fingerprint = container["autosomal_fingerprint"]
    autosomal_mask = autosomal_fingerprint[0]
    autosomal_is_alt = autosomal_fingerprint[1]
    autosomal_hom_alt = autosomal_fingerprint[2]
    autosomal_aafs = container["autosomal_aafs"]
    is_het = autosomal_mask & autosomal_is_alt & ~autosomal_hom_alt
    sqrt_var_het = autosomal_aafs[is_het] - 0.5
    return np.sum(sqrt_var_het * sqrt_var_het) / sqrt_var_het.shape[0]


@attrs.frozen
class SampleStats:
    sample_name: str
    release: str
    hets: int
    hom_alts: int
    hom_refs: int
    mask_ones: int
    var_het: typing.Optional[float]
    chrx_het_hom: typing.Optional[float]
    chrx_frac: typing.Optional[float]
    chry_frac: typing.Optional[float]


def compute_sample_stats(container) -> SampleStats:
    header = extract_header(container)

    autosomal_fingerprint = container["autosomal_fingerprint"]
    autosomal_mask = autosomal_fingerprint[0]
    autosomal_is_alt = autosomal_fingerprint[1]
    autosomal_hom_alt = autosomal_fingerprint[2]

    if "autosomal_aafs" in header.fields:
        var_het = compute_autosomal_aafs(container)
    else:
        var_het = None

    if "chrx_aafs" in header.fields:
        chrx_het_hom = compute_chrx_het_hom(container)
    else:
        chrx_het_hom = None

    if "samtools_idxstats" in header.fields:
        chrx_frac, chry_frac = parse_samtools_idxstats(str(container["samtools_idxstats"]))
    else:
        chrx_frac = None
        chry_frac = None

    return SampleStats(
        sample_name=header.sample,
        release=header.release,
        hets=np.count_nonzero(autosomal_is_alt & autosomal_mask),
        hom_alts=np.count_nonzero(autosomal_hom_alt & autosomal_mask),
        hom_refs=np.count_nonzero(~autosomal_is_alt & autosomal_mask),
        mask_ones=np.count_nonzero(autosomal_mask),
        var_het=var_het,
        chrx_het_hom=chrx_het_hom,
        chrx_frac=chrx_frac,
        chry_frac=chry_frac,
    )


def run(config: Config):
    logger.info("Writing statistics file...")
    with open(config.output, "wt") as outputf:
        col_headers = [
            "sample_name",
            "hets",
            "hom_alts",
            "hom_refs",
            "mask",
            "var_het",
            "chrx_het_hom",
            "chrx_frac",
            "chry_frac",
        ]

        print("\t".join(col_headers), file=outputf)
        for container in map(load_fingerprint_all, tqdm(config.fingerprints)):
            stats = compute_sample_stats(container)

            row = []
            for key in col_headers:
                if getattr(stats, key, None) is None:
                    row.append("-")
                else:
                    row.append(getattr(stats, key))

            print("\t".join(map(str, row)), file=outputf)
