import typing

import attrs
from logzero import logger
import numpy as np
from tqdm import tqdm


@attrs.frozen
class Config:
    verbosity: int
    fingerprints: typing.List[str]
    output: str


@attrs.frozen
class Header:
    magic_string: str
    version: str
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
    if version != "3":
        raise Exception(f"Wrong version {version}, must be '3'")
    if len(container["header"]) < 5:
        raise Exception("Expected at least 5 fields in header")

    return Header(
        magic_string=container["header"][0],
        version=container["header"][1],
        release=container["header"][2],
        sample=container["header"][3],
        fields=tuple(container["header"][4].split(",")),
    )


def run(config: Config):
    logger.info("Writing statistics file...")
    with open(config.output, "wt") as outputf:
        header_lines = ["sample", "hets", "hom_alts", "hom_refs", "mask", "var_het", "chrx_het_hom"]

        print("\t".join(header_lines), file=outputf)
        for container in map(load_fingerprint_all, tqdm(config.fingerprints)):
            header = extract_header(container)

            autosomal_fingerprint = container["autosomal_fingerprint"]
            autosomal_mask = autosomal_fingerprint[0]
            autosomal_is_alt = autosomal_fingerprint[1]
            autosomal_hom_alt = autosomal_fingerprint[2]

            if "autosomal_aafs" in header.fields:
                autosomal_aafs = container["autosomal_aafs"]
                is_het = autosomal_mask & autosomal_is_alt & ~autosomal_hom_alt
                sqrt_var_het = autosomal_aafs[is_het] - 0.5
                var_het = np.sum(sqrt_var_het * sqrt_var_het) / sqrt_var_het.shape[0]
            else:
                var_het = None

            if "chrx_aafs" in header.fields:
                chrx_fingerprint = container["chrx_fingerprint"]
                chrx_mask = chrx_fingerprint[0]
                chrx_is_alt = chrx_fingerprint[1]
                chrx_hom_alt = chrx_fingerprint[2]
                num_homs = np.count_nonzero(chrx_hom_alt & chrx_mask) + np.count_nonzero(
                    ~chrx_is_alt & chrx_mask
                )
                if num_homs > 0:
                    chrx_het_hom = np.count_nonzero(chrx_is_alt & chrx_mask) / num_homs
                else:
                    chrx_het_hom = None
            else:
                chrx_het_hom = None

            row = [
                header.sample,
                np.count_nonzero(autosomal_is_alt & autosomal_mask),
                np.count_nonzero(autosomal_hom_alt & autosomal_mask),
                np.count_nonzero(~autosomal_is_alt & autosomal_mask),
                np.count_nonzero(autosomal_mask),
                var_het or "-",
                chrx_het_hom or "-",
            ]
            print("\t".join(map(str, row)), file=outputf)
