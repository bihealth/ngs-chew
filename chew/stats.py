import typing

import attrs
from logzero import logger
import numpy as np
from numpy import ma
from tqdm import tqdm

from chew.compare import load_fingerprint_with_aafs


@attrs.frozen
class Config:
    verbosity: int
    fingerprints: typing.List[str]
    output: str


def run(config: Config):
    logger.info("Writing statistics file...")
    with open(config.output, "wt") as outputf:
        header = ["sample", "hets", "hom_alts", "hom_refs", "mask", "var_het"]
        print("\t".join(header), file=outputf)
        for name, fp, aafs in map(load_fingerprint_with_aafs, tqdm(config.fingerprints)):
            mask = fp[0]
            is_alt = fp[1]
            hom_alt = fp[2]
            if aafs is not None:
                is_het = mask & is_alt & ~hom_alt
                sqrt_var_het = aafs[is_het] - 0.5
                var_het = np.sum(sqrt_var_het * sqrt_var_het) / sqrt_var_het.shape[0]
            else:
                var_het = None
            row = [
                name,
                np.count_nonzero(is_alt & mask),
                np.count_nonzero(hom_alt & mask),
                np.count_nonzero(~is_alt & mask),
                np.count_nonzero(mask),
                var_het or "-",
            ]
            print("\t".join(map(str, row)), file=outputf)
