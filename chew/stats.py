from logzero import logger
import numpy as np
import numpy.ma as ma
from tqdm import tqdm

from .compare import load_fingerprint


def run(args):
    logger.info("Writing statistics file...")
    with open(args.output, "wt") as outputf:
        print("\t".join(["sample", "hets", "hom_alts", "hom_refs", "mask", "mean_af"]), file=outputf)
        for name, fp, allelic_fractions in map(load_fingerprint, tqdm(args.fingerprints)):
            mask = fp[0]
            is_alt = fp[1]
            hom_alt = fp[2]
            row = [
                name,
                np.count_nonzero(is_alt & mask),
                np.count_nonzero(hom_alt & mask),
                np.count_nonzero(~is_alt & mask),
                np.count_nonzero(mask),
                ma.masked_array(allelic_fractions, mask=mask).mean()
            ]
            print("\t".join(map(str, row)), file=outputf)
