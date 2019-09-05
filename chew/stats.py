from logzero import logger
import numpy as np
from tqdm import tqdm

from .compare import load_fingerprint


def run(args):
    logger.info("Writing statistics file...")
    with open(args.output, "wt") as outputf:
        print("\t".join(["sample", "hets", "hom_alts", "hom_refs", "mask"]), file=outputf)
        for name, fp in map(load_fingerprint, tqdm(args.fingerprints)):
            mask = fp[0]
            is_alt = fp[1]
            hom_alt = fp[2]
            row = [
                name,
                np.count_nonzero(is_alt & mask),
                np.count_nonzero(hom_alt & mask),
                np.count_nonzero(~is_alt & mask),
                np.count_nonzero(mask),
            ]
            print("\t".join(map(str, row)), file=outputf)
