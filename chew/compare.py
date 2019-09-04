import math

import numpy as np

from logzero import logger


def load_fingerprint(path):
    nparr = np.load(path)
    return nparr["header"][3], nparr["fingerprint"]


def relatedness(lhs, rhs):
    # Obtain shortcuts...
    lhs_mask = lhs[0]
    lhs_is_alt = lhs[1]
    lhs_hom_alt = rhs[2]
    rhs_mask = rhs[0]
    rhs_is_alt = rhs[1]
    rhs_hom_alt = rhs[2]
    mask = lhs_mask & rhs_mask
    # Compute bit array
    lhs_ref = ~lhs_is_alt & mask
    rhs_ref = ~rhs_is_alt & mask
    arr_i = lhs_is_alt & ~lhs_hom_alt & mask
    arr_j = rhs_is_alt & ~rhs_hom_alt & mask
    arr_ij = arr_i & arr_j & mask
    # Compute partial terms
    het_i = np.count_nonzero(arr_i)
    het_j = np.count_nonzero(arr_j)
    het_ij = np.count_nonzero(arr_ij)
    n_ibs0 = np.count_nonzero(((lhs_ref & rhs_hom_alt) | (rhs_ref & lhs_hom_alt)) & mask)
    # Compute relateness
    return (het_ij - 2 * n_ibs0) / (0.5 * math.sqrt(het_i * het_j))


def run(args):
    logger.info("Loading fingerprints...")
    fps = {name: fingerprint for name, fingerprint in map(load_fingerprint, args.fingerprints)}
    logger.info("Loaded %d fingerprints", len(fps))
    print("\t" + "\t".join(fps.keys()))
    for lhs in fps.keys():
        row = [lhs] + [relatedness(fps[lhs], fps[rhs]) for rhs in fps.keys()]
        print("\t".join(map(str, row)))
