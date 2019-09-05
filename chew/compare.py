import math

import numpy as np
from logzero import logger
from tqdm import tqdm


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
    rel = (het_ij - 2 * n_ibs0) / (0.5 * math.sqrt(het_i * het_j))
    return n_ibs0, rel, np.count_nonzero(mask)


def run(args):
    logger.info("Loading fingerprints...")
    fps = {
        name: fingerprint for name, fingerprint in map(load_fingerprint, tqdm(args.fingerprints))
    }
    logger.info("Loaded %d fingerprints", len(fps))
    if args.min_mask_ones or args.max_mask_ones:
        fps = {
            name: fp
            for name, fp in fps.items()
            if (not args.min_mask_ones or np.count_nonzero(fp[0]) >= args.min_mask_ones)
            and (not args.max_mask_ones or np.count_nonzero(fp[0]) <= args.max_mask_ones)
        }
    logger.info("Filtered to %d fingerprints", len(fps))
    header = ["name"] + list(fps)
    rows_n_ibs0 = []
    rows_rel = []
    rows_mask = []
    for lhs in tqdm(fps):
        row_n_ibs0 = [lhs]
        row_rel = [lhs]
        row_mask = [lhs]
        for rhs in fps:
            n_ibs0, rel, mask = relatedness(fps[lhs], fps[rhs])
            row_n_ibs0.append(n_ibs0)
            row_rel.append(rel)
            row_mask.append(mask)
        rows_n_ibs0.append(row_n_ibs0)
        rows_rel.append(row_rel)
        rows_mask.append(row_mask)

    with open("%s.relatedness.txt" % args.output_prefix, "wt") as outputf:
        for row in [header] + rows_rel:
            print("\t".join(map(str, row)), file=outputf)
    with open("%s.ibs0.txt" % args.output_prefix, "wt") as outputf:
        for row in [header] + rows_n_ibs0:
            print("\t".join(map(str, row)), file=outputf)
    with open("%s.mask.txt" % args.output_prefix, "wt") as outputf:
        for row in [header] + rows_mask:
            print("\t".join(map(str, row)), file=outputf)
