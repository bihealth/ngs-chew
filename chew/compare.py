import math

import attr
import numpy as np
from logzero import logger
from tqdm import tqdm


#re-work
@attr.s(frozen=True, auto_attribs=True)
class Fingerprint:
    """Store information from a fingerprint file."""
    DTYPE: type
    #: Genome release
    genome_release: str
    #: Name of the sample
    sample_name: str
    #: Genotype:  array for fingeprint
    data: np.array
    #: Allelic fraction array

    def maxval_to_frac(self,val):
        bits = t = np.iinfo(self.DTYPE).bits -3
        return float(2**(bits)-1 & val) / float(2**(bits)-1)

    def frac_to_maxval(self,frac:float):
        bits = t = np.iinfo(self.DTYPE).bits -3
        return round(frac * (2**bits-1))

    def get_bit(self,v,k:int):
        t = np.iinfo(self.DTYPE).bits
        if 0 <= k < t:
            return type(v)((v & 2**(t-k-1)) >> (t-k-1))
        else:
            raise ValueError("k must be 0 <= k < t")

    def set_bit(self,v,k:int,x=True):
        t = np.iinfo(self.DTYPE).bits
        if 0 <= k < t:
            return type(v)(2**(t-k-1) | v) if x else type(v)(~2**(t-k-1) & v)
        else:
            raise ValueError("k must be 0 <= k < t")

    def __init__(self,sample_name:str,genome_release:str,mask,is_alt,hom_alt,al_frac,dtype=np.uint16):
        self.DTYPE = dtype
        if not (len(mask) == len(is_alt) == len(hom_alt) == len(al_frac)):
            raise ValueError("All given lists or vectors must have the same length.")
        if al_frac.min() < 0 or al_frac.max() > 1.0:
            raise ValueError("Allelic fractions x must satisfy 0 <= x <= 1.")
        self.data = np.vectorize(self.frac_to_maxval,otypes=[self.DTYPE])(al_frac)
        self.data = np.vectorize(self.set_bit,otypes=[self.DTYPE])(self.data,2,hom_alt)
        self.data = np.vectorize(self.set_bit,otypes=[self.DTYPE])(self.data,1,is_alt)
        self.data = np.vectorize(self.set_bit,otypes=[self.DTYPE])(self.data,0,mask)
        self.genome_release = genome_release
        self.sample_name = sample_name

    def __getitem__(self, key:int):
        if key == 3:
            return np.vectorize(self.maxval_to_frac)(self.data).astype(float)
        if key >= 0:
            return np.vectorize(self.get_bit)(self.data,key).astype(bool)
        raise ValueError("Index i must satisfy 0 <= i <= %d."%3)


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
