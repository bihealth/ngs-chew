"""Commonly used code"""

import enum
import gzip
import os
import typing

import attrs
from logzero import logger

#: Chromosome lengths in GRCh37.
CHROM_LENS_GRCH37 = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
}

#: Chromosome lengths in hg19.
CHROM_LENS_HG19 = {f"chr{name}": length for name, length in CHROM_LENS_GRCH37.items()}

#: Chromosome lengths in GRCh38.
CHROM_LENS_GRCH38 = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}


class Sex(enum.Enum):
    UNKNOWN = "unknown"
    MALE = "male"
    FEMALE = "female"


PED_SEX_MAP = {
    "0": Sex.UNKNOWN,
    "1": Sex.MALE,
    "2": Sex.FEMALE,
}


class DiseaseState(enum.Enum):
    UNKNOWN = "unknown"
    UNAFFECTED = "affected"
    AFFECTED = "unaffected"


PED_DISEASE_MAP = {
    "0": DiseaseState.UNKNOWN,
    "1": DiseaseState.UNAFFECTED,
    "2": DiseaseState.AFFECTED,
}


@attrs.frozen
class PedigreeMember:
    family_name: str
    name: str
    father: str
    mother: str
    sex: Sex
    disease_state: DiseaseState


def pedigree_member_from_tsv(arr: typing.List[str]) -> PedigreeMember:
    if len(arr) < 6:
        raise Exception("TSV array must have at least 6 fields")
    return PedigreeMember(
        family_name=arr[0],
        name=arr[1],
        father=arr[2],
        mother=arr[3],
        sex=PED_SEX_MAP.get(arr[4], Sex.UNKNOWN),
        disease_state=PED_DISEASE_MAP.get(arr[5], DiseaseState.UNKNOWN),
    )


@attrs.frozen(order=True)
class Site:
    chrom: str
    pos: int
    ref: str
    alt: str


def load_sites(genome_release: str) -> typing.List[Site]:
    logger.info("Loading sites .bed.gz for %s", genome_release)
    path_gz = os.path.join(os.path.dirname(__file__), "data", f"{genome_release}_sites.bed.gz")
    result = []
    with gzip.open(path_gz, "rt") as inputf:
        for line in inputf:
            arr = line.split("\t")
            # We construct sites with arbitrary REF/ALT as bcftools roh does not
            # interpret them but only GT.
            result.append(
                Site(
                    chrom=arr[0],
                    pos=int(arr[1]) + 1,
                    ref="N",
                    alt="A",
                )
            )
    logger.info("  finished reading %d sites", len(result))
    logger.info("  first = %s, last = %s", result[0], result[-1])
    return result
