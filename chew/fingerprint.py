import contextlib
import gzip
import os
import shlex
import subprocess
import tempfile
import typing
import warnings

import attrs
from logzero import logger
import numpy as np
import pysam
import vcfpy

#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Template for creating ``bcftools mpileup`` call.
BCFTOOLS_PILEUP = (
    r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"
)

#: Template for creating ``bcftools call`` call.
BCFTOOLS_CALL = r"bcftools call -c -Oz -o %(calls)s"

#: Template for creating ``samtools idxstats`` call.
SAMTOOLS_IDXSTATS = r"samtools idxstats %(input_bam)s"

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


class InvalidInputWarning(UserWarning):
    """Raised when using input warning."""


@attrs.frozen
class Config:
    verbosity: int
    min_coverage: int
    reference: str
    output_fingerprint: str
    output_aafs: bool
    genome_release: typing.Optional[str]
    input_bam: str
    max_sites: int
    write_vcf: bool
    #: Whether to collect autosomal snp information.
    step_autosomal_snps: bool
    #: Whether to collect "samtools idxstats" output.  Only applicable when extracting information
    #: from BAM files.
    step_samtools_idxstats: bool


def analyze_bam_header(
    input_bam, genome_release: typing.Optional[str] = None
) -> typing.Tuple[str, str, str]:
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        if "RG" in samfile.header:
            sample_names = {rg["SM"] for rg in samfile.header["RG"] if "SM" in rg}
        else:
            sample_names = set()

        samfile_chrom_lens = {
            record["SN"]: record["LN"]
            for record in samfile.header.get("SQ")
            if "SN" in record
            and record["SN"] in CHROM_LENS_GRCH37
            or record["SN"] in CHROM_LENS_GRCH38
        }

    if len(sample_names) != 1:
        warnings.warn(
            InvalidInputWarning(
                "Could not determine sample name from BAM file's read groups. "
                "Inferring from file name."
            )
        )
        sample_name = os.path.basename(input_bam)[: -len(".bam")]
    else:
        sample_name = list(sample_names)[0]
        logger.info("Using sample name %s", sample_name)

    if genome_release:
        logger.info("Using genome release from command line: %s", genome_release)
    else:
        if samfile_chrom_lens == CHROM_LENS_GRCH37:
            genome_release = "GRCh37"
            prefix = ""
        elif samfile_chrom_lens == CHROM_LENS_HG19:
            genome_release = "GRCh37"
            prefix = "chr"
        elif samfile_chrom_lens == CHROM_LENS_GRCH38:
            genome_release = "GRCh38"
            prefix = "chr"
        else:
            raise Exception("")
        logger.info("Guessing genome release to be %s", genome_release)

    return sample_name, prefix, genome_release


def write_sites_bed(args, prefix, genome_release, tmp_dir):
    path_gz = os.path.join(os.path.dirname(__file__), "data", "%s_sites.bed.gz" % genome_release)
    path_bed = os.path.join(tmp_dir, "sites.bed")
    logger.info("Writing sites BED file to %s", path_bed)
    # TODO: handle "chr" prefix...
    with gzip.open(path_gz, "rt") as inputf:
        with open(path_bed, "wt") as outputf:
            for lineno, line in enumerate(inputf):
                if not args.max_sites or lineno < args.max_sites:
                    print("%s%s" % (prefix, line.strip()), file=outputf)
    logger.info("Wrote %s sites", "{:,}".format(lineno))
    return path_bed


def call_sites(config: Config, path_sites: str, tmp_dir: str):
    logger.info("Performing variant calling at sites")
    path_vcf = os.path.join(tmp_dir, "calls.vcf.gz")
    cmd_pileup = BCFTOOLS_PILEUP % {
        "sites": path_sites,
        "reference": config.reference,
        "input_bam": config.input_bam,
    }
    cmd_call = BCFTOOLS_CALL % {"calls": path_vcf}
    logger.info("  mpileup: %s", " ".join(shlex.split(cmd_pileup)))
    logger.info("  call:    %s", " ".join(shlex.split(cmd_call)))
    p_pileup = subprocess.Popen(shlex.split(cmd_pileup), stdout=subprocess.PIPE)
    p_call = subprocess.Popen(shlex.split(cmd_call), stdin=p_pileup.stdout)
    p_call.wait()
    p_pileup.wait()
    return path_vcf


def autosomal_snps_step_call(
    config: Config,
    chr_prefix: str,
    genome_release: str,
    path_calls: str,
    prefix_fingerprint: typing.Optional[str],
):
    logger.info("Reading sites BED...")
    path_gz = os.path.join(os.path.dirname(__file__), "data", "%s_sites.bed.gz" % genome_release)
    with gzip.open(path_gz, "rt") as inputf:
        sites = {}
        for line in inputf:
            arr = line.strip().split("\t")
            sites["%s%s:%s" % (chr_prefix, arr[0], int(arr[1]) + 1)] = (0, 0, float("nan"))
    logger.info("Converting VCF to fingerprint...")
    with vcfpy.Reader.from_path(path_calls) as vcf_reader:
        if prefix_fingerprint:
            logger.info("Writing VCF to %s", prefix_fingerprint + ".vcf.gz")
            out_vcf = vcfpy.Writer.from_path(prefix_fingerprint + ".vcf.gz", vcf_reader.header)
        else:
            logger.info("Not writing out VCF")
            out_vcf = contextlib.suppress()
        with out_vcf as vcf_writer:
            sample = vcf_reader.header.samples.names[0]
            for record in vcf_reader:
                if prefix_fingerprint:
                    vcf_writer.write_record(record)
                key = "%s%s:%s" % (chr_prefix, record.CHROM, record.POS)
                if key in sites:
                    call = record.call_for_sample[sample]
                    cov = sum(call.data["AD"])
                    if cov:
                        aab_at_site = min(
                            call.data["AD"][0] / cov,
                            1.0 - call.data["AD"][0] / cov,
                        )
                    else:
                        aab_at_site = float("nan")
                    sites[key] = (record.INFO["DP"], call.gt_type, aab_at_site)
    depths = [dp for dp, _, _ in sites.values()]
    genotypes = [gt for _, gt, _ in sites.values()]
    aafs = np.array([aab for _, _, aab in sites.values()], dtype="float32")
    fingerprint = np.array(
        [
            [dp > config.min_coverage for dp in depths],
            [gt != vcfpy.HOM_REF for gt in genotypes],
            [gt == vcfpy.HOM_ALT for gt in genotypes],
        ],
        dtype=bool,
    )
    return fingerprint, aafs


def autosomal_snps_step(
    config: Config,
    chr_prefix: str,
    genome_release: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        path_sites = write_sites_bed(config, chr_prefix, genome_release, tmp_dir)
        path_calls = call_sites(config, path_sites, tmp_dir)
        autosomal_fingerprint, autosomal_aafs = autosomal_snps_step_call(
            config,
            chr_prefix,
            genome_release,
            path_calls,
            config.output_fingerprint if config.write_vcf else None,
        )
        if not config.output_aafs:
            autosomal_aafs = None  # only pass to write_fingerprint if asked to do so
    return autosomal_fingerprint, autosomal_aafs


def write_fingerprint(
    config: Config,
    genome_release: str,
    sample: str,
    fingerprint: typing.Optional[np.array],
    aafs: typing.Optional[typing.List[float]],
    samtools_idxstats: typing.Optional[str]
):
    logger.info("Writing fingerprint to %s.npz ...", config.output_fingerprint)
    sections = "fingerprint,aafs" if aafs is not None else "fingerprint"
    header = np.array(
        [
            "ngs_chew_fingerprint",  # file identifier
            "2",  # file format version
            genome_release,  # genome release
            sample,  # sample name
            sections,  # sections in the file
        ]
    )
    np.savez_compressed(
        config.output_fingerprint, header=header, fingerprint=fingerprint, aafs=aafs,
        samtools_idxstats=samtools_idxstats
    )


def samtools_idxstats_step(config: Config):
    logger.info("Running samtools idxstats")
    cmd_idxstats = SAMTOOLS_IDXSTATS % {"input_bam": config.input_bam}
    logger.info("  command: %s", " ".join(shlex.split(cmd_idxstats)))
    result = subprocess.check_output(shlex.split(cmd_idxstats), text=True)
    result_log = result
    if len(result_log) > 100:
        result_log = result_log[:100] + "[...]"
    logger.info("Result is (truncated after 100 chars for log) \n%s", result_log)
    return result


def run(config: Config):
    if config.output_fingerprint.endswith(".npz"):
        config = attrs.evolve(config, output_fingerprint=config.output_fingerprint[: -len(".npz")])
    logger.info("Chewing NGS at %s", config.input_bam)

    sample, chr_prefix, genome_release = analyze_bam_header(config.input_bam, config.genome_release)
    # TODO: compare contigs in BAM and reference

    if config.step_autosomal_snps:
        autosomal_fingerprint, autosomal_aafs = autosomal_snps_step(
            config, chr_prefix, genome_release
        )
    else:
        autosomal_fingerprint, autosomal_aafs = None, None

    if config.step_samtools_idxstats:
        samtools_idxstats_out = samtools_idxstats_step(config)
    else:
        samtools_idxstats_out = None

    write_fingerprint(config, genome_release, sample, autosomal_fingerprint, autosomal_aafs, samtools_idxstats_out)

    logger.info("All done. Have a nice day!")
