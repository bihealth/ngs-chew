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
import numpy.typing
import pysam
import vcfpy

import chew
from chew import roh
from chew.common import (
    CHROM_LENS_GRCH37,
    CHROM_LENS_GRCH38,
    CHROM_LENS_HG19,
    load_sites,
)

#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Template for creating ``bcftools mpileup`` call.
BCFTOOLS_PILEUP = (
    r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"
)

#: Template for creating ``bcftools call`` call.
BCFTOOLS_CALL = r"bcftools call --ploidy 2 -c -Oz -o %(calls)s"

#: Template for creating ``samtools idxstats`` call.
SAMTOOLS_IDXSTATS = r"samtools idxstats %(input_bam)s"


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
    #: Whether to collect autosomal SNP information.
    step_autosomal_snps: bool
    #: Whether to collect chrX SNP information.
    step_chrx_snps: bool
    #: Whether to collect "samtools idxstats" output.  Only applicable when extracting information
    #: from BAM files.
    step_samtools_idxstats: bool
    #: Whether to run "bcftools roh"
    step_bcftools_roh: bool


def analyze_bam_header(
    input_bam, genome_release: typing.Optional[str] = None
) -> typing.Tuple[str, str, str]:
    with pysam.AlignmentFile(input_bam, "rb") as samfile:
        header = samfile.header.to_dict()
        if "RG" in header:
            sample_names = {rg["SM"] for rg in header["RG"] if "SM" in rg}
        else:
            sample_names = set()

        samfile_chrom_lens = {
            record["SN"]: record["LN"]
            for record in (header.get("SQ") or [])
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
            raise InvalidInputWarning(
                "Could not determine genome release from BAM file's header. "
                "Please specify it using the --genome-release option."
            )
        logger.info("Guessing genome release to be %s", genome_release)

    if genome_release == "GRCh37":
        prefix = ""
    elif genome_release == "GRCh38":
        prefix = "chr"
    else:
        raise InvalidInputWarning(f"Invalid genome release {genome_release}")

    return sample_name, prefix, genome_release


def write_sites_bed(args, prefix, sites_suffix, genome_release, tmp_dir):
    path_gz = os.path.join(
        os.path.dirname(__file__), "data", f"{genome_release}_{sites_suffix}.bed.gz"
    )
    path_bed = os.path.join(tmp_dir, "sites.bed")
    logger.info("Writing sites BED file to %s", path_bed)
    arrs = []
    with gzip.open(path_gz, "rt") as inputf:
        for lineno, line in enumerate(inputf):
            if not args.max_sites or lineno < args.max_sites:
                arr = line.split("\t", 2)
                arrs.append((f"{prefix}{arr[0]}", int(arr[1]), arr[2]))
    arrs.sort()
    with open(path_bed, "wt") as outputf:
        for arr in arrs:
            print("\t".join(map(str, arr)).strip(), file=outputf)
    logger.info("Wrote %s sites", "{:,}".format(lineno))
    logger.info("  first = %s", "\t".join(map(str, arrs[0])).strip())
    logger.info("  last  = %s", "\t".join(map(str, arrs[-1])).strip())
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


def snps_step_call(
    config: Config,
    chr_prefix: str,
    genome_release: str,
    sites_suffix: str,
    path_calls: str,
    prefix_fingerprint: typing.Optional[str],
):
    bed_file = f"{genome_release}_{sites_suffix}.bed.gz"
    logger.info("Reading sites BED (%s)...", bed_file)
    sites = {
        "%s%s:%s" % (chr_prefix, site.chrom, site.pos): (0, 0, float("nan"))
        for site in load_sites(genome_release)
    }
    logger.info("Converting VCF to fingerprint...")
    with vcfpy.Reader.from_path(path_calls) as vcf_reader:
        if prefix_fingerprint:
            logger.info("Writing VCF to %s", f"{prefix_fingerprint}.vcf.gz")
            out_vcf = vcfpy.Writer.from_path(
                f"{prefix_fingerprint}.{sites_suffix}.vcf.gz", vcf_reader.header
            )
        else:
            logger.info("Not writing out VCF")
            out_vcf = contextlib.suppress()
        with out_vcf as vcf_writer:
            sample = vcf_reader.header.samples.names[0]
            for record in vcf_reader:
                if prefix_fingerprint:
                    vcf_writer.write_record(record)
                if chr_prefix == "chr" and not record.CHROM.startswith("chr"):
                    prefixed_chrom = f"chr{record.CHROM}"
                else:
                    prefixed_chrom = record.CHROM
                key = f"{prefixed_chrom}:{record.POS}"
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
    sites_suffix: str,
    genome_release: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        path_sites = write_sites_bed(config, chr_prefix, sites_suffix, genome_release, tmp_dir)
        path_calls = call_sites(config, path_sites, tmp_dir)
        fingerprint, aafs = snps_step_call(
            config,
            chr_prefix,
            genome_release,
            sites_suffix,
            path_calls,
            config.output_fingerprint if config.write_vcf else None,
        )
        if not config.output_aafs:
            aafs = None  # only pass to write_fingerprint if asked to do so
    return fingerprint, aafs


def write_fingerprint(
    config: Config,
    genome_release: str,
    sample: str,
    autosomal_fingerprint: typing.Optional[numpy.typing.NDArray],
    autosomal_aafs: typing.Optional[typing.List[float]],
    chrx_fingerprint: typing.Optional[numpy.typing.NDArray],
    chrx_aafs: typing.Optional[typing.List[float]],
    samtools_idxstats: typing.Optional[str],
    roh_txt: typing.Optional[roh.RohTxtContents],
):
    logger.info("Writing fingerprint to %s.npz ...", config.output_fingerprint)
    sections = []
    if autosomal_fingerprint is not None:
        sections.append("autosomal_fingerprint")
    if autosomal_aafs is not None:
        sections.append("autosomal_aafs")
    if chrx_fingerprint is not None:
        sections.append("chrx_fingerprint")
    if chrx_aafs is not None:
        sections.append("chrx_aafs")
    if samtools_idxstats is not None:
        sections.append("samtools_idxstats")
    if roh_txt is not None:
        sections.append("bcftools_roh")
    header = np.array(
        [
            "ngs_chew_fingerprint",  # file identifier
            "4",  # file format version
            chew.__version__,
            genome_release,  # genome release
            sample,  # sample name
            ",".join(sections),  # sections in the file
        ]
    )
    np.savez_compressed(
        config.output_fingerprint,
        header=header,
        autosomal_fingerprint=autosomal_fingerprint
        if autosomal_fingerprint is not None
        else np.zeros(0),
        autosomal_aafs=autosomal_aafs if autosomal_aafs is not None else np.zeros(0),
        chrx_fingerprint=chrx_fingerprint if chrx_fingerprint is not None else np.zeros(0),
        chrx_aafs=chrx_aafs if chrx_aafs is not None else np.zeros(0),
        samtools_idxstats=samtools_idxstats if samtools_idxstats is not None else np.zeros(0),
        roh_chroms=roh_txt.chroms if roh_txt is not None else np.zeros(0),
        roh_starts=roh_txt.starts if roh_txt is not None else np.zeros(0),
        roh_ends=roh_txt.ends if roh_txt is not None else np.zeros(0),
        roh_quals=roh_txt.quals if roh_txt is not None else np.zeros(0),
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


def bcftools_roh_step(sample: str, release: str, autosomal_fingerprint) -> roh.RohTxtContents:
    with tempfile.TemporaryDirectory() as tmpdir:
        logger.info("Prepare VCF files...")
        path_vcf = roh.write_vcf(str(tmpdir), sample, release, autosomal_fingerprint)
        path_txt = roh.run_roh(str(tmpdir), path_vcf)
        return roh.parse_roh_txt(path_txt)


def run(config: Config):
    if config.output_fingerprint.endswith(".npz"):
        config = attrs.evolve(config, output_fingerprint=config.output_fingerprint[: -len(".npz")])
    logger.info("Chewing NGS at %s", config.input_bam)

    sample, chr_prefix, genome_release = analyze_bam_header(config.input_bam, config.genome_release)
    # TODO: compare contigs in BAM and reference

    if config.step_autosomal_snps:
        autosomal_fingerprint, autosomal_aafs = autosomal_snps_step(
            config, chr_prefix, "sites", genome_release
        )
    else:
        autosomal_fingerprint, autosomal_aafs = None, None

    if config.step_chrx_snps:
        chrx_fingerprint, chrx_aafs = autosomal_snps_step(
            config, chr_prefix, "sitesX", genome_release
        )
    else:
        chrx_fingerprint, chrx_aafs = None, None

    if config.step_samtools_idxstats:
        samtools_idxstats_out = samtools_idxstats_step(config)
    else:
        samtools_idxstats_out = None

    if config.step_bcftools_roh:
        roh_txt_contents = bcftools_roh_step(
            sample=sample, release=genome_release, autosomal_fingerprint=autosomal_fingerprint
        )
    else:
        roh_txt_contents = None

    write_fingerprint(
        config,
        genome_release,
        sample,
        autosomal_fingerprint,
        autosomal_aafs,
        chrx_fingerprint,
        chrx_aafs,
        samtools_idxstats_out,
        roh_txt_contents,
    )

    logger.info("All done. Have a nice day!")
