import contextlib
import gzip
import os
import shlex
import subprocess
import tempfile

from logzero import logger
import numpy as np
import vcfpy

#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Template for creating ``bcftools mpileup`` call.
TPL_PILEUP = r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"

#: Template for creating ``bcftools call`` call.
TPL_CALL = r"bcftools call -c -Oz -o %(calls)s"


def guess_release(input_bam, genome_release=None):
    prefix = ""
    if genome_release:
        logger.info("Using genome release %s", genome_release)
        return prefix, genome_release
    else:
        genome_release = "GRCh37"  # TODO: actually implement!
        logger.info("Guessing genome release to be %s", genome_release)
        return prefix, genome_release


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


def call_sites(args, path_sites, tmp_dir):
    logger.info("Performing variant calling at sites")
    path_vcf = os.path.join(tmp_dir, "calls.vcf.gz")
    cmd_pileup = TPL_PILEUP % {"sites": path_sites, **vars(args)}
    cmd_call = TPL_CALL % {"calls": path_vcf}
    logger.info("  mpileup: %s", " ".join(shlex.split(cmd_pileup)))
    logger.info("  call:    %s", " ".join(shlex.split(cmd_call)))
    p_pileup = subprocess.Popen(shlex.split(cmd_pileup), stdout=subprocess.PIPE)
    p_call = subprocess.Popen(shlex.split(cmd_call), stdin=p_pileup.stdout)
    p_call.wait()
    p_pileup.wait()
    return path_vcf


def vcf_to_fingerprint(args, chr_prefix, genome_release, path_calls, prefix_fingerprint):
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
    aafs = [aab for _, _, aab in sites.values()]
    fingerprint = np.array(
        [
            [dp > args.min_coverage for dp in depths],
            [gt != vcfpy.HOM_REF for gt in genotypes],
            [gt == vcfpy.HOM_ALT for gt in genotypes],
        ],
        dtype=bool,
    )
    return sample, fingerprint, aafs


def write_fingerprint(args, genome_release, sample, fingerprint, aafs):
    logger.info("Writing fingerprint to %s.npz ...", args.output_fingerprint)
    sections = "fingerprint,aafs" if aafs else "fingerprint"
    header = np.array(
        [
            "ngs_chew_fingerprint",  # file identifier
            "2",  # file format version
            genome_release,  # genome release
            sample,  # sample name
            sections,  # sections in the file
        ]
    )
    kwargs = {"header": header, "fingerprint": fingerprint}
    if aafs:
        kwargs["aafs"] = aafs
    np.savez_compressed(args.output_fingerprint, **kwargs)


def run(args):
    if args.output_fingerprint.endswith(".npz"):
        args.output_fingerprint = args.output_fingerprint[: -len(".npz")]
    logger.info("Chewing NGS at %s", args.input_bam)
    with tempfile.TemporaryDirectory() as tmp_dir:
        chr_prefix, genome_release = guess_release(args.input_bam, args.genome_release)
        path_sites = write_sites_bed(args, chr_prefix, genome_release, tmp_dir)
        path_calls = call_sites(args, path_sites, tmp_dir)
        sample, fingerprint, aafs = vcf_to_fingerprint(
            args,
            chr_prefix,
            genome_release,
            path_calls,
            args.output_fingerprint if args.write_vcf else None,
        )
        if not args.output_aafs:
            aafs = None  # only pass to write_fingerprint if asked to do so
        write_fingerprint(args, genome_release, sample, fingerprint, aafs)
    logger.info("All done. Have a nice day!")
