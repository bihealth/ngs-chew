import contextlib
import gzip
import os
import pathlib
import shlex
import subprocess
import tempfile

from logzero import logger
import numpy as np
import vcfpy

from .common import guess_release_bam, write_sites_bed

#: Template for creating ``bcftools mpileup`` call.
TPL_PILEUP = r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input_bam)s"

#: Template for creating ``bcftools call`` call.
TPL_CALL = r"bcftools call -c -Oz -o %(calls)s"


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


def vcf_to_fingerprint(args, prefix, genome_release, path_calls, prefix_fingerprint):
    logger.info("Reading sites BED...")
    if args.sites_bed:
        path_bed = args.sites_bed
        fn_open = gzip.open if path_bed.endswith(".gz") else open
    else:
        path_bed = os.path.join(os.path.dirname(__file__), "data", "%s_sites.bed.gz" % genome_release)
        fn_open = gzip.open
    with fn_open(path_bed, "rt") as inputf:
        sites = {}
        for line in inputf:
            arr = line.strip().split("\t")
            sites["%s%s:%s" % (prefix, arr[0], int(arr[1]) + 1)] = (0, 0)
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
                key = "%s%s:%s" % (prefix, record.CHROM, record.POS)
                if key in sites:
                    sites[key] = (record.INFO["DP"], record.call_for_sample[sample].gt_type)
    depths = [dp for dp, _ in sites.values()]
    genotypes = [gt for _, gt in sites.values()]
    fingerprint = np.array(
        [
            [dp > args.min_coverage for dp in depths],
            [gt != vcfpy.HOM_REF for gt in genotypes],
            [gt == vcfpy.HOM_ALT for gt in genotypes],
        ],
        dtype=bool,
    )
    return sample, fingerprint


def write_fingerprint(args, genome_release, sample, fingerprint):
    logger.info("Writing fingerprint to %s.npz ...", args.output_fingerprint)
    header = np.array(
        [
            "ngs_chew_fingerprint",  # file identifier
            "1",  # file format version
            genome_release,  # genome release
            sample,  # sample name
        ]
    )
    np.savez_compressed(args.output_fingerprint, header=header, fingerprint=fingerprint)


def run(args):
    if args.output_fingerprint.endswith(".npz"):
        args.output_fingerprint = args.output_fingerprint[: -len(".npz")]
    logger.info("Chewing NGS at %s", args.input_bam)
    if not (pathlib.Path(args.reference) / ".fai").exists():
        logger.error("FAI file for %s does not exit!", args.reference))
        return 1
    with tempfile.TemporaryDirectory() as tmp_dir:
        prefix, genome_release = guess_release_bam(args.input_bam, args.genome_release)
        path_sites = write_sites_bed(args, prefix, genome_release, tmp_dir)
        path_calls = call_sites(args, path_sites, tmp_dir)
        sample, fingerprint = vcf_to_fingerprint(
            args,
            prefix,
            genome_release,
            path_calls,
            args.output_fingerprint if args.write_vcf else None,
        )
        write_fingerprint(args, genome_release, sample, fingerprint)
    logger.info("All done. Have a nice day!")
