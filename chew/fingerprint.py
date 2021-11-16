import contextlib
import gzip
import os
import shlex
import subprocess
import tempfile

from logzero import logger
from tqdm import tqdm
import numpy as np
import vcfpy

#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Template for creating ``bcftools mpileup`` call.
TPL_PILEUP = r"bcftools mpileup -a AD,DP --threads 2 -I -R %(sites)s -f %(reference)s %(input)s"

#: Template for creating ``bcftools call`` call.
TPL_CALL = r"bcftools call -c -Oz -o %(calls)s"


def guess_release(input, genome_release=None):
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


def vcf_to_fingerprint(args, prefix, genome_release, path_calls, prefix_fingerprint):
    logger.info("Reading sites BED...")
    path_gz = os.path.join(os.path.dirname(__file__), "data", "%s_sites.bed.gz" % genome_release)
    with gzip.open(path_gz, "rt") as inputf:
        sites = {}
        for line in inputf:
            arr = line.strip().split("\t")
            sites["%s%s:%s" % (prefix, arr[0], int(arr[1]) + 1)] = (0, 0, 0, 0)
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
                    ad = record.call_for_sample[sample].data['AD']
                    sites[key] = (
                        record.INFO["DP"],
                        record.call_for_sample[sample].gt_type,
                        ad[0],
                        sum(ad))
    depths = [dp for dp,_,_,_ in sites.values()]
    genotypes = [gt for _,gt,_,_ in sites.values()]
    allelic_fractions = np.array([(adsum-ad) / adsum if adsum else 0.0 for _,_,ad,adsum in sites.values()],dtype=float)
    fingerprint = np.array(
        [
            [dp > int(args.min_coverage) for dp in depths],
            [gt != vcfpy.HOM_REF for gt in genotypes],
            [gt == vcfpy.HOM_ALT for gt in genotypes],
        ],dtype=bool
    )
    return sample, fingerprint, allelic_fractions

def fastq_to_fingerprint(args, prefix, genome_release, path_fastq, prefix_fingerprint):

    return sample, fingerprint, allelic_fractions

def write_fingerprint(args, genome_release, sample, fingerprint, allelic_fractions):
    logger.info("Writing fingerprint to %s.npz ...", args.output_fingerprint)
    header = np.array(
        [
            "ngs_chew_fingerprint",  # file identifier
            "2",  # file format version
            genome_release,  # genome release
            sample,  # sample name
        ]
    )
    np.savez_compressed(args.output_fingerprint, header=header, fingerprint=fingerprint, allelic_fractions=allelic_fractions)

def load_fingerprint(path):
    nparr = np.load(path)
    return nparr["header"][3], nparr["fingerprint"], nparr["allelic_fractions"]

def load_fingerprints(paths):
    logger.info("Loading fingerprints...")
    fps = {
        name: (fingerprint, allelc_fraction) for name, fingerprint, allelc_fraction in map(load_fingerprint, tqdm(paths))
    }
    logger.info("Loaded %d fingerprints", len(fps))
    return fps

def generate_fingerprints_from_bam(args):
    logger.info("Chewing NGS at %s", args.input)
    with tempfile.TemporaryDirectory() as tmp_dir:
        prefix, genome_release = guess_release(args.input, args.genome_release)
        path_sites = write_sites_bed(args, prefix, genome_release, tmp_dir)
        path_calls = call_sites(args, path_sites, tmp_dir)
        sample, fingerprint, allelic_fractions = vcf_to_fingerprint(
            args,
            prefix,
            genome_release,
            path_calls,
            args.output_fingerprint if args.write_vcf else None,
        )
        write_fingerprint(args, genome_release, sample, fingerprint, allelic_fractions)

def generate_fingerprints_from_vcf(args):
    logger.info("Chewing NGS at %s", args.input)
    prefix, genome_release = guess_release(args.input, args.genome_release)
    path_calls = args.input
    sample, fingerprint, allelic_fractions = vcf_to_fingerprint(
        args,
        prefix,
        genome_release,
        path_calls,
        args.output_fingerprint if args.write_vcf else None,
    )
    write_fingerprint(args, genome_release, sample, fingerprint, allelic_fractions)

def run(args):
    if args.output_fingerprint.endswith(".npz"):
        args.output_fingerprint = args.output_fingerprint[: -len(".npz")]
    if args.input.endswith(".bam"):
        generate_fingerprints_from_bam(args)
    elif args.input.endswith(".vcf.gz"):
        generate_fingerprints_from_vcf(args)
    else:
        logger.error("Input format unknown. Accepted formats are: bam, vcf.gz.")

    logger.info("All done. Have a nice day!")
