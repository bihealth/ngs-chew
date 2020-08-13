import gzip
import pathlib

import pysam
import vcfpy


#: Key to use for GRCh37 release.
RELEASE_37 = "GRCh37"

#: Key to use for GRCh38 release.
RELEASE_38 = "GRCh38"

#: Length of GRCh37 chr1.
LEN_CHR1_RELEASE_37 = 123 # XXX

#: Length of GRCh38 chr1.
LEN_CHR1_RELEASE_38 = 123 # XXX


class CouldNotGuessGenomeRelease(Exception):
    """Raised when the genome release could not be guessed."""


def guess_release_fasta(input_fasta, genome_release=None):
    prefix = ""
    release = ""
    with open(input_fasta + ".fai", "rt") as inputf:
        for line in inputf:
            arr = line.strip().split("\t")
            name = arr[0]
            length = int(arr[1])
            if name in ("chr1", "1"):
                prefix = "chr" if name.startswith("chr") else ""
                if length == LEN_CHR1_RELEASE_37:
                    release = RELEASE_37
                    break
                elif length == LEN_CHR1_RELEASE_38:
                    release = RELEASE_38
                    break
                else:
                    raise CouldNotGuessGenomeRelease("Could not guess genome release from length of %s: %s", name, ref_len)
        else:  # did not break out above
            raise CouldNotGuessGenomeRelease("Could not guess genome release, contig 1/chr1 not found in BAM header.")
    logger.info("Guessing genome release to be %s, prefix is %s", releae, repr(prefix))
    return prefix, release



def guess_release_bam(input_bam, genome_release=None):
    prefix = ""
    release = ""
    if genome_release:
        release = genome_release
        with pysam.AlignmentFile(input_bam, "rb") as bamf:
            if bamf.reference_names and bam.reference_names[0].startswith("chr"):
                prefix = "chr"
            else:
                prefix = ""
        logger.info("Using genome release %s, prefix is %s", genome_release, repr(prefix))
    else:
        with pysam.AlignmentFile(input_bam, "rb") as bamf:
            for i, name in enumerate(bamf.reference_names):
                if name in ("chr1", "1"):
                    prefix = "chr" if name.startswith("chr") else ""
                    ref_len = bamf.reference_lengths[i]
                    if ref_len == LEN_CHR1_RELEASE_37:
                        release = RELEASE_37
                        break
                    elif ref_len == LEN_CHR1_RELEASE_38:
                        release = RELEASE_38
                        break
                    else:
                        raise CouldNotGuessGenomeRelease("Could not guess genome release from length of %s: %s", name, ref_len)
            else:  # did not break out above
                raise CouldNotGuessGenomeRelease("Could not guess genome release, contig 1/chr1 not found in BAM header.")
        logger.info("Guessing genome release to be %s, prefix is %s", releae, repr(prefix))
        return prefix, release


def write_sites_bed(args, prefix, genome_release, tmp_dir):
    if args.sites_vcf:
        path_vcf = args.sites_vcf
    else:
        path_vcf = os.path.join(os.path.dirname(__file__), "data", "%s_sites.vcf.gz" % genome_release)
    path_bed = os.path.join(tmp_dir, "sites.bed")
    logger.info("Writing sites BED file to %s", path_bed)
    with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
        # First, detect whether to add or strip prefix.
        for header_line in vcf_reader.header.get_lines("contig"):
            offset = 0
            prefix = ""
            if line.startswith("chr") and not prefix == "chr":
                offset = 3
            elif not line.startswith("chr") and prefix != "chr":
                prefix = "chr"
        # Then, write out sites BED file.
        with open(path_bed, "wt") as outputf:
            for lineno, record in enumerate(vcf_reader):
                if not args.max_sites or lineno < args.max_sites:
                    print("%s%s\t%s\t%s" % (prefix, record.CHROM[offset:], record.POS - 1, record.POS), file=outputf)
                else:
                    break
    logger.info("Wrote %s sites", "{:,}".format(lineno))
    return path_bed
