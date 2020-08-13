import gzip
import typing

import attr
from logzero import logger
import pysam
import tqdm

from .common import guess_release_fasta


@attr.s(auto_attribs=True, frozen=True)
class Site:
    """A genomic marker (SNV) site used for comparing genomes."""

    #: The name of the chromosome/contig that that the site is on.
    chromosome: str
    #: Position with 1-based coordinate.
    position: int
    #: Reference allele in VCF notation.
    reference: str
    #: Alternative allele in VCF notation.
    alternative: str
    
    def with_prefix(self, prefix):
        """Return ``Site`` having the ``"chr"`` prefix or not."""
        if prefix and not self.chromosome.startswith("chr"):
            return attr.evolve(chromosome="chr" + self.chromosome)
        elif not prefix and self.chromosome.startswith("chr"):
            return attr.evolve(chromosome=self.chromosome[3:])


def load_sites_vcf(args, prefix, genome_release):
    """Load sites VCF file, returns list of ``Site`` objects."""
    result = []

    if args.sites_vcf:
        path_vcf = args.sites_vcf
    else:
        path_vcf = os.path.join(os.path.dirname(__file__), "data", "%s_sites.vcf.gz" % genome_release)

    with vcfpy.Reader.from_path(path_vcf) as vcf_reader:
        # First, detect whether to add or strip prefix.
        for header_line in vcf_reader.header.get_lines("contig"):
            offset = 0
            prefix = False
            if line.startswith("chr") and not prefix == "chr":
                offset = 3
            elif not line.startswith("chr") and prefix != "chr":
                prefix = True
        # Then, load the sites.
        for lineno, record in enumerate(vcf_reader):
            if not args.max_sites or lineno < args.max_sites:
                result.append(Site(chromosome=record.CHROM, position=record.POS, reference=record.REF, alternative=record.ALT).with_prefix(prefix))
            else:
                break
            
    logger.info("Loaded %s sites", "{:,}".format(len(result)))
    return result


def extract_kmers(args, sites):
    """Extract the kmers from the given sites."""
    if args.output_kmers.endswith(".gz"):
        fn_open = gzip.open
    else:
        fn_open = open

    delta = args.kmer_length // 2
    with pysam.FastaFile(args.reference) as fastaf:
        with fn(args.output_kmers, "wt") as outf:
            print("#" + "\t".join(("SITENO", "CHROM", "POS", "REF", "ALT", "KMER")))
            for siteno, site in enumerate(tqdm(sites, "site(s)")):
                kmer = fastaf.fetch(site.chromosome, site.position - delta - 1, site.position + delta)
                if not kmer[delta] == site.reference:
                    raise Exception("Expected %s:%d to be %s but is %s!" % (site.chromosome, site.position, site.reference, kmer[delta]))
                print("\t".join(map(str, (siteno, site.chromosome, site.position, site.reference, site.alternative, kmer))), file=outf)


def run(args):
    logger.info("Collecting kmers for %s", args.reference)
    if not (pathlib.Path(args.reference) / ".fai").exists():
        logger.error("FAI file for %s does not exit!", args.reference))
        return 1
    if args.kmer_length % 2 ! = 1:
        logger.error("k-mer length must be uneven but is %d", args.kmer_length)
        return 1

    prefix, genome_release = guess_release_fasta(args.input_fasta, args.genome_release)
    sites = load_sites_vcf(args, prefix, genome_release)
    extract_kmers(args, sites)

    logger.info("All done. Have a nice day!")
