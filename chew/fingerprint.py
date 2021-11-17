import contextlib
import gzip
import os
import shlex
import subprocess
import tempfile

import pandas as pd
from collections import defaultdict
import pathlib
import tempfile as tmp
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

def canonical_kmer(x):
    return min(x,x.translate(x.maketrans('ACGTNtacgtn','TGCANtgcan'))[::-1])

def fastq_to_fingerprint(args, prefix, genome_release, paths_reads, path_reference):
    SNPs_bed_path = pathlib.Path("data") / f"{genome_release}.SNPs.bed.gz"
    # prepare h-mer sequences (h = 4*k-1)
    D = pd.read_csv(SNPs_bed_path,header=None,compression='gzip',sep='\t')
    D[0] = D[0].apply(lambda x: prefix+str(x))
    index_original = [str(D.loc[i,0])+':'+str(D.loc[i,1])+'-'+str(D.loc[i,2]) for i in D.index]
    size_k = args.k
    D[1] -= 2 * size_k
    D[2] += 2 * size_k
    tmp_bed_file = tmp.NamedTemporaryFile()
    D.loc[:,:2].to_csv(tmp_bed_file.name,index=False,header=False,sep='\t')
    D[1] += 2 * size_k
    D[2] -= 2 * size_k
    hmers_file = tmp.NamedTemporaryFile()
    #### hmers_file = pathlib.Path("/home/memsonmi/Desktop/hmers.tsv")
    cmd = shlex.split(f"bedtools getfasta -tab -fi {path_reference} -bed {tmp_bed_file.name} -fo {hmers_file.name}")
    subprocess.check_call(cmd)
    # read h-mer sequences
    hmers_df = pd.read_csv(hmers_file.name,sep='\t',header=None)
    # dictionary h-mer: SNP id
    hmer_to_kmer_dict = {hmers_df.iloc[i,0]:str(D.iloc[i,0])+':'+str(D.iloc[i,1])+'-'+str(D.iloc[i,2]) for i in D.index}
    hmers_distinct_index = hmers_df.drop_duplicates(subset=[0]).index
    # Open new bed-files for writing.
    kmers_ref_file = tmp.NamedTemporaryFile()
    kmers_alt_file = tmp.NamedTemporaryFile()
    with open(kmers_ref_file.name,'w') as f:
        with open(kmers_alt_file.name,'w') as g:
            dict_kmer_id_ref = defaultdict(list)
            for key in hmers_distinct_index:
                id = hmers_df.loc[key,0]
                line = hmers_df.loc[key,1]
                #  create 3 kmers surrounding SNP
                kmers = (line[2*size_k:].rstrip(),line[size_k:3*size_k+1].rstrip(),line[:2*size_k+1].rstrip())
                kmers = map(canonical_kmer,kmers)
                #  write all 3 canonical k-mers to fasta formatted file
                #  id can be used later to order all k-mer counts to the original ordering
                for i,kmer in enumerate(kmers):
                    dict_kmer_id_ref[kmer].append(id)
                    print('>'+id+'_'+['l','c','r'][i],file=f)
                    print(kmer,file=f)
            dict_kmer_id_alt = defaultdict(list)
            for key in hmers_df.index:
                id = hmers_df.loc[key,0]
                line = hmers_df.loc[key,1]
                alt_base = D.iloc[key,3]
                Y = (alt_base+line[2*size_k+1:].rstrip(),
                    line[size_k:2*size_k].rstrip()+alt_base+line[2*size_k+1:3*size_k+1].rstrip(),
                    line[:2*size_k].rstrip()+alt_base)
                Y = map(canonical_kmer,Y)
                for i,y in enumerate(Y):
                    dict_kmer_id_alt[y].append(id)
                    print('>'+id+'_'+['l','c','r'][i],file=g)
                    print(y,file=g)
            # use jellyfish to count k-mers
            # ref first
            df_counts = pd.DataFrame(np.zeros((len(index_original),2),dtype=int),index = index_original, columns = ['ref','alt'])
            # --- count ref --- #
            tmpfcounts_file = tmp.NamedTemporaryFile()
            kmercounts_file = tmp.NamedTemporaryFile()
            # open reads
            cmd_zcat = shlex.split(f"zcat -f {' '.join([str(r) for r in paths_reads])}")
            cmd_count = shlex.split(f"jellyfish count -m {2*size_k+1} -s {args.k2mer_size} --bf-size {args.bf_size} -C -t {args.cores} -o {tmpfcounts_file.name} --if {kmers_ref_file.name} {paths_reads}")
            p_zcat = subprocess.Popen(cmd_zcat, stdout=subprocess.PIPE)
            p_count = subprocess.Popen(cmd_count, stdin=p_zcat.stdout)
            p_count.communicate()
            #logger.info('calling:',' '.join(cmd_count))
            cmd_dump = shlex.split(f"jellyfish dump -t -c -o {kmercounts_file.name} {tmpfcounts_file.name}")
            #logger.info('calling:',' '.join(cmd_dump))
            subprocess.check_call(cmd_dump)
            kmercounts_table = pd.read_csv(kmercounts_file.name,sep='\t',header=None,index_col=0)
            # count_table per original id (chr:start-end), hold counts: ref, alt
            for kmer in (kmercounts_table[kmercounts_table[1] > 0]).index.intersection(pd.Index(dict_kmer_id_ref.keys())):
                hmers = dict_kmer_id_ref[kmer]
                count = kmercounts_table.loc[kmer,1]
                for hmer in hmers:
                    id = hmer_to_kmer_dict[hmer]
                    df_counts.loc[id,'ref'] += count
            # --- count alt --- #
            tmpfcounts_file = tmp.NamedTemporaryFile()
            kmercounts_file = tmp.NamedTemporaryFile()
            cmd_zcat = shlex.split(f"zcat -f {' '.join([str(r) for r in paths_reads])}")
            cmd_count = shlex.split(f"jellyfish count -m {2*size_k+1} -s {args.k2mer_size} --bf-size {args.bf_size} -C -t {args.cores} -o {tmpfcounts_file.name} --if {kmers_ref_file.name} {paths_reads}")
            p_zcat = subprocess.Popen(cmd_zcat, stdout=subprocess.PIPE)
            p_count = subprocess.Popen(cmd_count, stdin=p_zcat.stdout)
            p_count.communicate()
            #logger.info('calling:',' '.join(cmd_count))
            cmd_dump = shlex.split(f"jellyfish dump -t -c -o {kmercounts_file.name} {tmpfcounts_file.name}")
            #logger.info('calling:',' '.join(cmd_dump))
            subprocess.check_call(cmd_dump)
            kmercounts_table = pd.read_csv(kmercounts_file.name,sep='\t',header=None,index_col=0)
            # count_table per original id (chr:start-end), hold counts: ref, alt
            for kmer in (kmercounts_table[kmercounts_table[1] > 0]).index.intersection(pd.Index(dict_kmer_id_alt.keys())):
                hmers = dict_kmer_id_alt[kmer]
                count = kmercounts_table.loc[kmer,1]
                for hmer in hmers:
                    id = hmer_to_kmer_dict[hmer]
                    df_counts.loc[id,'alt'] += count
            df_counts_array = df_counts.to_numpy().T
            # genotype calling - very simple majority decision
            fingerprint = np.array(
            [
                np.sum(df_counts_array,axis=0) > 5, # args.min_coverage
                df_counts_array[0] > df_counts_array[1]*2,
                df_counts_array[1] > df_counts_array[0]*2
            ],dtype=bool)
            a,b = df_counts_array[1],np.sum(df_counts_array,axis=0)
            allelic_fractions = np.divide(a, b,out=np.zeros(len(b),dtype=float), where=b != 0)
            return fingerprint, allelic_fractions

def generate_fingerprints_from_fastq(args):
    logger.info("Chewing NGS at %s", args.infq)


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
    logger.info("Chewing NGS at %s", args.inbam)
    with tempfile.TemporaryDirectory() as tmp_dir:
        prefix, genome_release = guess_release(args.inbam, args.genome_release)
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
    logger.info("Chewing NGS at %s", args.invcf)
    prefix, genome_release = guess_release(args.invcf, args.genome_release)
    path_calls = args.invcf
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
    if args.inbam:
        generate_fingerprints_from_bam(args)
    elif args.invcf:
        generate_fingerprints_from_vcf(args)
    elif args.infq:
        generate_fingerprints_from_fastq(args)
    else:
        logger.error("Input format unknown. Accepted formats are: bam, vcf.gz, fastq, fastq.gz.")

    logger.info("All done. Have a nice day!")
