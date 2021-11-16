# fingerprints from fastq
#%%
"""
SNPs=$(realpath /fast/projects/cubit/work/20.05/static_data/db/dbSNP/b150/GRCh37/All_20170710.vcf.gz)
SITES=$(realpath chew/data/GRCh37_sites.bed.gz)
OUTPUT=$(realpath ~/work/development/ngs-chew/chew/data/GRCh37.SNPs.bed.gz)
bedtools intersect -wa -a <(bcftools norm -m- $SNPs) -b <(zcat $SITES) | awk -v OF='\t' '{if (length($4)==1 && length($5)==1) print $1,$2,$2+1,$5}' | gzip > $OUTPUT
"""

import subprocess, shlex, pathlib, pandas as pd, numpy as np, tempfile as tmp, os
from collections import defaultdict
from logzero import logger
from chew.fingerprint import write_fingerprint

from argparse import Namespace

# load BED = GRCh37.SNPs.bed.gz
# expand for left, center, right
# shift frame of BED
# write to tmp file
# outf = bedtools getfasta
# split outf in left, center, right, write to tmp files (bed files)
# keep alt bases

def canonical_kmer(x):
    return min(x,x.translate(x.maketrans('ACGT','TGCA'))[::-1])

#%%
#  input:
SNPs_bed_path  = pathlib.Path("/home/memsonmi/DEV/ngs-chew/chew/data/GRCh37.SNPs.bed.gz")
reference_path = pathlib.Path("/media/memsonmi/ext_drive/hs37d5/hs37d5.fa")
read_path      = pathlib.Path("/media/memsonmi/ext_drive/ngs-chew/tests/ddata/fasta/bwa.cfs_1264-N1-DNA1-WES1.fasta")
# in args:
# cores = 2
# bf_size = "1G"
# k2mer_size = "200M" #roughly the number of targets?
# k
#%%
genome_release = "GRCh37"
def fastq_to_fingerprint(args, prefix, genome_release, paths_reads, prefix_fingerprint, path_reference):
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
    cmd = shlex.split(f"bedtools getfasta -tab -fi {reference_path} -bed {tmp_bed_file.name} -fo {hmers_file.name}")
    subprocess.check_call(cmd)
    # read h-mer sequences
    hmers_df = pd.read_csv(hmers_file.name,sep='\t',header=None)
    # dictionary h-mer: SNP id
    hmer_to_kmer_dict = {hmers_df.iloc[i,0]:str(D.iloc[i,0])+':'+str(D.iloc[i,1])+'-'+str(D.iloc[i,2]) for i in D.index}
    hmers_distinct_index = hmers_df.drop_duplicates(subset=[0]).index

    # Open new bed-file for writing.
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
            cmd_zcat = shlex.split(f"zcat {' '.join([str(r) for r in paths_reads])}")
            cmd_count = shlex.split(f"jellyfish count -m {2*size_k+1} -s {args.k2mer_size} --bf-size {args.bf_size} -C -t {args.cores} -o {tmpfcounts_file.name} --if {kmers_ref_file.name} {read_path}")
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
            cmd_zcat = shlex.split(f"zcat {' '.join([str(r) for r in paths_reads])}")
            cmd_count = shlex.split(f"jellyfish count -m {2*size_k+1} -s {args.k2mer_size} --bf-size {args.bf_size} -C -t {args.cores} -o {tmpfcounts_file.name} --if {kmers_ref_file.name} {read_path}")
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
            write_fingerprint(args,genome_release,args.sample,fingerprint)
#%%

