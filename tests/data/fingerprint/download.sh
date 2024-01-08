#!/usr/bin/bash

# Obtain data for the fingerprint test.

set -x
set -euo pipefail

if [[ ! -e GRCh37.sites.bed ]]; then
    gzip -cd ../../../chew/data/GRCh37*.bed.gz \
    > GRCh37.sites.bed
fi

if [[ ! -e GRCh38.sites.bed ]]; then
    gzip -cd ../../../chew/data/GRCh38*.bed.gz \
    | sed -e 's/^/chr/g' \
    > GRCh38.sites.bed
fi

# Genome References
#
# GRCh37
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
# GRCh38
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

if [[ ! -e hs37d5.fa.fai ]]; then
    rm -f hs37d5.fa.gz
    aria2c -x 8 -s 8 -t 8 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
    pigz -c -d hs37d5.fa.gz \
    > hs37d5.fa
    samtools faidx hs37d5.fa.gz
    samtools faidx hs37d5.fa || rm -f hs37d5.fa.fai
fi
if [[ ! -e GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ]]; then
    rm -f GRCh38_full_analysis_set_plus_decoy_hla.fa*
    aria2c -x 8 -s 8 -t 8 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    bgzip -c GRCh38_full_analysis_set_plus_decoy_hla.fa \
    > GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
    samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
    samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa || rm -f GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
fi

# HG00138: British in England and Scotland, European Ancestry
#
# https://www.internationalgenome.org/data-portal/sample/HG00138
# 
# WES GRCh37
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00138/exome_alignment/HG00138.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam

if [[ ! -e HG00138_WES_GRCh37.excerpt.bam.bai ]]; do
    rm -f HG00138_WES_GRCh37.excerpt.bam*
    samtools view \
        --regions-file GRCh38.sites.bed \
        -O BAM \
        -o HG00138_WES_GRCh37.excerpt.bam \
        -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
        http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00138/exome_alignment/HG00138.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam
    samtools index HG00138_WES_GRCh37.excerpt.bam || rm -f HG00138_WES_GRCh37.excerpt.bam.bai
done

# WES GRCh38
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00138/exome_alignment/HG00138.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram

if [[ ! -e HG00138_WES_GRCh38.excerpt.bam.bai ]]; do
    rm -f HG00138_WES_GRCh38.excerpt.bam*
    samtools view \
        --regions-file GRCh38.sites.bed \
        -O BAM \
        -o HG00138_WES_GRCh38.excerpt.bam \
        -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
        http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00138/exome_alignment/HG00138.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
    samtools index HG00138_WES_GRCh38.excerpt.bam || rm -f HG00138_WES_GRCh38.excerpt.bam.bai
done

# WGS GRCh38
# ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240144/HG00138.final.cram

if [[ ! -e HG00138_WGS_GRCh38.excerpt.bam.bai ]]; do
    rm -f HG00138_WGS_GRCh38.excerpt.bam*
    samtools view \
        --regions-file GRCh38.sites.bed \
        -O BAM \
        -o HG00138_WGS_GRCh38.excerpt.bam \
        -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
        http://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240144/HG00138.final.cram
    samtools index HG00138_WGS_GRCh38.excerpt.bam || rm -f HG00138_WGS_GRCh38.excerpt.bam.bai
done
