#!/usr/bin/bash

# Generate fingerprint data.

set -x
set -euo pipefail

ngs-chew fingerprint \
    --reference GRCh38_full_analysis_set_plus_decoy_hla.fa.gz \
    --genome-release GRCh38 \
    --input-bam HG00138_WES_GRCh38.excerpt.bam \
    --output-fingerprint HG00138_WES_GRCh38.fingerprint.npz

ngs-chew fingerprint \
    --reference GRCh38_full_analysis_set_plus_decoy_hla.fa.gz \
    --genome-release GRCh38 \
    --input-bam HG00138_WGS_GRCh38.excerpt.bam \
    --output-fingerprint HG00138_WGS_GRCh38.fingerprint.npz

ngs-chew fingerprint \
    --reference hs37d5.fa.gz \
    --genome-release GRCh37 \
    --input-bam HG00138_WES_GRCh37.excerpt.bam \
    --output-fingerprint HG00138_WES_GRCh37.fingerprint.npz

ngs-chew compare --by-path *.npz
