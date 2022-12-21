# NGS Chew

> "Tasty, tasty NGS Data!"

NGS Chew is a growing toolbox of software for running quality control and sanity checks on NGS.
NGS chew can:

- Create a variant fingerprint file from BAM and VCF/BCF files.
    - The fingerprint files are store efficiently as compressed `numpy` arrays.
    - Optionally, allele balance information can be stored which enables advanced analysis downstream.
- Compare such fingerprint files to others to detect sample swaps and cryptic relationships.
- Analyze balance-enhanced fingerprint files for detecting cross-sample contamination.

## Quickstart

The following will create a `sample.npz` fingerprint file from the given BAM file.

```bash
ngs-chew fingerprint \
    --reference REFERENCE.fasta \
    --output-fingerprint sample.npz \
    --input-bam INPUT.bam \
    --genome-release GRCh37
```
