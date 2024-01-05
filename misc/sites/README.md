# NGS Chew Sites Related Computation

This directory contains the documentation on the `ngs-chew` sites are computed and verified as well as the generation of the PCA projection.

## Site Selection

The details can be found in the `select-auto` and `select-xy` directories.
To rerun it, obtain dbSNP GRCh37/GRCh38 files and `pip -r install /requirements/misc.txt`.
Then run the following:

```
PATH_DBSNP37=path/to/dbsnp-b151-grch37.vcf.gz \
PATH_DBSNP38=path/to/dbsnp-b151-grch38.vcf.gz \
    snakemake -c -p work/label/grch37.common.bed
```

## Autosomal Site Selection

We base the selection of our autosomal sites on the GRCh37 sites from [peddy](https://github.com/brentp/peddy).
However, we limit the sites to those that have an RSID and are also in GRCh38.
The sites are in the order of the rsid.
This allows to re-use the fingerprints for both releases with only small limitations.

## Gonosomal Site Selection

See [this issue](https://github.com/bihealth/ngs-chew/issues/20) for the original selection of GRCh37 sites.
We use a similar approach as for autosomal sites.
However, these sites are only used for prediction of gonomosomal karyotypes and not fingerprints.

## Population Group Analysis

The folder `pca` contains a Snakemake workflow that infers population groups / clusters from the thousand genotypes data at our selected sites.
