# Changelog

### [0.9.2](https://www.github.com/bihealth/ngs-chew/compare/v0.9.1...v0.9.2) (2024-01-10)


### Bug Fixes

* release-please action ([#60](https://www.github.com/bihealth/ngs-chew/issues/60)) ([718bd1d](https://www.github.com/bihealth/ngs-chew/commit/718bd1d5fa8f38d6e61e2495870e40fa4ff4f221))

### [0.9.1](https://www.github.com/bihealth/ngs-chew/compare/v0.9.0...v0.9.1) (2024-01-10)


### Bug Fixes

* remove versioneer (Python 3.12+) ([#58](https://www.github.com/bihealth/ngs-chew/issues/58)) ([0b3876e](https://www.github.com/bihealth/ngs-chew/commit/0b3876e598207ee992ea2cec311353e86a08fef2))

## [0.9.0](https://www.github.com/bihealth/ngs-chew/compare/v0.8.1...v0.9.0) (2024-01-09)


### Features

* add "serve" command ([#24](https://www.github.com/bihealth/ngs-chew/issues/24)) ([#35](https://www.github.com/bihealth/ngs-chew/issues/35)) ([56f7078](https://www.github.com/bihealth/ngs-chew/commit/56f70788292a675ea5ee262ca344287d17fab94b))
* tests with larger (excerpt) datasets ([#55](https://www.github.com/bihealth/ngs-chew/issues/55)) ([8ee5421](https://www.github.com/bihealth/ngs-chew/commit/8ee54218fe95f453a75cca391b0a6bd3ee75f627))
* updating sites files for stable RSIDs ([#42](https://www.github.com/bihealth/ngs-chew/issues/42)) ([eaef39d](https://www.github.com/bihealth/ngs-chew/commit/eaef39dd1588939687f98d1daae8310ab344d62b))


### Bug Fixes

* Fills the stats line ([#57](https://www.github.com/bihealth/ngs-chew/issues/57)) ([f74bb70](https://www.github.com/bihealth/ngs-chew/commit/f74bb70cd1f54a42ce9187a053be4f9e228ae689))


### Documentation

* adding badges to README ([#52](https://www.github.com/bihealth/ngs-chew/issues/52)) ([68c5c82](https://www.github.com/bihealth/ngs-chew/commit/68c5c828b434c61454ae55c79b50160561a5c672))
* adding sphinx-based documentation ([#51](https://www.github.com/bihealth/ngs-chew/issues/51)) ([18ab48b](https://www.github.com/bihealth/ngs-chew/commit/18ab48b3b4be4501a216383b42dc9ac85a1ce7c7))

### [0.8.1](https://www.github.com/bihealth/ngs-chew/compare/v0.8.0...v0.8.1) (2023-07-18)


### Bug Fixes

* sites for GRCh38 ([#24](https://www.github.com/bihealth/ngs-chew/issues/24)) ([#39](https://www.github.com/bihealth/ngs-chew/issues/39)) ([79773ef](https://www.github.com/bihealth/ngs-chew/commit/79773ef0e829184e763aee61373ae46f6ea83ed5))

## [0.8.0](https://www.github.com/bihealth/ngs-chew/compare/v0.7.1...v0.8.0) (2023-01-16)


### Features

* roh calling with "bcftools roh" ([#25](https://www.github.com/bihealth/ngs-chew/issues/25)) ([a4af62c](https://www.github.com/bihealth/ngs-chew/commit/a4af62ced20bee37b4342d58aca48aaf9d269e68))

### [0.7.1](https://www.github.com/bihealth/ngs-chew/compare/v0.7.0...v0.7.1) (2023-01-10)


### Bug Fixes

* fixing chrX sites BED files ([#22](https://www.github.com/bihealth/ngs-chew/issues/22)) ([c0fb33e](https://www.github.com/bihealth/ngs-chew/commit/c0fb33e1a7261cfe234dba4938725be6755fb5f1))

## [0.7.0](https://www.github.com/bihealth/ngs-chew/compare/v0.6.0...v0.7.0) (2023-01-10)


### Features

* collect chrX SNP information for sex identification ([#12](https://www.github.com/bihealth/ngs-chew/issues/12)) ([#15](https://www.github.com/bihealth/ngs-chew/issues/15)) ([b210312](https://www.github.com/bihealth/ngs-chew/commit/b210312b40a3aad8e524a8765a8e8b2ef8d0aa8f))
* gather samtools idxstats output ([#13](https://www.github.com/bihealth/ngs-chew/issues/13)) ([aafa3cf](https://www.github.com/bihealth/ngs-chew/commit/aafa3cf7ca94005828c5ce2dd9927d3454291d81))
* interpret samtools idxstats output ([#18](https://www.github.com/bihealth/ngs-chew/issues/18)) ([50cd8fd](https://www.github.com/bihealth/ngs-chew/commit/50cd8fdd72d48bc7330ca5afdb5b108210cc3f75))
* using coding regions for chrX sites ([#20](https://www.github.com/bihealth/ngs-chew/issues/20)) ([#21](https://www.github.com/bihealth/ngs-chew/issues/21)) ([6c50688](https://www.github.com/bihealth/ngs-chew/commit/6c506887f809f9fc834e00f290a17097e6486a67))
* write out ngs-chew version to header ([#19](https://www.github.com/bihealth/ngs-chew/issues/19)) ([c0f2de5](https://www.github.com/bihealth/ngs-chew/commit/c0f2de5ff310410c37d3000533cf45d1c062a520))


### Bug Fixes

* use actual peddy formulat for relatedness ([#17](https://www.github.com/bihealth/ngs-chew/issues/17)) ([cea0b6c](https://www.github.com/bihealth/ngs-chew/commit/cea0b6cda3cfcec795172ef95571458ed74d3cab))

## [0.6.0](https://www.github.com/bihealth/ngs-chew/compare/v0.5.1...v0.6.0) (2023-01-02)


### Features

* allow writing allele fraction to fingerprint file ([#7](https://www.github.com/bihealth/ngs-chew/issues/7)) ([76f1511](https://www.github.com/bihealth/ngs-chew/commit/76f1511e2816ad08e37d76a35a0de02ba9e74c51))


### Bug Fixes

* formula for relatedness was off by *2 ([#9](https://www.github.com/bihealth/ngs-chew/issues/9)) ([3550dfb](https://www.github.com/bihealth/ngs-chew/commit/3550dfb0f35ae85b0e30de74cfda6c8db577bd94))

### [0.5.1](https://www.github.com/bihealth/ngs-chew/compare/v0.5.0...v0.5.1) (2022-12-21)


### Bug Fixes

* fix for the setup.py based build ([#4](https://www.github.com/bihealth/ngs-chew/issues/4)) ([255c848](https://www.github.com/bihealth/ngs-chew/commit/255c8482d1c9d14aadf15de95afaf97140e79205))

## [0.5.0](https://www.github.com/bihealth/ngs-chew/compare/v0.4.0...v0.5.0) (2022-12-21)


### Features

* adding README file ([#1](https://www.github.com/bihealth/ngs-chew/issues/1)) ([#2](https://www.github.com/bihealth/ngs-chew/issues/2)) ([620d487](https://www.github.com/bihealth/ngs-chew/commit/620d48747b845e93533a9f84aff082cc03cb2448))

## v0.4.0

- Adding `plot_var_het`.

## v0.3.0

- Adding option for writing out VCF file in fingerprinting.
- Adding command `plot_aab` command for visualizing B allele frequency from VCF files.

## v0.2.0

- Adding `plot_compare` command.

## v0.1.1

- Adding `--version` argument.

## v0.1.0

- Everything is new.
