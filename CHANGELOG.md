# nf-core/nascent: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0 - 2024-02-16

### Added

- [1494fff](https://github.com/nf-core/nascent/commit/1494fff2ecd8b498e19d7d0fa3f7ee0f71088ab4) - Update pipeline template to [nf-core/tools 2.12](https://github.com/nf-core/tools/releases/tag/2.12)
- [ec24cfb](https://github.com/nf-core/nascent/commit/ec24cfb2646904bcc78379a071b50b53a855b9a9) - Update pipeline template to [nf-core/tools 2.11.1](https://github.com/nf-core/tools/releases/tag/2.11.1)
- [15b8ff3](https://github.com/nf-core/nascent/commit/15b8ff31cd43ec2e330ce7958bd6d7d65529b4d3) - Update pipeline template to [nf-core/tools 2.11](https://github.com/nf-core/tools/releases/tag/2.11)
- [9712163](https://github.com/nf-core/nascent/commit/97121638eb77e175b912ff45a669426e532c5d7f) - Update pipeline template to [nf-core/tools 2.10](https://github.com/nf-core/tools/releases/tag/2.10)
- [[#130](https://github.com/nf-core/nascent/pull/130)] - AWS Mega Tests
- [[#132](https://github.com/nf-core/nascent/pull/132)] - Add Bowtie 2 as an aligner option

### Changed

- [[#117](https://github.com/nf-core/nascent/pull/117)] - Use fromSamplesheet

## v2.1.1 - 2023-02-17

### Fixed

- [[#109](https://github.com/nf-core/nascent/pull/109)] - Using AWS igenomes or passing a pre-made index fails because of an update to the module expects them to have a metamap. Generation was still working.

## v2.1.0 - 2023-02-15

### Added

- [[#94](https://github.com/nf-core/nascent/pull/94)] - Add a second BEDTools intersect step to allow filtering and overlapping in the same workflow.
- [[#101](https://github.com/nf-core/nascent/pull/101)] - Initialized nf-test

### Changed

- [[#103](https://github.com/nf-core/nascent/pull/103)] - Updated Modules

### Fixed

- [[841ae62](https://github.com/nf-core/nascent/commit/841ae62)] - Updated PINTS version from 1.1.6 to 1.1.8 ([Fixes an issue where PINTS fails if one of the predictions was empty](https://github.com/hyulab/PINTS/issues/12))
- [[#97](https://github.com/nf-core/nascent/pull/97)] - Add HOMER channels to fix error about "Missing workflow output parameter: homer_peaks" when homer is skipped
- Add missing DOIs (@apeltzer)

## v2.0.0 - 2022-10-24

### Added

- DSL2 conversion
- [[#28](https://github.com/nf-core/nascent/issues/28)] - Add DRAGMAP alignment
- [[#64](https://github.com/nf-core/nascent/pull/64)] - Add CHM13 igenomes config
- [[#39](https://github.com/nf-core/nascent/issues/39)] - Add PINTS for TSS identification
- [[#71](https://github.com/nf-core/nascent/issues/71)] - Add FASTP for adapter trimming
- [[#77](https://github.com/nf-core/nascent/issues/77)] - Add dedup subworkflow

### Fixed

- [[#33](https://github.com/nf-core/nascent/issues/33)] - groHMM works on full runs. Added the keep standard chromosomes function to standardize bam files.

### Dependencies

- Updated Nextflow version to `v21.10.6`

## v1.0.1 - 2020-03-03

Update to the container to meet the latest template requirements, and dependencies for new features in an upcoming PR (R and Picard tools).

## v1.0 - 2019-04-16

Initial release of nf-core/nascent, created with the [nf-core](http://nf-co.re/) template.
