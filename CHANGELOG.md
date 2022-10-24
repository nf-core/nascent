# nf-core/nascent: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - 2022-10-20

### `Added`

- DSL2 conversion
- [[#28](https://github.com/nf-core/nascent/issues/28)] - Added DRAGMAP alignment
- [[#64](https://github.com/nf-core/nascent/pull/64)] - Added CHM13 igenomes config
- [[#39](https://github.com/nf-core/nascent/issues/39)] - Add PINTS for TSS identification
- [[#71](https://github.com/nf-core/nascent/issues/71)] - Add FASTP for adapter trimming
- [[#77](https://github.com/nf-core/nascent/issues/77)] - Add dedup subworkflow

### `Fixed`

- [[#33](https://github.com/nf-core/nascent/issues/33)] - groHMM works on full runs. Added the keep standard chromosomes function to standardize bam files.

### `Dependencies`

- Updated Nextflow version to `v21.10.6`

## v1.0.1 - 2020-03-03

Update to the container to meet the latest template requirements, and dependencies for new features in an upcoming PR (R and Picard tools).

## v1.0 - 2019-04-16

Initial release of nf-core/nascent, created with the [nf-core](http://nf-co.re/) template.
