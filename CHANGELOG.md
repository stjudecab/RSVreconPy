# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - v0.2

### Changed
- Updated tree builder from FastTree to IQtree3 for more accurate phylogenetic analysis.
- Updated the QC calling function (`qc_call`) for more nuanced quality assessment.

## [0.1.2] - 2025-04-01

### Changed
- Updated the reference database for genotyping.
- Removed G-protein based genotyping to remain consistent with current research standards. More details can be found in the Nextclade data repository changelog.

## [0.1.1]

### Added
- Implemented a version control function to track pipeline and database versions.
- Added two new mutations to the RSV-B F-protein mutation list: I64T+K65E and N208D.

## [0.1.0]

### Added
- Initial release of the RSVrecon pipeline.