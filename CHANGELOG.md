# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.2.0] - 2026-03-07
### Added
- Co-infection detection logic to identify multiple RSV strains (different subtypes or clades) in a single sample.
- Automated read binning using BWA and Samtools to separate co-infection components into distinct read sets.
- Component-specific assembly and genotyping pipelines, creating separate output directories (e.g., `-comp1`, `-comp2`).
- Real-time console notifications in the main pipeline when co-infections are identified.
- Highlighted and styled co-infection notifications in individual sample report sections.

### Changed
- Refactored `Mapping.py` to support multi-component processing and standardized read naming.
- Updated `Report_functions.py` to correctly identify references and subtypes for co-infection components by inspecting the reference directory.
- Silenced IQTREE3 STDOUT by redirecting output to log files to keep the terminal output clean.
- Updated tree builder from FastTree to IQtree3 for more accurate phylogenetic analysis.
- Updated the QC calling function (`qc_call`) for more nuanced quality assessment.
- Improved PDF and HTML reports with more detailed QC status in summary tables.
- Refined QC calling logic to provide more concise "Good" status messages.
- Optimized PDF summary table column widths to prevent clipping and ensure all data fits within page margins.

### Fixed
- Resolved `NameError` in `Mapping.py` related to KMA logging variables.
- Resolved `FileNotFoundError` in reporting when processing co-infection components that utilized different reference genomes.
- Resolved `ValueError` in percentage formatting for the `QC rate` column in reports.
- Fixed logic to hide co-infection notifications for single-infection samples in the detailed report sections.

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