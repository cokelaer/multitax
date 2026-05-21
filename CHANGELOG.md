# Changelog

All notable changes to this project will be documented in this file.

## [0.15.1] - 2026-05-21

### Added
- Improved test coverage with 8 new CLI tests (up from 3 to 11 total)
  - 5 new CLI config verification tests: paired reads, multiple databases, store-unclassified, kraken-confidence, blast options
  - Tests follow lora pipeline pattern using yaml.safe_load for config verification
  - Tests for paired-end data handling and multiple database configuration

### Fixed
- Fix multitax.rules: ensure unclassified.fastq output always exists
  - Wrapper shell script now touches unclassified.fastq file after sequana_taxonomy completes
  - Fixes missing output errors when store_unclassified=False

## [0.15.0] - 2024-08-07

### Added
- Initial stable release with Kraken2 and multiple database support
