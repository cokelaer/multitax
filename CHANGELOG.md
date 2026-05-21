# Changelog

All notable changes to this project will be documented in this file.

## [0.15.1] - 2026-05-21

### Added
- Comprehensive test coverage improvements (11 total tests)
  - 5 new CLI config verification tests for paired reads, multiple databases, store-unclassified, kraken-confidence, and blast options
  - 3 new full integration tests (test_full, test_full_paired, test_full_multiple_dbs) that execute the pipeline end-to-end
  - New test data directories: test/data/simple/ (unpaired) and test/data/paired/ (paired-end)

### Fixed
- Fix multitax.rules: make unclassified.fastq output conditional on store_unclassified config
  - Previously unconditionally declared as output, now uses dict unpacking for conditional Snakemake output declaration
  - Prevents pipeline failures when store_unclassified=False

## [0.15.0] - 2024-08-07

### Added
- Initial stable release with Kraken2 and multiple database support
