.. image:: https://badge.fury.io/py/sequana-multitax.svg
     :target: https://pypi.python.org/pypi/sequana_multitax

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

.. image:: https://github.com/sequana/multitax/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/multitax/actions/workflows/main.yaml


MULTITAX — Multi-database Taxonomic Classification pipeline
============================================================

:Overview: Runs taxonomic analysis on a set of samples using sequana_taxonomy
           (Kraken2 under the hood), optionally followed by BLAST on
           unclassified reads.
:Input: A set of FastQ files (paired or single-end).
:Output: HTML report for each sample and a summary HTML report for all samples.
:Status: Production
:Citation: Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines,
           Journal of Open Source Software, 2(16), 352,
           `doi:10.21105/joss.00352 <https://doi.org/10.21105/joss.00352>`_

.. image:: https://raw.githubusercontent.com/sequana/multitax/main/sequana_pipelines/multitax/dag.png
   :alt: Pipeline DAG


Installation
------------

::

    pip install sequana-multitax

To upgrade an existing installation::

    pip install sequana-multitax --upgrade


Quick Start
-----------

**Step 1 — prepare the working directory**::

    sequana_multitax \
        --input-directory /path/to/reads \
        --databases /path/to/krakendb

This creates a ``multitax/`` working directory containing ``config.yaml`` and a
``multitax.sh`` launch script.

**Step 2 — review the configuration** (optional but recommended)::

    cd multitax
    cat config.yaml   # adjust parameters as needed

**Step 3 — run the pipeline**::

    sh multitax.sh


Taxonomic database
------------------

You will need one or more Kraken2 databases. You can download a toy database
for testing::

    sequana_taxonomy --download toydb

The pipeline also requires a taxonomy file stored in
``~/.config/sequana/taxonomy.dat``. Download it once with::

    sequana_multitax --update-taxonomy

Call this command again from time to time when unknown taxon IDs appear in the
HTML reports.

Multiple databases can be passed to run iterative classification::

    sequana_multitax \
        --input-directory /path/to/reads \
        --databases /path/to/virusdb /path/to/bacteriadb


Apptainer / Singularity
-----------------------

Every tool runs inside a pre-built container. Point ``--apptainer-prefix`` to a
shared directory so images are downloaded once and reused across projects::

    sequana_multitax \
        --input-directory /path/to/reads \
        --databases /path/to/krakendb \
        --apptainer-prefix ~/.sequana/apptainers

Pass extra bind mounts with ``--apptainer-args`` if your data lives outside
``$HOME``::

    --apptainer-args "-B /data:/data"

When running snakemake manually, include the singularity options::

    snakemake -s multitax.rules -c config.yaml --cores 4 \
        --use-singularity \
        --singularity-prefix ~/.sequana/apptainers \
        --singularity-args "-B /home:/home"


HPC / SLURM cluster
-------------------

On a cluster with SLURM, pass ``--profile slurm``::

    sequana_multitax \
        --input-directory /path/to/reads \
        --databases /path/to/krakendb \
        --profile slurm \
        --slurm-queue fast \
        --jobs 40 \
        --apptainer-prefix /shared/containers


BLAST on unclassified reads
---------------------------

Reads that remain unclassified after Kraken can optionally be BLASTed against a
local database::

    sequana_multitax \
        --input-directory /path/to/reads \
        --databases /path/to/krakendb \
        --store-unclassified \
        --do-blast-unclassified

This requires a local BLAST+ installation and a downloaded ``nt`` database.


Pipeline overview
-----------------

1. **Kraken2** — classify reads against one or more databases sequentially.
2. **Krona** — interactive pie charts per sample.
3. **[Optional] BLAST** — align unclassified reads against a nucleotide DB.
4. **MultiQC** — aggregated summary report across all samples.

Each sample produces an HTML report with a static pie chart (species
distribution; grey = unclassified) that links to an interactive Krona chart.

.. image:: https://raw.githubusercontent.com/sequana/multitax/main/doc/images/piechart.png
   :alt: Sample pie chart

When multiple databases are provided they are applied sequentially. The order
matters: reads classified by the first database are removed before the second
database is run.


Configuration file
------------------

After running ``sequana_multitax``, a ``config.yaml`` is created in the working
directory. Key sections:

- ``sequana_taxonomy`` — databases, confidence threshold, store_unclassified
- ``blast`` — enable/disable BLAST on unclassified reads
- ``multiqc`` — aggregated report settings

Full reference:
`config.yaml <https://raw.githubusercontent.com/sequana/multitax/main/sequana_pipelines/multitax/config.yaml>`_


Requirements
------------

- kraken2
- sequana_taxonomy
- krona


Changelog
---------

========= ====================================================================
Version   Description
========= ====================================================================
0.14.1    * fix dict-style config assignment (use dot-notation on _Namespace)
          * update README to follow sequana pipeline conventions
0.14.0    * updated container and sequana to fix issue with sequential
            analysis (several DBs)
0.13.0    * new containerisaton
0.12.2    * switch apptainer for sequana_taxonomy to the apptainer sequana
          * Fix version of sequana_wrappers to v23.12.5
          * add precommit
          * Fix dag to fix multiqc when using apptainers
0.12.1    * update apptainers
0.12.0    * Refactor to use new Click framework
0.11.1    * add missing import in the main script
          * add wrapper version in config
0.11.0    * use latest wrappers and graphivz apptainer
          * create and use a sequana-wrappers for the sequana_taxonomy rule
          * fix type when downloading taxonomy.dat
0.10.2    * add singularity containers
0.10.1    * fix blast run when no taxid is found and HTML report
0.10.0    * uses new sequana wrappers and framework
          * add ability to run blast on unclassified reads
          * handle case of empty FastQ files
0.9.2     * add --update-taxonomy DB option
          * add --store-unclassified option
0.9.1     * fix a logger issue
0.9.0     * fix plot summary dbs (sample names). Add options in schema+config
            file to tune the image if required.
          * HTML now includes links towards data that generates the top plots
          * fix case where zero sequences are found
          * check existence of input databases
          * add the --run argument
          * add multitax version in the header
          * add search box (Sequana feature) in the CSV tables
0.8.7     * Update HTML report: fix the title of images. include table with DB
            proportion. Text to explain images and reports
0.8.6     * A better report with new features from sequana.taxonomy
0.8.5     * fix typo in doc, factorise multiqc rule
0.8.4     * implement the --from-project option
0.8.3     * add the confidence option in sequana_taxonomy rule
          * improve html report
          * uses new sequana framework to speed up --help calls
0.8.2     * less stringent on requirements (mode warning)
          * fix input of the multiqc rule
0.8.1     Fix requirements.
0.8.0     **First release.**
========= ====================================================================
