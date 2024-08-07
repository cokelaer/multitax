
.. image:: https://badge.fury.io/py/sequana-multitax.svg
     :target: https://pypi.python.org/pypi/sequana_multitax


.. image:: https://github.com/sequana/multitax/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/multitax/actions/workflows/main.yaml

.. image:: https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C3.10-blue.svg
    :target: https://pypi.python.org/pypi/sequana
    :alt: Python 3.8 | 3.9 | 3.10

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

This is is the **multitax** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project

:Overview: Runs taxonomic analysis on a set of samples using sequana_taxonomy (kraken behing the scene)
:Input: A set of Fastq files
:Output: HTML report for each sample and a summary HTML report for all samples.
:Status: Production
:Citation: Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~

If you already have all requirements, you can install the packages using pip::

    pip install sequana_multitax --upgrade


Usage
~~~~~

In order to use this pipeline, you will need Kraken databases. Please see Kraken pages for help.

We will also need a taxonomic databases. In principle this pipeline will download the file in your home, in
/home/user/.config/sequana/taxonomy.dat one for all. Note, that with singularity, this file will be downloaded locally in your working directory for each analysis,except if it is found in your home, in which case a simple copy is performed.

::

    sequana_multitax --help
    sequana_multitax --input-directory DATAPATH  --databases toydb

For the database, you will need to provide your own databases. You can check out
the documentation of kraken. The toydb here above is shipped with sequana and
should work for demo. See sequana_taxonomy standalone for more help and
information. You can also checkout the sequana documentation (kraken module).


The Kraken final report and blast analysis (if set) will need a taxonomic file
stored in the sequana config directory (HOME/.config/sequana/taxonomy.dat). If
not already done, type this command::

    sequana_multitax --update-taxonomy

You may need to call this command from time to time if unknown taxon appears in
the HTML reports.


This creates a directory with the pipeline and configuration file. You will then need
to execute the pipeline::

    cd multitax
    sh multitax.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s multitax.rules -c config.yaml --cores 4 --stats stats.txt

Or use `sequanix <https://sequana.readthedocs.io/en/main/sequanix.html>`_ interface.


Usage with singularity::
~~~~~~~~~~~~~~~~~~~~~~~~~

With singularity, initiate the working directory as follows::

    sequana_multitax --use-singularity ...

Images are downloaded in the working directory but you can store then in a directory globally (e.g.)::

    sequana_multitax --use-singularity --singularity-prefix ~/.sequana/apptainers

and then::

    cd multitax
    sh multitax.sh

if you decide to use snakemake manually, do not forget to add singularity options::

    snakemake -s multitax.rules -c config.yaml --cores 4 --stats stats.txt --use-singularity --singularity-prefix ~/.sequana/apptainers --singularity-args "-B /home:/home"


Requirements
~~~~~~~~~~~~

This pipelines requires the following executable(s):

- kraken2
- sequana_taxonomy
- krona


.. image:: https://raw.githubusercontent.com/sequana/multitax/main/sequana_pipelines/multitax/dag.png

You can download databases from kraken website. We provide some databases on
https://github.com/sequana/resources. You can download a toy database as follows::

    sequana_taxonomy --download toydb

The first time, a taxonomic database will be downloaded and stored locally in
.config/sequana/taxonomy.dat file. You can update it from time to time using::

    sequana_taxonomy --update-taxonomy


Details
~~~~~~~~~

This pipeline runs **sequana_taxonomy** (based on kraken) in parallel on the input fastq files (paired or not).
A brief sequana summary report is also produced. For each sample, a HTML page is
reported with the following kind of image. This pie chart is a static image
summarizing the species found in your sample. Unclassified reads are in grey.
Colors correspond to a kingdom (green for viruses). If you click on the image,
you will be redirect to a more precise pie chart base on Krona pie chart, which
is more interactive.

.. image:: https://raw.githubusercontent.com/sequana/multitax/main/doc/images/piechart.png

The analysis is enterily based on Kraken tool. If several databases are
provided, they are run sequentially. This requires a careful interpretation of
the results. Indeed analysing data with viruses then bacteria may give different
results as compare to analysing with bacteria then viruses.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/multitax/main/sequana_pipelines/multitax/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file.

Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
0.12.2    * switch apptainer for sequana_taxonomy to the apptainer sequana
          * Fix version of sequana_wrappers to v23.12.5
          * add precommit
          * Fix dag to fix multiqc when using apptainers
0.12.1    * update apptainers
0.12.0    * Refactor to use new Click framework
0.11.1    * add missing import in the main script
          * add wrapper version in config
0.11.0    * use latest wrappers and graphivz apptainer
          * create and use a sequana-wrappers for the sequana_taxonomy ruke
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
