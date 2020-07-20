This is is the **multitax** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project

:Overview: Runs taxonomic analysis on a set of samples using sequana_taxonomy (kraken behing the scene)
:Input: A set of Fastq files
:Output: HTML report for each sample and a summary HTML report for all (multiqc +  dendogram)
:Status: draft
:Citation: Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~

You must install Sequana first::

    pip install sequana

Then, just install this package::

    pip install sequana_multitax


Usage
~~~~~

::

    sequana_pipelines_multitax --help
    sequana_pipelines_multitax --input-directory DATAPATH  --databases toydb

For the database, you will need to provide your own databases. You can check out
the documentation of kraken. The toydb here above is shipped wit sequana and
should work for demo. See sequana_taxonomy standalone for more help and
information. You can also checkout the sequana documentation (kraken module) 

This creates a directory with the pipeline and configuration file. You will then need 
to execute the pipeline::

    cd multitax
    sh multitax.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can 
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s multitax.rules -c config.yaml --cores 4 --stats stats.txt

Or use `sequanix <https://sequana.readthedocs.io/en/master/sequanix.html>`_ interface.

Requirements
~~~~~~~~~~~~

This pipelines requires the following executable(s):

- kraken and/or kraken2
- sequana_taxonomy


You cannot install both kraken1 and kraken2 together. We recommende to use the
latest version::

    conda install kraken2

.. image:: https://raw.githubusercontent.com/sequana/sequana_multitax/master/sequana_pipelines/multitax/dag.png


You can download databases from kraken website. We provide some databases on
github.com/sequana/resources. You can download a toy database as follows::

    sequana_taxonomy --download toydb

A more complete database is available here::

    sequana_taxonomy --download kraken_db1

The first time, a taxonomic database will be downloaded and stored locally in
.config/sequana/taxonomy.data file. You can updateit from time to time using::

    sequana_taxonomy --update-taxonomy


Details
~~~~~~~~~

This pipeline runs **sequana_taxonomy** (based on kraken) in parallel on the input fastq files (paired or not). 
A brief sequana summary report is also produced.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/sequana_multitax/master/sequana_pipelines/multitax/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 

Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
0.8.8     * fix plot summary dbs (sample names). Add options in schema+config
            file to tune the image if required.
          * HTML now includes links towards data that generates the top plots
          * fix case where zero sequences are found
          * check existence of input databases
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


