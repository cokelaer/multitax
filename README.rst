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

- kraken
- kraken2
- sequana_taxonomy


.. image::   https://github.com/sequana/sequana_multitax/raw/master/sequana_pipelines/multitax/dag.png

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
0.8.4     * implement the --from-project option
0.8.3     * add the confidence option in sequana_taxonomy rule
          * improve html report
          * uses new sequana framework to speed up --help calls
0.8.2     * less stringent on requirements (mode warning)  
          * fix input of the multiqc rule
0.8.1     Fix requirements.
0.8.0     **First release.**
========= ====================================================================


