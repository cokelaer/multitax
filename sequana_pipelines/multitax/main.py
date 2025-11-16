# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import ast
import os
import subprocess
import sys

import click_completion
import rich_click as click
from sequana_pipetools import SequanaManager
from sequana_pipetools.options import *

click_completion.init()

NAME = "multitax"


help = init_click(
    NAME,
    groups={
        "Pipeline Specific": [
            "--databases",
            "--kraken-confidence",
            "--store-unclassified",
            "--do-blast-unclassified",
            "--kraken-level",
            "--update-taxonomy",
        ],
    },
)


# callback for --update-taxonomy option
def update_taxonomy(ctx, param, value):
    if value:
        cmd = "sequana_taxonomy --update-taxonomy"
        p = subprocess.Popen(cmd.split())
        p.wait()
        sys.exit(0)
    return value


# callback for --databases multiple arguments
def check_databases(ctx, param, value):
    if value:
        # click transform the input databases  (tuple) into a string
        # we need to convert it back to a tuple before checking the databases
        values = ast.literal_eval(value)
        for db in values:
            if not os.path.exists(db):
                click.echo(f"{db} does not exists. Check its path name")
                sys.exit(1)
    return ast.literal_eval(value)


@click.command(context_settings=help)
@include_options_from(ClickSnakemakeOptions, working_directory=NAME)
@include_options_from(ClickSlurmOptions)
@include_options_from(ClickInputOptions)
@include_options_from(ClickGeneralOptions)
@click.option(
    "--update-taxonomy",
    is_flag=True,
    callback=update_taxonomy,
    is_eager=True,
    help="""To set the lineage of taxon ID, you need to update the taxonomic DB used internally from time to time. """,
)
@click.option(
    "--kraken-confidence",
    "kraken_confidence",
    type=click.FLOAT,
    default=0.05,
    show_default=True,
    help="""confidence parameter used with kraken2 databases only""",
)
@click.option(
    "--databases",
    "databases",
    type=click.STRING,
    cls=OptionEatAll,
    callback=check_databases,
    required="--update-taxonomy" not in sys.argv,
    help="""Path to a valid Kraken database(s). See sequana_taxonomy
            standaline to download some. You may use several, in which case, an
            iterative taxonomy is performed as explained in online sequana
            documentation. You may mix kraken1 and kraken2 databases""",
)
@click.option(
    "--store-unclassified", default=False, is_flag=True, help="Unclassified reads are stored in the output directories"
)
@click.option(
    "--do-blast-unclassified",
    default=False,
    show_default=True,
    is_flag=True,
    help="""blast unclassified read from kraken. Requires a local Blast DB and --stored-unclassified""",
)
@click.option(
    "--keep-kraken-output",
    default=False,
    show_default=True,
    is_flag=True,
    help="""By default kraken output are deleted (to save space). use this flag to keep these files.""",
)
def main(**options):

    if options["from_project"]:
        click.echo("--from-project Not yet implemented")
        sys.exit(1)

    # the real stuff is here
    manager = SequanaManager(options, NAME)
    manager.setup()

    # fill the config file with input parameters

    cfg = manager.config.config

    # fills input_data, input_directory, input_readtag
    manager.fill_data_options()

    # specific options to the pipeline
    def fill_databases():
        if options["databases"]:
            for db in options["databases"]:
                if os.path.exists(db) is False:
                    click.echo(f"{db} not found. check its path.")
                    sys.exit(1)
            cfg["sequana_taxonomy"]["databases"] = [os.path.abspath(x) for x in options["databases"]]

    def fill_kraken_confidence():
        cfg["sequana_taxonomy"]["confidence"] = options["kraken_confidence"]

    def fill_store_unclassified():
        cfg["sequana_taxonomy"]["store_unclassified"] = options["store_unclassified"]
        cfg["sequana_taxonomy"]["keep_kraken_output"] = options["keep_kraken_output"]

    def fill_do_blast_unclassified():
        if options["do_blast_unclassified"]:
            cfg["sequana_taxonomy"]["store_unclassified"] = True
            cfg["blast"]["do"] = True

    if options["from_project"]:
        if "--databases" in sys.argv:
            fill_databases()
        if "--kraken-confidence" in sys.argv:
            fill_kraken_confidence()
        if "--store-unclassified" in sys.argv:
            fill_store_unclassified()
        if "--do_blast_unclassified" in sys.argv:
            fill_do_blast_unclassified()
    else:
        fill_databases()
        fill_kraken_confidence()
        fill_store_unclassified()
        fill_do_blast_unclassified()

    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()


if __name__ == "__main__":
    main()
