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
import sys
import os
import argparse
import shutil
import subprocess

from sequana_pipetools.options import *
from sequana_pipetools.misc import Colors
from sequana_pipetools.info import sequana_epilog, sequana_prolog
from sequana_pipetools import SequanaManager

col = Colors()

NAME = "multitax"


class Options(argparse.ArgumentParser):
    def __init__(self, prog=NAME, epilog=None):
        usage = col.purple(sequana_prolog.format(**{"name": NAME}))
        super(Options, self).__init__(
            usage=usage,
            prog=prog,
            description="",
            epilog=epilog,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        # add a new group of options to the parser
        so = SlurmOptions()
        so.add_options(self)

        # add a snakemake group of options to the parser
        so = SnakemakeOptions(working_directory=NAME)
        so.add_options(self)

        so = InputOptions()
        so.add_options(self)

        so = GeneralOptions()
        so.add_options(self)

        pipeline_group = self.add_argument_group("pipeline")

        pipeline_group.add_argument('--update-taxonomy', action="store_true", 
            help="""To set the lineage of taxon ID, you need to update the taxonomic 
                DB used internally from time to time. """)

        pipeline_group.add_argument('--kraken-level', dest="kraken_level",
            default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

        pipeline_group.add_argument('--kraken-confidence', dest="kraken_confidence",
            type=float,
            default=0.05, help="""confidence parameter used with kraken2 databases only""")

        pipeline_group.add_argument("--databases", dest="databases", type=str,
            nargs="+", required="--update-taxonomy" not in sys.argv,
            help="""Path to a valid Kraken database(s). See sequana_taxonomy
                standaline to download some. You may use several, in which case, an
                iterative taxonomy is performed as explained in online sequana
                documentation. You may mix kraken1 and kraken2 databases""")

        pipeline_group.add_argument("--store-unclassified", default=False, action="store_true",
            help="Unclassified reads are stored in the output directories")

        pipeline_group.add_argument("--do-blast-unclassified", default=False,
            action="store_true", help="""blast unclassified read from kraken.
            Requires a local Blast DB and --stored-unclassified""")

        self.add_argument("--run", default=False, action="store_true",
            help="execute the pipeline directly")



    def parse_args(self, *args):
        args_list = list(*args)
        if "--from-project" in args_list:
            if len(args_list)>2:
                msg = "WARNING [sequana]: With --from-project option, " + \
                        "pipeline and data-related options will be ignored."
                print(col.error(msg))
            for action in self._actions:
                if action.required is True:
                    action.required = False
        options = super(Options, self).parse_args(*args)
        return options


def main(args=None):

    if args is None:
        args = sys.argv


    # whatever needs to be called by all pipeline before the options parsing
    from sequana_pipetools.options import before_pipeline
    before_pipeline(NAME)

    # option parsing including common epilog
    options = Options(NAME, epilog=sequana_epilog).parse_args(args[1:])

    if options.update_taxonomy:
        cmd = "sequana_taxonomy --update-taxonomy"
        p = subprocess.Popen(cmd.split())
        p.wait()
        sys.exit(0)


    # the real stuff is here
    manager = SequanaManager(options, NAME)

    # create the beginning of the command and the working directory
    manager.setup()
    from sequana import logger

    logger.setLevel(options.level)
    logger.name = "sequana_multitax"
    logger.info(f"Welcome to sequana_multitax  pipeline.")

    # fill the config file with input parameters
    if options.from_project is None:
        cfg = manager.config.config
        cfg.input_pattern = options.input_pattern
        cfg.input_readtag = options.input_readtag
        cfg.input_directory = os.path.abspath(options.input_directory)

        manager.exists(cfg.input_directory)

        cfg['sequana_taxonomy']['level'] = options.kraken_level
        cfg['sequana_taxonomy']['databases'] = options.databases
        for db in options.databases:
            manager.exists(db)
        cfg['sequana_taxonomy']['confidence'] = options.kraken_confidence
        cfg['sequana_taxonomy']['store_unclassified'] = options.store_unclassified

        if options.do_blast_unclassified:

            cfg['sequana_taxonomy']['store_unclassified'] = True
            cfg['blast']['do'] = True


    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()

    if options.run:
        subprocess.Popen(["sh", "{}.sh".format(NAME)], cwd=options.workdir)


if __name__ == "__main__":
    main()
