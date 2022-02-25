#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
import pandas as pd
from collections import defaultdict
from sequana.taxonomy import Taxonomy


def parse_blast(filename):
    """Parse the output of a blast file

    Note that we hard-code the requested columns both here and in the pipeline

    """

    names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "taxids",
    ]

    df = pd.read_csv(filename, sep="\t", header=None, names=names)

    # cast the taxids column is a string
    df["taxids"] = df["taxids"].astype(str)

    # On vérifie si les 1000 reads initiaux ont été blastés
    #len_diff = (Ntotal - len(df["qseqid"].unique())) / 10

    # Percentage of unclassified reads not in blast
    #print(f"percentage of input reads not classified by Blast: {len_diff}% ")

    # keep only the first taxids if there are several in the same row knowing
    # they belong to the same lineage
    df["taxids"] = df["taxids"].str.partition(";")[0]

    # Distribution of the bitscores for each read
    # print("Distribution des bitscores pour chaque read:")
    df_2 = pd.crosstab(df["qseqid"], df["bitscore"])

    # Keeping only the best bitscore for each read
    qseqid_set = set(df["qseqid"])
    for qseqid in qseqid_set:
        bitscore_max = df.loc[(df["qseqid"] == str(qseqid))]["bitscore"].max()
        df.drop(
            df[
                ((df["qseqid"] == str(qseqid)) & (df["bitscore"] < int(bitscore_max)))
            ].index,
            inplace=True,
        )

    return df



def read_taxonomy(filename=None):
    # Reads ID, parent ID, scientific names
    # Read the sequana taxonomy file

    if filename is None:
        tax = Taxonomy()
        filename = tax.database

    with open(filename, "r") as f:
        current_key = None
        TAXID = defaultdict(list)

        for line in f.readlines():

            if " : " in line:
                deb = line.split(" : ")
                (key, value) = (deb[0].strip(), deb[1].strip())

                if key == "ID":
                    current_key = value
                else:
                    TAXID[current_key].append(value)

    return TAXID


# OPTIONNEL : RECUPERATION DES TAXIDS DEPUIS LES ACCESSION NUMBERS

"""def acctotaxids(filename):

    df = parse_blast(filename)
    TAXID = read_taxonomy("taxonomy.dat")

    acc_set = set(df['Sacc'])
    #On utlise le fichier de correspondance entre les accession numbers et les taxids disponible sur le site de la ncbi
    with open ("nucl_gb.accession2taxid") as ACC2TAXID:
        ACC2TAXID = {}
        for line in ACC2TAXID:
            div = line.split()
            accnb = div[1]
            if accnb in acc_set:
                ACC2TAXID[accnb] = div[2]

    df['taxids'] = df.apply(lambda row: ACC2TAXID[row.Sacc], axis=1)
    df['PARENT_ID'] = df.apply(lambda row: TAXID[row.taxids][0], axis=1)
    df['RANK'] = df.apply(lambda row: TAXID[row.taxids][1], axis=1)
    df['SC_NAME'] = df.apply(lambda row: TAXID[row.taxids][2], axis=1)
    return df
"""
# RECUPERATION DU LINEAGE COMPLET A PARTIR DES TAXIDS


def taxidstolineage(taxid_set):

    TAXID2LIN = defaultdict(list)

    TAXID = read_taxonomy()

    # On reconstruit le lineage complet pour chaque taxid
    head_ranks = (
        "strain",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
    )  # Peut se modifier en fonction de ce que l'on veut

    for taxid in taxid_set:

        if taxid not in TAXID or taxid is None:
            TAXID2LIN[taxid] = ["None"] * 16
        elif taxid == "nan":
            TAXID2LIN["nan"] = ["None"] * 16
        else:
            # ranks attribue à chaque head_rank son sub_group de ranks
            taxid_sup = taxid
            for head_rank in head_ranks:
                if (
                    TAXID[str(taxid_sup)][1] == head_rank
                ):  # On est déjà au niveau du head_rank qui nous intéresse
                    value_rank = TAXID[str(taxid_sup)][
                        2
                    ]  # On prend le nom scientifique au niveau du head rank
                    taxid_rank = taxid_sup
                    taxid_sup = TAXID[str(taxid_sup)][
                        0
                    ]  # On remonte de head rank en head rank

                else:  # Il faut prendre en compte que certains rangs ne sont pas classés dans la hiérarchie des rangs mais que l'on peut tout de mâme éventuellement remonter et tomber sur le head_rank d'intérêt. Il s'agit des rangs suivants : 'no_rank', 'clade', 'genotype', 'pathogroup', 'serotype', 'serogroup'.
                    taxid_bis = TAXID[str(taxid_sup)][0]

                    for i in range(0, 20):

                        if (
                            TAXID[str(taxid_bis)][1] != head_rank
                        ):  # On peut remonter à l'infini sans retomber sur aucun head_rank (notammment pour les 'no_rank') donc pas de while ici
                            taxid_bis = TAXID[str(taxid_bis)][0]
                            i = i + 1
                            value_rank = "None"
                            taxid_rank = "None"

                        elif TAXID[str(taxid_bis)][1] == head_rank:
                            taxid_rank = taxid_bis
                            value_rank = TAXID[str(taxid_bis)][2]
                            taxid_sup = TAXID[str(taxid_bis)][0]
                            break

                TAXID2LIN[str(taxid)].append(taxid_rank)
                TAXID2LIN[str(taxid)].append(value_rank)

    return TAXID2LIN


def get_LCA(filename):
    """Get the least common ancestor including duplicated"""

    # scan the blast results
    df = parse_blast(filename)
    qseqid_set = set(df["qseqid"])
    taxid_set = set(df["taxids"])

    # convert the taxon to lineage
    TAXID2LIN = taxidstolineage(taxid_set)

    # We get back the full lineages for each hit to check the LCA
    df["Strain"] = df.apply(lambda row: TAXID2LIN[row.taxids][0], axis=1)
    df["Species"] = df.apply(lambda row: TAXID2LIN[row.taxids][2], axis=1)
    df["Genus"] = df.apply(lambda row: TAXID2LIN[row.taxids][4], axis=1)
    df["Family"] = df.apply(lambda row: TAXID2LIN[row.taxids][6], axis=1)
    df["Order"] = df.apply(lambda row: TAXID2LIN[row.taxids][8], axis=1)
    df["Class"] = df.apply(lambda row: TAXID2LIN[row.taxids][10], axis=1)
    df["Phylum"] = df.apply(lambda row: TAXID2LIN[row.taxids][12], axis=1)
    df["Superkingdom"] = df.apply(lambda row: TAXID2LIN[row.taxids][14], axis=1)

    # deal with uncultured bacteria 
    df["SpeciesName"] = df.apply(lambda row: TAXID2LIN[row.taxids][3], axis=1)

    Rank_Names = [
        "Strain",
        "Species",
        "Genus",
        "Family",
        "Order",
        "Class",
        "Phylum",
        "Superkingdom",
    ]

    SEQID2LCA = defaultdict(list)

    # We get back the LCA for each subset of duplicated entryies (gp de qseqid)
    # and store them in SEQID2LCA
    for qseqid in qseqid_set:

        for rank in Rank_Names:

            taxid_unique = df.loc[(df["qseqid"] == str(qseqid))][rank].unique()

            # deal with uncultured bacteria
            name_species_unique = df.loc[(df["qseqid"] == str(qseqid))][
                "SpeciesName"
            ].unique()

            if len(taxid_unique) == 1:

                if taxid_unique != "None":
                    taxid_LCA = taxid_unique[0]
                    break

            elif len(taxid_unique) > 1:

                # special case of eurkaryotic synthetic construct and human
                if "111789" in taxid_unique and "9606" in taxid_unique:
                    taxid_LCA = "9604"
                    break

                # special case of uncultered bacteria --> we keep the bacteria
                # at the level of superkingdom
                elif rank == "Superkingdom" and "uncultured" in str(
                    name_species_unique
                ):
                    taxid_LCA = "2"
                    break

                else:
                    taxid_LCA = "None"

        SEQID2LCA[str(qseqid)].append(str(taxid_LCA))

    df["taxid_LCA"] = df.apply(lambda row: SEQID2LCA[row.qseqid][0], axis=1)

    return df


def remove_duplicates(filename):
    """supress duplicated and save results

    We remove the duplicates. First, we keep the hits with best bitscore for
    each read. We consider only the LCA of the remaining bitscore

    """

    df = get_LCA(filename)

    # we remove the duplicates assuming the LCA is the best choice
    df = df.drop_duplicates(subset=["qseqid"])

    df = df.drop(
        columns=[
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "taxids",
            "Strain",
            "Species",
            "Genus",
            "Family",
            "Order",
            "Class",
            "Phylum",
            "Superkingdom",
        ]
    )

    # On récupère le lineage entier pour chaque résultat de chaque read
    taxid_set = set(df["taxid_LCA"])
    TAXID2LIN = taxidstolineage(taxid_set)

    df["Strain"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][1], axis=1)
    df["Species"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][3], axis=1)
    df["Genus"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][5], axis=1)
    df["Family"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][7], axis=1)
    df["Order"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][9], axis=1)
    df["Class"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][11], axis=1)
    df["Phylum"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][13], axis=1)
    df["Superkingdom"] = df.apply(lambda row: TAXID2LIN[row.taxid_LCA][15], axis=1)

    # On exporte les résultats au format .csv
    df.to_csv("{}_Results.csv".format(filename), index=False)

    return df


def krona(filename, output):
    """krona visulation"""

    df = pd.read_csv(filename, sep=",", header="infer")  # output of remove_duplicates

    # We keep only the first ranks
    df = df.drop(columns=["qseqid", "bitscore", "taxid_LCA", "SpeciesName", "Strain"])

    # Se the correct format for krona (ktImportText)
    df = (
        df.groupby(
            ["Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        )
        .size()
        .reset_index(name="Count")
    )
    df = df.reindex(
        columns=[
            "Count",
            "Superkingdom",
            "Phylum",
            "Class",
            "Order",
            "Family",
            "Genus",
            "Species",
        ]
    )

    df.to_csv("{}".format(filename), index=False)

    # Set the correct format
    with open(filename, "r") as input_file:
        text_list = []
        for line in input_file.readlines():
            lines = line.split(",", 6)
            text_list.append("".join(line))

    with open(output, "w") as output_file:
        for line in text_list:
            output_file.write(line.replace(",", "\t"))
