#! /usr/bin/env python3

import pandas as pd
import io
from io import StringIO
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict


#PREMIER NETTOYAGE DU FICHIER BLAST 

def parse_blast(filename):

    names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","taxids"]
    df = pd.read_csv(filename, sep="\t", header=None, names=names)
    #On vérfie que toute la colonne 'taxids' est bien de type string
    df['taxids'] = df['taxids'].astype(str)

    #On vérifie si les 1000 reads initiaux ont été blastés 
    len_diff = (1000 - len(df['qseqid'].unique()))/10
    #On calcule le % de reads non classés et ne figurant pas dans le fichier de sortie de Blast
    print(len_diff, "% des reads (des 1000 unclassified mis en input)  ne figurent pas dans les résultats de BLASTN")

    #Si l'on a plusieurs taxids pour une même ligne dans la colonne 'taxids', on ne garde que le premier (sachant qu'ils correspondent tous au même lineage)
    df['taxids'] = df['taxids'].str.partition(";")[0]

    #NDistribution des bitscores pour chaque read
    print("Distribution des bitscores pour chaque read:")
    df_2 = pd.crosstab(df['qseqid'],df['bitscore'])
    print(df_2)
    print("Distribution des valeurs de la colonne bitscore avant suppression d'une partie des doublons sur le critère du bitscore:")
    print(df['bitscore'].describe())

    #On ne garde que les hits avec le meilleur score (bitscore) pour chaque read (suppression d'une partie des doublons)
    qseqid_set = set(df['qseqid'])
    for qseqid in qseqid_set:
        bitscore_max = df.loc[(df['qseqid'] == str(qseqid))]['bitscore'].max()
        df.drop(df[((df['qseqid'] == str(qseqid)) & (df['bitscore'] < int(bitscore_max)))].index, inplace=True)

    #On vérifie que l'on a bien augmenté la valeur moyenne des bitscores puisque l'on n'a gardé que les hits avec les bitscores max pour chaque read
    print("Distribution des valeurs de la colonne bitscore après suppression d'une partie des doublons sur le critère du bitscore:")
    print(df['bitscore'].describe())

    return df


#RECUPERATION DU PARENT_ID, DU RANG ET DU NOM_SCIENTIFIQUE POUR CHAQUE TAXID EXISTANT

def dico(filename):

    #On prendra en entrée le fichier taxonomy.dat de sequana
    with open(filename, 'r') as f:
        lines = f.readlines()
        current_key = None
        Dico_TAXID = defaultdict(list)

        for line in lines:

            if " : " in line:
                deb = (line.split(" : "))
                (key, value) = (deb[0].strip(), deb[1].strip())

                if key == 'ID':
                    current_key = value
                else:
                    Dico_TAXID[current_key].append(value)

    return Dico_TAXID

 
#OPTIONNEL : RECUPERATION DES TAXIDS DEPUIS LES ACCESSION NUMBERS

def acctotaxids(filename):

    df = parse_blast(filename)
    Dico_TAXID = dico("/pasteur/zeus/projets/p02/Biomics/Users/dgarnier/Scripts_Analyses_BLAST_Output/taxonomy.dat")

    acc_set = set(df['Sacc'])
    #On utlise le fichier de correspondance entre les accession numbers et les taxids disponible sur le site de la ncbi
    with open ("/pasteur/zeus/projets/p02/Biomics/Users/dgarnier/Scripts_Analyses_BLAST_Output/nucl_gb.accession2taxid") as ACC2TAXID:
        Dico_ACC2TAXID = {}
        for line in ACC2TAXID:
            div = line.split()
            accnb = div[1]
            if accnb in acc_set:
                Dico_ACC2TAXID[accnb] = div[2]

    import pdb; pdb.set_trace()
    df['taxids'] = df.apply(lambda row: Dico_ACC2TAXID[row.Sacc], axis=1)

    df['PARENT_ID'] = df.apply(lambda row: Dico_TAXID[row.taxids][0], axis=1)
    df['RANK'] = df.apply(lambda row: Dico_TAXID[row.taxids][1], axis=1)
    df['SC_NAME'] = df.apply(lambda row: Dico_TAXID[row.taxids][2], axis=1)

    return df
 
 
#RECUPERATION DU LINEAGE COMPLET A PARTIR DES TAXIDS

def taxidstolineage(taxid_set):

    Dico_TAXID2LIN = defaultdict(list)

    Dico_TAXID = dico("/pasteur/zeus/projets/p02/Biomics/Users/dgarnier/Scripts_Analyses_BLAST_Output/taxonomy.dat")

    #On reconstruit le lineage complet pour chaque taxid
    for taxid in taxid_set:

        head_ranks = ("strain", "species", "genus", "family", "order", "class", "phylum", "superkingdom" ) #Peut se modifier en fonction de ce que l'on veut

        if taxid == "None":
            Dico_TAXID2LIN["None"]=["None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]

        elif taxid == "nan":
            Dico_TAXID2LIN["nan"]=["None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None", "None"]

        else:

            #Dico_ranks attribue à chaque head_rank son sub_group de ranks

            taxid_sup = taxid
            for head_rank in head_ranks:

                if Dico_TAXID[str(taxid_sup)][1] == head_rank:       #On est déjà au niveau du head_rank qui nous intéresse
                    value_rank = Dico_TAXID[str(taxid_sup)][2]    #On prend le nom scientifique au niveau du head rank
                    taxid_rank = taxid_sup
                    taxid_sup = Dico_TAXID[str(taxid_sup)][0]    #On remonte de head rank en head rank

                else:                                            #Il faut prendre en compte que certains rangs ne sont pas classés dans la hiérarchie des rangs mais que l'on peut tout de mâme éventuellement remonter et tomber sur le head_rank d'intérêt. Il s'agit des rangs suivants : 'no_rank', 'clade', 'genotype', 'pathogroup', 'serotype', 'serogroup'.
                    taxid_bis = Dico_TAXID[str(taxid_sup)][0]

                    for i in range(0,20):

                        if Dico_TAXID[str(taxid_bis)][1] != head_rank:   #On peut remonter à l'infini sans retomber sur aucun head_rank (notammment pour les 'no_rank') donc pas de while ici
                            taxid_bis = Dico_TAXID[str(taxid_bis)][0]
                            i = i+1
                            value_rank = "None"
                            taxid_rank = "None"

                        elif Dico_TAXID[str(taxid_bis)][1] == head_rank:
                            taxid_rank = taxid_bis
                            value_rank = Dico_TAXID[str(taxid_bis)][2]
                            taxid_sup = Dico_TAXID[str(taxid_bis)][0]
                            break

                Dico_TAXID2LIN[str(taxid)].append(taxid_rank)
                Dico_TAXID2LIN[str(taxid)].append(value_rank)

    return Dico_TAXID2LIN


#RECUPERATION DU LEAST COMMON ANCETSOR POUR TOUS LES GPS DE DOUBLONS

def get_LCA(filename):

    df = parse_blast(filename)
    qseqid_set = set(df['qseqid'])
    taxid_set = set(df['taxids'])
    Dico_TAXID2LIN = taxidstolineage(taxid_set)

    #On récupère le lineage complet pour chaque hit pour voir à quel rang se trouve le LCA
    df['Strain'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][0], axis=1)
    df['Species'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][2], axis=1)
    df['Genus'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][4], axis=1)
    df['Family'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][6], axis=1)
    df['Order'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][8], axis=1)
    df['Class'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][10], axis=1)
    df['Phylum'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][12], axis=1)
    df['Superkingdom'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][14], axis=1)
    #Pour le cas particulier des uncultured bacteria de toutes sortes (qui sont tjs au niveau de l'espèce)
    df['Species_Name'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxids][3], axis=1)

    Rank_Names = ['Strain', 'Species',  'Genus', 'Family',  'Order',  'Class', 'Phylum', 'Superkingdom']

    Dico_SEQID2LCA = defaultdict(list)

    #On récupère le LCA pour chaque subset de doublons (gp de qseqid) et on les stocke dans le Dico_SEQID2LCA
    for qseqid in qseqid_set:

        for rank in Rank_Names:

            taxid_unique = df.loc[(df['qseqid'] == str(qseqid))][rank].unique()
            #Pour le cas particulier des uncultured bacteria de toutes sortes (qui sont tjs au niveau de l'espèce)
            name_species_unique = df.loc[(df['qseqid'] == str(qseqid))]['Species_Name'].unique()

            if (len(taxid_unique) == 1):

                if taxid_unique != "None":
                    taxid_LCA = taxid_unique[0]
                    break

            elif (len(taxid_unique) > 1):

                #Cas particulier d'eurkaryotic synthetic construct et de l'homme --> on conserve l'homme
                if "111789" in taxid_unique and "9606" in taxid_unique:
                    taxid_LCA = "9604"
                    break

                #Cas des uncultered bacteria --> on conserve juste bacteria au niveau du superkingdom
                elif rank == "Superkingdom" and "uncultured" in str(name_species_unique):
                    taxid_LCA = "2"
                    break

                else:
                    taxid_LCA = "None"

        Dico_SEQID2LCA[str(qseqid)].append(str(taxid_LCA))

    df['taxid_LCA'] = df.apply(lambda row: Dico_SEQID2LCA[row.qseqid][0], axis=1)

    return df
    

#SUPPRESSION DES DOUBLONS ET EXPORTATION DES RESULTATS

def remove_duplicates(filename):

    df = get_LCA(filename)
    #On enlève les doublons en considérant que le LCA des doublons est le meilleur résultat que l'on puisse obtenir pour chaque read
    df = df.drop_duplicates(subset=['qseqid'])

    df = df.drop(columns=["sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","taxids","Strain","Species","Genus","Family","Order","Class","Phylum","Superkingdom"])

    #On récupère le lineage entier pour chaque résultat de chaque read
    taxid_set = set(df['taxid_LCA'])
    Dico_TAXID2LIN = taxidstolineage(taxid_set)

    df['Strain'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][1], axis=1)
    df['Species'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][3], axis=1)
    df['Genus'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][5], axis=1)
    df['Family'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][7], axis=1)
    df['Order'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][9], axis=1)
    df['Class'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][11], axis=1)
    df['Phylum'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][13], axis=1)
    df['Superkingdom'] = df.apply(lambda row: Dico_TAXID2LIN[row.taxid_LCA][15], axis=1)

    #On exporte les résultats au format .csv
    df.to_csv("{}_Results.csv".format(filename), index=False)
    print("Resultat final après suppression de tous les doublons : on a d'abord conservé les hits avec le bitscore le plus élevé pour chaque read puis on a considéré que le LCA des doublons restants est le meilleur résultat que l'on puisse obtenir pour chaque read:")
    print(df)

    return df
    

#VISUALISATION GRAPHIQUE AVEC KRONA

def krona(filename):

    #df = remove_duplicates(filename)
    df = pd.read_csv(filename,sep=",", header="infer") #si on reprend le fichier de sortie au format .csv de la fonction remove_duplicates
    #On ne garde que 5 rangs d'intérêt pour la visualisation graphique avec Krona (nombre de rangs max pour Krona = 5)
    df = df.drop(columns=["qseqid","bitscore","taxid_LCA","Species_Name","Phylum","Order","Strain"])

    #On met au bon format pour krona (ktImportText)
    df = df.groupby(['Superkingdom','Class','Family','Genus','Species']).size().reset_index(name='Count')
    df = df.reindex(columns=['Count','Superkingdom','Class','Family','Genus','Species'])
    print("Résultats au bon format pour Krona:")
    print(df)

    df.to_csv("{}".format(filename), index=False)

    #On met au format .txt
    with open("{}".format(filename), "r") as input_file:
        text_list = []
        for line in input_file.readlines():
            lines = line.split(",",6)
            text_list.append("".join(line))
    with open("{}_Results.txt".format(filename), "w") as output_file:
        for line in text_list:
            output_file.write(line.replace(",","\t"))


{}_1000_unclassified_1.tsv_Results.csv.format(filename) = remove_duplicates(filename)
{}_1000_unclassified_1.tsv_Results.txt.format(filename) = krona(filename)


