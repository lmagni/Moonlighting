# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 12:47:04 2021

Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it

This scrips perform a rough assessment of protein conservation across representative species spanning the tree of life.
The number of successful Blast queries is used as a proxy for the number of homologs.
Note: results are not meaningful for single entries (a reciprocal best hit approach should be used), but can be indicative for the whole dataset
"""


import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


"""customise paths
note: avoid input paths with spaces in name"""

path_to_blast = "path to blastp"
path_to_db = "path to database to run blast query"
   

datasets = {"moonlighting" : "path to file with list of moonlighting protein",
            "non-moonlighting" : "path to file with list of moonlighting protein"
            }

blast_outdir = "directory for blast outputs"
outdir = "directory for text output"
plotdir = "directory for output plots"


ref_proteomes = {
"Mus musculus": 10090,
"Gallus gallus": 9031,
"Danio rerio": 7955,
"Drosophila melanogaster": 7227,
"Caernorhabtitis elegans": 6239,
"Arabidopsis thaliana": 3702,
"Oryza sativa": 39947,
"Zea mays": 4577,
"Dyctostelyum discoideum": 44689,
"Paramecium tetraurelia": 5888,
"Neurospora crassa": 367110,
"saccharomyces cerevisiae": 559292,
"Clamydomonas reinhardtii": 272561,
"Psichomitrella patens": 3218,
"Escherichia coli": 83333,
"Bacillus subtilis": 224308,
"Pseudomonas aeruginosa": 208964,
"Mycobacterium tubercolosis": 83332,
"Geobacter solphurreducens": 243231,
"Mycoplasma genitalium": 243273,
"Methanococcus jannaschii": 243232,
"Halobacterium salinarium": 64091,
"Sulpholobus solfataricum": 273057    
}

groups = {
"metazoans" : range(0,5),
"unicellular_eukaryotes" : range(5,9),
"plants" : range(11,14),
"eukaryotes" : range(0,14),
"Bacteria" : range(14,20),
"Archaea" : range(20,23)
}


# blast search
def runblast(query, taxid, out):
    outfmt = "\"7 qacc sacc staxids\""
    max_target_seqs = 30
    evalue = 0.005
    commandline = path_to_blast + " -db " + path_to_db + " -query " + query + " -out " + out + " -taxids " + str(taxid) + " -evalue " + str(evalue) + " -parse_deflines -outfmt " + outfmt + " -max_target_seqs " + str(max_target_seqs)
    os.system(commandline)
    return


def parse_blast(out):
    # returns dict query: number of hits
    successful_queries = []
    with open(out, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                query = line.strip().split("\t")[0]
                successful_queries.append(query)               
    return dict(Counter(successful_queries))
"""select one of the following two functions: the first if the input for blast is a list of accessions
the second if the input is a fasta file"""


# if blast input is a list of accessions
def get_accessions(accession_list):
    accessions = []
    with open(accession_list, 'r') as file:
        for line in file:
            line = line.strip()
            if line:                
                accessions.append(line)
    return accessions

"""
# if blast input is a fasta file
def get_accessions(fasta):
    accessions = []
    with open(fasta, 'r') as file:
        for line in file:
            if line.startswith(">"):
                accession = line.split("|")[1]
                accessions.append(accession)
    return accessions
"""


"""select one of the following two functions: the first if you have not run blast already
otherwise choose the second function"""

# the input "fasta" can either be a fasta file or a list of accessions
def make_conservation_mat(fasta, ref_proteomes):
    queries = get_accessions(fasta)   # rows
    taxids = list(ref_proteomes.values())   # columns
    out_mat = np.zeros((len(queries), len(taxids)))    
    for j, taxid in enumerate(taxids):
        print("taxid: " + str(taxid))
        out = os.path.join(blast_outdir, "blastp_" + name + "_" + str(taxid) + ".txt")
        runblast(fasta, taxid, out)
        tax_dict = parse_blast(out)        
        for query, n_hits in tax_dict.items():
            i = queries.index(query)
            out_mat[i, j] = n_hits
    return queries, taxids, out_mat

"""
# if blast outputs are already available
def make_conservation_mat(fasta, ref_proteomes):
    queries = get_accessions(fasta)   # rows
    taxids = list(ref_proteomes.values())   # columns
    out_mat = np.zeros((len(queries), len(taxids)))    
    for j, taxid in enumerate(taxids):
        out = os.path.join(blast_outdir, "blastp_" + name + "_" + str(taxid) + ".txt")        
        tax_dict = parse_blast(out)
        for query, n_hits in tax_dict.items():
            i = queries.index(query)
            out_mat[i, j] = n_hits
    return queries, taxids, out_mat
"""

def conserved_in_group(mat, group_idx, homolog_threshold, species_threshold):
    group_mat = mat[:, group_idx]
    sp_consv = np.sum(group_mat >= homolog_threshold, axis=1)
    return sp_consv >= species_threshold

homologs_per_species = []
per_species = []
per_group = []
for name, dataset in datasets.items():
    queries, taxids, out_mat = make_conservation_mat(dataset, ref_proteomes)
    outfile = os.path.join(outdir, name + "_conservation_mat.txt")
    with open(outfile, 'w') as file:
        header = "#query/taxid:\t" + "\t".join([str(t) for t in taxids])
        file.write(header + "\n")
        for i, row in enumerate(list(out_mat)):
            file.write(queries[i] + "\t" + "\t".join([str(x) for x in row]) + "\n")
    
    for i, species in enumerate(list(ref_proteomes.keys())):
        homologs_per_species.append(pd.DataFrame.from_dict({'dataset': name, 'species' : species, 'value': np.ravel(out_mat[:, i])}))
        per_species.append(pd.DataFrame.from_dict({'dataset': name, 'species' : species, 'value': np.ravel(out_mat[:, i] > 0)}))
    
    for group, group_idx in groups.items():
        group_consv = conserved_in_group(out_mat, group_idx, 1, 1)
        per_group.append(pd.DataFrame.from_dict({'dataset': name, 'group' : group, 'value': np.ravel(group_consv)}))
    
        
per_species_df = pd.concat(axis=0, ignore_index=True, objs=per_species)
per_group_df = pd.concat(axis=0, ignore_index=True, objs=per_group)

homologs_per_species_df = pd.concat(axis=0, ignore_index=True, objs=homologs_per_species)


# make plots               
plot_1 = plt.figure(figsize=(1.5 * len(ref_proteomes), 5))
sns.barplot(data=per_species_df, x='species', y='value', hue='dataset', capsize=.2)
plt.savefig(os.path.join(plotdir, "species_barplot.png"), bbox_inches='tight')

plot_2 = plt.figure(figsize=(2.5 * len(groups), 10))
sns.barplot(data=per_group_df, x='group', y='value', hue='dataset', capsize=.2)
plt.savefig(os.path.join(plotdir, "group_barplot.png"), bbox_inches='tight')



plot_3 = plt.figure()
g = sns.catplot(data=homologs_per_species_df, kind="box", col='species', col_wrap=11, x='dataset', y='value', sharey=False, sharex=False, showfliers=False, height=5, aspect=1)
# change labels
[plt.setp(ax.texts, text="") for ax in g.axes.flat] 
g.set_titles(row_template = '{row_name}', col_template = '{col_name}')
g.add_legend()
plt.suptitle('GO pairwise correlations')
g.tight_layout(pad=2.0)    
plt.savefig(os.path.join(plotdir, "species_boxplot.png"), bbox_inches='tight')
