#!/data/SW/anaconda3/envs/xli/bin/python
"""
Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it


Some theory:
    
The correlation between two GO annotations is defined here as the probability of finding them together 
in less proteins than they actually are, assuming hat they were independent.

I.e. the correlation between two GO annotation GO1 and GO2 is the hypergeometric cdf H(k, M, n, N),
where:
    k = number of proteins annotated with GO1 and GO2
    M = number of all proteins         
    n = number of proteins annotated with GO1
    N = number of proteins annotated with GO2

This approach was proposed in: https://www.frontiersin.org/articles/10.3389/fgene.2015.00200/full
However, the clusterisation approach was not used by these authors
"""


import parse_uniprot_tsv as prs
import os
from time import time
from itertools import chain, combinations, combinations_with_replacement
import numpy as np
from scipy.stats import hypergeom


# paths
"""customise paths
note: path_to_intact mast be the path to a file with interacting protein pairs each on a separate line spaced by '\t' """

path_to_uniprot_tsv = "path to uniprot tab separated file"
path_to_intact = "path to list of intercting pairs"
outdir = "directory for output matrix"




### FUNCTIONS #################################################################


def zeros(nrows, ncolumns):    
    row = [0 for i in range(ncolumns)]    
    M = [list(row) for j in range(nrows)]
    return(M)

def ctime(start, end):    
    total = end-start
    days = total // 86400
    hours = total // 3600 % 24
    minutes = total // 60 % 60
    seconds = total % 60
    return str(days) + "d " + str(hours) + "h " + str(minutes) + "m " + str(seconds) + "s"


# clustering ------------------------------------------------------------------
def load_interactors(path_to_intact):
    
    # output list of interactor pairs
    interactors = []
    
    with open(path_to_intact, 'r') as file:       
        for line in file:
            line = line.strip()
            if line:
                split = line.strip().split("\t")
                #interactor_A = split[0]
                #interactor_B = split[1]
                # avoid repeated ones
                #if not (interactor_A, interactor_B) in interactors and not (interactor_B, interactor_A) in interactors:            
                interactors.append((split[0], split[1]))
    interactors = list(set(interactors))
    return interactors

# return dict entry : primary interactors
def batch_primary(entries, interactor_list):
    primary_interactors = {entry : set() for entry in entries}
    for (interactor_A, interactor_B) in interactor_list:        
        # do not append entry if it interacts with itself
        if interactor_A in entries and interactor_A != interactor_B:                
            primary_interactors[interactor_A].add(interactor_B)
        if interactor_B in entries and interactor_B != interactor_A:        
            primary_interactors[interactor_B].add(interactor_A)    
    return primary_interactors

# return nested dict entry : primary interactor : secondary interactor
def batch_secondary(entries, interactor_list):
    primary_interactors = batch_primary(entries, interactor_list)    
    all_primary = set()
    for interactors in primary_interactors.values():
        all_primary |= set(interactors)
    all_secodary_interactors = batch_primary(all_primary, interactor_list)
    secondary_interactors = {entry : {primary : all_secodary_interactors[primary] for primary in primary_interactors[entry]} for entry in entries}
    return secondary_interactors

def up_to_secondary(entries, interactor_list):
    primary_and_secondary_interactors = {entry : set() for entry in entries}
    secondary_dict = batch_secondary(entries, interactor_list)    
    for entry, struct in secondary_dict.items():
        for primary_interactor, secondary_interactors in struct.items():
            primary_and_secondary_interactors[entry].add(primary_interactor)
            primary_and_secondary_interactors[entry] |= secondary_interactors
            primary_and_secondary_interactors[entry].remove(entry)
    return primary_and_secondary_interactors

# getting GO annotations ------------------------------------------------------


# gene ontology annotatons to be considered
myfields = ["Entry", "Gene ontology (GO)", "Gene ontology (cellular component)", "Gene ontology (biological process)", "Gene ontology (molecular function)"]

myGOs = myfields[2:]  # to reduce time work only on GO subsets

proteome = prs.parse_uniprot_tsv(path_to_uniprot_tsv)

def get_GOs(entries):    
    raw_data = proteome.find(*entries, fields=myGOs)
    data = [set()] * len(myGOs)
    for y in raw_data:
        GOs = [x.split("; ") for x in y]   # to reduce time work only on GO subsets
        for i, GO_type in enumerate(GOs):
            data[i] |= set(GO_type)
    return [list(y) for y in data]    

def get_cluster_GOs(entries, interactor_list):
    interactors_dict = batch_primary(entries, interactor_list)
    #interactors_dict = up_to_secondary(entries, interactor_list)
    clusters = [{entry} | set(interactors) for entry, interactors in interactors_dict.items()]
    cluster_GOs = []
    for cluster in clusters:
        cluster_GOs.append(get_GOs(cluster))
    return cluster_GOs
    

### MAIN ######################################################################

print("start")
t_0 = time()

interactor_list = load_interactors(path_to_intact)

human_entries = proteome.list_entries()

cluster_all_GOs = get_cluster_GOs(human_entries, interactor_list)

for i, myGO in enumerate(myGOs):
    outname = "cluster_" + myGO.split("(")[1].split(")")[0].replace(" ", "_")
    out_counts_mat = os.path.join(outdir, outname + "_counts.mat")
    out_corr_mat = os.path.join(outdir, outname + "_correlations.mat")
    out_indexing = os.path.join(outdir, outname + "_indexing.mat")
    
    # make dictionary of all GOs per protein
    cluster_GOs = []
    for data in cluster_all_GOs:        
        if len(data[i]) > 1:
            cluster_GOs.append(data[i])
    
    # number of proteins with at least two annotations
    M = len(cluster_GOs)
   
    # make list of all GO entries sorted in alphabetic order
    #all_GOs = list(set(chain.from_iterable(GO_groupped)))
    all_GOs = []
    for subset in cluster_GOs:
        all_GOs += subset
    all_GOs = list(set(all_GOs))
    #all_GOs.sort()
    
    # make annotation to index dictionary to access numpy array
    # at the same time save indexing
    with open(out_indexing, 'w') as file:        
        idx = {}
        for k, GO in enumerate(all_GOs):
            idx[GO] = k
            file.write(str(k) + "\t" + GO + "\n")
    
    # number of possible annotations
    L = len(all_GOs)
        
    
    # COUNTS
    """
    count the number each pair of annotations is found within the same cluster
    """
    
    # counts will be a symmetric matrix
    counts = zeros(L,L)
    #counts = np.zeros((L,L))
            
    for GOs in cluster_GOs:
        cur_idx = [idx[GO] for GO in GOs]
        for i, j in combinations_with_replacement(cur_idx, 2):            
            counts[i][j] += 1
            counts[j][i] += 1
            #counts[i, j] = counts[j, i] = counts[i, j] + 1
        
    # save count matrix
    with open(out_counts_mat, 'w') as file:
        file.write("\n".join(["\t".join([str(x) for x in row]) for row in counts]))
    
    
    
    
    # CORRELATIONS
    # correlation matrix will be a symmetric matrix with diagonal 1      
    corrs = zeros(L,L)
    #corrs = np.zeros((L,L))
    for i, j in combinations_with_replacement(range(L), 2):
        corrs[i][j] = corrs[j][i] = hypergeom.cdf(counts[i][j], M,  counts[i][i], counts[j][j])
    
    
    # save correlation matrix
    with open(out_corr_mat, 'w') as file:
        file.write("\n".join(["\t".join([str(x) for x in row]) for row in corrs]))
    
    

t_1 = time()
print("finished\ntook: " + ctime(t_0, t_1))
