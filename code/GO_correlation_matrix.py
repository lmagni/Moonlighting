#!/data/SW/anaconda3/envs/xli/bin/python
"""
Created on Fri Aug 13 16:49:47 2021
@author: Lorenzo

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
"""





import parse_uniprot_tsv as prs
import os
from time import time
from itertools import combinations_with_replacement
from scipy.stats import hypergeom



"""customise paths"""
# paths
path_to_uniprot_tsv = "path to uniprot tab separated file"
outdir = "directory for output matrix"


# gene ontology annotatons to be considered
myfields = ["Entry", "Gene ontology (GO)", "Gene ontology (cellular component)", "Gene ontology (biological process)", "Gene ontology (molecular function)"]
go_idx = 2 # to reduce time work only on GO subsets
myGOs = myfields[go_idx:]  




# numpy cannot be used because the matrices are too big
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



print("start")
t_0 = time()

# get proteome-wide info
proteome = prs.parse_uniprot_tsv(path_to_uniprot_tsv)
all_raw_data = proteome.find(fields=myfields)
protein_data_dict = {}
for data in all_raw_data:
    entry = data[0]
    #GOs = [x.split("; ") for x in data[1:]]
    GOs = [x.split("; ") for x in data[go_idx:]]
    protein_data_dict[entry] = GOs


for i, myGO in enumerate(myGOs):
    outname = myGO.split("(")[1].split(")")[0].replace(" ", "_")
    out_counts_mat = os.path.join(outdir, outname + "_counts.mat")
    out_corr_mat = os.path.join(outdir, outname + "_correlations.mat")
    out_indexing = os.path.join(outdir, outname + "_indexing.mat")
    
    # make dictionary of all GOs per protein
    protein_GOs = []
    for entry, data in protein_data_dict.items():        
        if len(data[i]) > 1:
            protein_GOs.append(data[i])
    
    # number of proteins with at least two annotations
    M = len(protein_GOs)
   
    # make list of all GO entries sorted in alphabetic order
    #all_GOs = list(set(chain.from_iterable(GO_groupped)))
    all_GOs = []
    for subset in protein_GOs:
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
        
    
    # counts will be a symmetric matrix
    counts = zeros(L,L)
            
    for GOs in protein_GOs:
        cur_idx = [idx[GO] for GO in GOs]
        for i, j in combinations_with_replacement(cur_idx, 2):            
            counts[i][j] += 1
            counts[j][i] += 1
            #counts[i, j] = counts[j, i] = counts[i, j] + 1
        
    # save count matrix
    with open(out_counts_mat, 'w') as file:
        for row in counts:
            file.write("\t".join([str(element) for element in row]) + "\n")
    
    
    
    
    # CORRELATIONS
    # correlation matrix will be a symmetric matrix with diagonal 1      
    corrs = zeros(L,L)
    for i, j in combinations_with_replacement(range(L), 2):
        corrs[i][j] = corrs[j][i] = hypergeom.cdf(counts[i][j], M, counts[i][i], counts[j][j])
        
        
    # save correlation matrix
    with open(out_corr_mat, 'w') as file:
        for row in corrs:
            file.write("\t".join([str(element) for element in row]) + "\n")
   
    
t_1 = time()
print("finished\ntook: " + ctime(t_0, t_1))
