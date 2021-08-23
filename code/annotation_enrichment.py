# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 14:04:32 2021

Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it
"""

import parse_uniprot_tsv as prs
import os
from collections import Counter
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



"""customise paths and datasets"""

### PATHS #####################################################################
path_to_uniprot_tsv = "path to uniprot tab separated file"

# datasets

datasets = {"moonlighting" : "path to file with list of moonlighting protein",
            "non-moonlighting" : "path to file with list of moonlighting protein"
            }

# output
outdir = "directory where you want your outputs"
plot_name = "enrichment.png"
"""select the number of annotations to display in plot"""
num_annotations = 10

text_name = "enrichment.txt"
path_to_outfile = os.path.join(outdir, text_name)

                               
# annotatons to be considered
"""you can chose different sets of GO or ChEBI annotations"""

myfields = ["Gene ontology (GO)", "Gene ontology (cellular component)", "Gene ontology (biological process)", "Gene ontology (molecular function)"]
#myfields = ["ChEBI (Cofactor)", "ChEBI (Catalytic activity)"]




### MAIN ######################################################################

# retrieve entries and data from uniprot
proteome = prs.parse_uniprot_tsv(path_to_uniprot_tsv)

dataset_raw_data = {}
for name, path_to_dataset in datasets.items():
    entries = []
    with open(path_to_dataset, 'r') as dataset:
        for line in dataset:
            if not line.startswith("#"):
                entries.append(line.replace("\n", ""))    
    dataset_raw_data[name] = proteome.find(*entries, fields=myfields)


          

# create output file
outfile = open(os.path.join(outdir, text_name), 'w')
outfile.write("# Gene Ontology enrichment of human moonlighting and non-moonlighting proteins\n\n")
   
# create plot
ncolumns = 1
nrows = len(myfields)
myplot = plt.figure(figsize=(6.5 * ncolumns, 5 * nrows))



for i, myfield in enumerate(myfields):
    field_idx = myfields.index(myfield)    
     
    enrichment_dict = {}
    
    for j, (name, raw_data) in enumerate(dataset_raw_data.items()):
        outfile.write("# GO type: " + myfield + "\n dataset: " + name + "\n")
        # restrict data to the specific field
        data = [x[field_idx].split("; ") for x in raw_data]
        ref_data = [x[field_idx].split("; ") for x in raw_data if x[field_idx]]   # only consider cofactor binding proteins
        all_annotations = []
        for annotations in ref_data:
            all_annotations += annotations
       
        # count occurrences
        GO = dict(Counter(all_annotations))
        GO_ordered = dict(sorted(GO.items(), key=lambda item: item[1], reverse = True))
        
        # count dataset size
        N = len(data)
        
        dataset_enrichment = {}              
        for go_annotation, count in GO_ordered.items():
            enrichment = float(count) / N * 100
            dataset_enrichment[go_annotation] = enrichment
            outfile.write(go_annotation + ":\t" + str(count) + "\t" + str(enrichment) +"%\n")
            #print(go_annotation + ":\t" + str(count) + "\t" + str(enrichment) +"%")       
        
        enrichment_dict[name] = dict(dataset_enrichment)
        
    objects = []  
    
    annotations = list(enrichment_dict["moonlighting"].keys())[:num_annotations]
    for name, dataset_enrichment in enrichment_dict.items():
        for annotation, enrichment in dataset_enrichment.items():
            if annotation in annotations:
                objects.append(pd.DataFrame.from_dict({'annotation' : [annotation], 'enrichment': [enrichment], 'dataset': name}))
    
    # make subplots
    df = pd.concat(axis=0, ignore_index=True, objs=objects)
    plt.subplot(nrows, ncolumns, i + 1)
    sns.barplot(data=df, x='enrichment', y='annotation', hue='dataset')
    plt.title(myfield, fontsize=20)
    plt.xlabel('% of entries\n', fontsize=15)
    
    

myplot.tight_layout(pad=2.0)    
plt.savefig(os.path.join(outdir, plot_name), bbox_inches='tight')

outfile.close()
