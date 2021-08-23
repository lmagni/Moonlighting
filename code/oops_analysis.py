# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 00:54:57 2021

Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it
"""

import parse_uniprot_tsv as prs
import os
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""the OOPS and oligo(dT)_RBP_Capture databases are available in the supplementary material of:
    https://doi.org/10.1038/s41587-018-0001-2
    you can dowload them in xml format and save them as a tsv file from Excell"""
    
path_to_oops = "path to OOPS database in text format"
path_to_capture = "path to oligo(dT)_RBP_Capture file in text format"


path_to_uniprot_tsv = "path to uniprot tab separated file"
path_to_moonlighting_entries = "path to file with list of moonlighting protein"


# datasets
datasets = {"moonlighting" : "path to file with list of moonlighting protein",
            "non-moonlighting" : "path to file with list of moonlighting protein",
            "any other dataset": "path to dataset"
            }

"""if you also ant a random dataset from Uniprot set random=True and customise size"""
random=False
dataset_size = 2000


outdir = "directory for text output"
out_name = "oops_moonprot_multitask.txt"
plotdir = "directory for output plots"
enrichment = "oops_enrichment_datasets.png"
heatmap = "oops_heatmap.png"

oops_dict = {}
with open(path_to_oops) as oops:
    # skip first two lines
    next(oops)
    next(oops)
    for line in oops:
        # file headers Uniprot_ID	HEK293	MCF10A	U2OS	GO_RBP	RBP_Capture        
        split = line.strip().split("\t")
        if len(split) == 6:
            oops_dict[split[0]] = [int(x) for x in split[1:4]]

capture_dict = {}
with open(path_to_capture) as capture:
    # skip first two lines
    next(capture)
    next(capture)
    for line in capture:
        # file headers master_protein	CL	NC        
        split = line.strip().split("\t")
        if len(split) == 3:
            capture_dict[split[0]] = [float(x)/3 for x in split[1:]]

moonlighting_entries = []
with open(path_to_moonlighting_entries, 'r') as file:
    moonlighting_entries = [line.replace("\n", "") for line in file.readlines()]



# merge dicts into array moon_dataframe
fields = ["OOPS HEK293", "OOPS MCF10A", "OOPS U2OS", "oligo(dT)_RBP_Capture CL", "oligo(dT)_RBP_Capture NC"]
moon_data = []
# get moon_data froom oops dict
for protein in moonlighting_entries:
    if protein in oops_dict.keys():
        oops_moon_data = list(oops_dict[protein])
    else:
        oops_moon_data = [0, 0, 0]
    if protein in capture_dict.keys():
        capture_moon_data = list(capture_dict[protein])
    else:
        capture_moon_data = [0, 0]
    moon_data.append(oops_moon_data + capture_moon_data)

moon_data = np.array(moon_data)


### heatmap ###################################################################

moon_data_binary = moon_data > 0
moon_data_rot = np.rot90(moon_data_binary)
dataframe = pd.DataFrame(moon_data_rot, index = fields[::-1], columns = moonlighting_entries)

heatmap_fig = plt.figure(figsize=(0.5 * len(moonlighting_entries), 0.7 * len(fields)))
#ax = sns.heatmap(dataframe, annot=True, cmap='vlag')
ax = sns.heatmap(dataframe, annot=True, vmin=0, vmax=1, cmap='vlag', linewidths=.5)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize = 10)
ax.set_yticklabels(ax.get_yticklabels(), fontsize = 15)

plt.savefig(os.path.join(plotdir, heatmap), bbox_inches='tight')


### make datasets  ############################################################

# retrieve entries from datasets
dataset_entries = {}
for name, path_to_dataset in datasets.items():
    entries = []
    with open(path_to_dataset, 'r') as dataset:
        for line in dataset:
            if not line.startswith("#"):
                entries.append(line.replace("\n", ""))
    dataset_entries[name] = entries

# create random dataset from uniprot from uniprot
if random:
    proteome = prs.parse_uniprot_tsv(path_to_uniprot_tsv)
    #all_entries = proteome.list_entries()
    dataset_entries["random"] = proteome.sample_field(dataset_size)


### bootstrap datasets  #######################################################
sample_size = len(moonlighting_entries)
sample_number = 1000

# function to bootstrap datasets
def sample(mylist, sample_size, sample_number):
        all_idx = np.arange(len(mylist))        
        samples = []
        for k in range(sample_number):
            sample_idx = np.random.choice(all_idx, size = sample_size)
            samples.append([mylist[i] for i in sample_idx])               
        return samples

dataset_samples = {}
for name, entries in dataset_entries.items():        
    dataset_samples[name] = sample(entries, sample_size, sample_number)


dataset_data = {}
for name, boot_samples in dataset_samples.items():        
    all_boot_data = []
    for sample in boot_samples:
        boot_data = []
        # get moon_data froom oops dict
        for protein in sample:
            if protein in oops_dict.keys():
                oops_boot_data = list(oops_dict[protein])
            else:
                oops_boot_data = [0, 0, 0]
            if protein in capture_dict.keys():
                capture_boot_data = list(capture_dict[protein])
            else:
                capture_boot_data = [0, 0]
            boot_data.append(oops_boot_data + capture_boot_data)    
        all_boot_data.append(np.array(boot_data))
    dataset_data[name] = all_boot_data


### enrichment ################################################################

# make plot
n = len(fields)
ncolumns = n
nrows = int(np.ceil(n/ncolumns))
myplot = plt.figure(figsize=(5 * ncolumns, 5 * nrows))


for i, field in enumerate(fields):
    field_moon_data = moon_data[:, i]
    # for oligo capture do not consider the number of replicates
    moon_value = np.sum(field_moon_data > 0) / sample_size * 100
    objects = []
    for name, all_boot_data in dataset_data.items():        
        boot_values = []
        for boot_data in all_boot_data:
            field_boot_data = boot_data[:, i]
            # for oligo capture do not consider the number of replicates
            boot_value = np.sum(field_boot_data > 0) / sample_size * 100
            boot_values.append(boot_value)
        objects.append(pd.DataFrame.from_dict({'value': boot_values, 'dataset': name}))
    
    # make dataframe of dataset values for plotting
    df = pd.concat(axis=0, ignore_index=True, objs=objects)
    
    # make subplots   
    plt.subplot(nrows, ncolumns, i + 1)
    sns.histplot(data=df, x="value", hue="dataset", element="step")
    #sns.displot(data=df, x="value", hue="dataset", kind="kde")
    plt.title(field)
    plt.xlabel('% of enriched proteins')
    plt.ylabel('number of bootstrap samples')
    plt.axvline(x=moon_value, color='r')
    
    

myplot.tight_layout(pad=2.0)    
plt.savefig(os.path.join(plotdir, enrichment))

