# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 15:20:00 2021

Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it

"""

import parse_uniprot_tsv as prs
import os
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#paths
"""customise paths"""
path_to_uniprot_tsv = "path to uniprot tab separated file"
path_to_moonlighting_entries = "path to file with list of moonlighting protein"
outfile = "path to text output"
plot = "path to output plot"


proteome = prs.parse_uniprot_tsv(path_to_uniprot_tsv)

moolighting_entries = []
with open(path_to_moonlighting_entries, 'r') as file:
    moolighting_entries = [line.replace("\n", "") for line in file.readlines()]

# defined parameters: mean, variance, count, mean count, mean domain count
mydata = [("Mass", "mean"), 
            ("Mass", "variance"), 
            ("Length", "mean"), 
            ("Length", "variance"), 
            ("Cofactor", "count"),
            ("Calcium binding", "count"),
            ("Glycosylation", "count"),
            ("Lipidation", "count"),
            ("Propeptide", "count"),
            ("Nucleotide binding", "count"), 
            ("Post-translational modification", "count"), 
            ("Gene ontology (cellular component)", "mean count"),
            ("Interacts with", "mean count"),
            ("Domain [FT]", "mean domain count"),
            ("PubMed ID", "mean count"),
            ("Involvement in disease", "count"),
            ("Cross-reference (PDB)", "mean count"),
            ("Cross-reference (PDB)", "count")
            ]

#myfields =list(set([x[0] for x in mydata]))
myfields =[x[0] for x in mydata]
print(myfields)
# process data
def get_mean(data):    
    data = [float(x.replace(",","")) for x in data]        
    mean = np.mean(data)
    return mean

def get_variance(data):
    data = [float(x.replace(",","")) for x in data]
    variance = np.var(data)
    return variance

# count number of entries with non-empty field
def get_count(data):        
    count = 0
    for x in data:
        if x:
            count += 1
    return count

# count number of domains for a single entry
def domain_count(domains):
    count = 0
    if domains:
        count = domains.count("DOMAIN")
    return count

# calculate mean number of domains
def mean_domain_count(data):
    counts = [domain_count(x) for x in data]
    mean = np.mean(counts)
    return mean

def mean_count(data):
    counts = [len(x.split(";")) for x in data]
    mean = np.mean(counts)
    return mean

# number of isoforms
def get_number_of_isoforms(field):    
    try:
        i = field.index("Named isoforms=") + len("Named isoforms=")
        n = int(field[i])
    except:
        n = 0    
    return n

def mean_number_of_isoforms(data):
    counts = [get_number_of_isoforms(x) for x in data]
    mean = np.mean(counts)
    return mean

def var_number_of_isoforms(data):
    counts = [get_number_of_isoforms(x) for x in data]
    mean = np.var(counts)
    return mean

def process(data, **kwargs):
    par = kwargs["par"]
    if par == "mean":
        return get_mean(data)
    if par == "variance":
        return get_variance(data)
    if par == "count":
        return get_count(data)    
    if par == "mean count":
        return mean_count(data)
    if par == "mean domain count":
        return mean_domain_count(data)
    if par == "mean number of isoforms":
        return mean_number_of_isoforms(data)
    if par == "variance in number of isoforms":
        return var_number_of_isoforms(data)
   

# create output file
outfile.write("# Bootstrapping analysis of Uniprot reviewed proteins in comparison with moonlighting entries from MoonProt and Multitask\n\n")
   
# create plot
n = len(mydata)
ncolumns = 4
nrows = int(np.ceil(n/ncolumns))
myplot = plt.figure(figsize=(5 * ncolumns, 5 * nrows))



# find data for moonlighting proteins
moon_all_data = proteome.find(*moolighting_entries, fields=myfields)

# bootstrapping
sample_size = len(moolighting_entries)
sample_number = 1000
boot_all_data = []
for k in range(sample_number):
    boot_all_data.append(proteome.sample(sample_size, fields=myfields))
    




for i, (myfield, par) in enumerate(mydata):    
    field_idx = myfields.index(myfield)
    outfile.write("examined field: " + myfield + "\nexamined paramenter: " + par +"\n")
    print("examined field: " + myfield + "\nexamined parameter: " + par)
    
    # restrict data to the specific field
    moon_data = [x[field_idx] for x in moon_all_data]    
        
    # bootstrapping    
    boot_data = []
    for k in range(sample_number):
        boot_data.append([x[field_idx] for x in boot_all_data[k]])
                
    # get info for moonlighting proteins
    moon_value = process(moon_data, par=par)
    print(moon_value)
    outfile.write("value in moonlighting proteins: " + str(moon_value) +"\n")
    
    # get infos about bootstrap data
    boot_values = [process(data, par=par) for data in boot_data]
    outfile.write("mean value across bootstrap samples: " + str(np.mean(boot_values)) +"\n")
    outfile.write("variance across bootstrap samples: " + str(np.var(boot_values)) +"\n\n")
        
    # make subplots   
    plt.subplot(nrows, ncolumns, i + 1)
    sns.histplot(data=boot_values)
    plt.title(myfield)
    plt.xlabel(par)
    plt.ylabel('bootstrap samples')
    plt.axvline(x=moon_value, color='r')

myplot.tight_layout(pad=2.0)    
plt.savefig(plot)

outfile.close()
