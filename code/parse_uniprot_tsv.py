# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 17:44:40 2021

Project: Moonlighting Proteins, Prof Lilley's Lab, Amgen Scolar Program

Author: Lorenzo Magni

Contact: lorenzo.magni@sns.it

"""

import numpy as np   

class parse_uniprot_tsv:
    
    def __init__(self, path_to_file):
        #initialise attributes
        # dictionary mapping attributes of data to index (to speed up retrival of specific attributes)
        self.fields = {}
        # dictionary mapping entry to index (to speed up retrival of individual entries)
        self.entries = {}
        # list of lists corresponding to each line in the file
        self.data = []        
        # parse file
        with open(path_to_file, 'r') as file:
            lines = file.readlines()
            fields = lines[0].replace("\n", "").split("\t")
            for i, field in enumerate(fields):
                self.fields[field] = i
            for j, line in enumerate(lines[1:]):
                entry_ID = line.split("\t")[0]
                self.entries[entry_ID] = j          
                self.data.append(line.split("\t"))
        return
    
    # trivial functions to get list of fields or entries
    def list_fields(self):
        return list(self.fields.keys())    
    def list_entries(self):
        return list(self.entries.keys())
    
    
    # find specified attributes for specified entries
    # by default gives all attributes and all entries                
    def find(self, *entries, **kwargs):
       
        # retrieve indices of required entries
        entry_idx = []
        if entries:
            for entry_ID in entries:
                entry_idx.append(self.entries[entry_ID])
        else:
            entry_idx = list(self.entries.values())
        
        # retrieve required attributes
        try:
            fields = kwargs["fields"]
        except KeyError:
            # by default output all fields
            fields = self.fields.keys()
        
        # retrieve indices of required fields
        key_idx =[]
        for field in fields:
            if field in self.fields.keys():
                key_idx.append(self.fields[field])
            
        # return required fields
        out_data = []
        for i in entry_idx:
            out_data.append([self.data[i][j] for j in key_idx])
        
        return out_data

    # returns single required field (kwarg field) about specific entries
    # by default returns ids
    # by default considers all entries
    def find_field(self, *entries, **kwargs):
       
        # retrieve indices of required entries
        entry_idx = []
        if entries:
            for entry_ID in entries:
                entry_idx.append(self.entries[entry_ID])
        else:
            entry_idx = list(self.entries.values())
        
        # retrieve required field
        try:
            field = kwargs["field"]
            # retrieve index of required fields
            key_idx = self.fields[field]
        except KeyError:
            # by default output all fields
            key_idx = 0
               
            
        # return required fields
        out_data = []
        for i in entry_idx:
            out_data.append(self.data[i][key_idx])
        
        return out_data    

    # returns all fields if specified attributes are not empty
    def isit(self, *fields):
        field_idx = [self.fields[field] for field in fields]
        itis = []
        for element in self.data:
            for idx in field_idx:
                if element[idx]:
                    itis.append(element)        
        return itis
    
    # returns entry if specified fields are not empty
    def which_entries(self, *fields):
        field_idx = [self.fields[field] for field in fields]
        itis = []
        for element in self.data:
            for idx in field_idx:
                if element[idx]:
                    itis.append(element[0])
        return itis
        
              
    # output specified fields of a random sample of a given size
    def sample(self, sample_size, **kwargs):
        all_idx = np.arange(len(self.entries))
        sample_idx = np.random.choice(all_idx, size = sample_size)
        
        # retrieve required fields
        try:
            fields = kwargs["fields"]
        except KeyError:
            # by default output all fields
            fields = self.fields.keys()
        
        # retrieve indices of required attributes
        key_idx =[]
        for field in fields:
            if field in self.fields.keys():
                key_idx.append(self.fields[field])
            
        # return required fields
        out_data = []
        for i in sample_idx:
            out_data.append([self.data[i][j] for j in key_idx])
        
        return out_data

    # samples a specific single field
    def sample_field(self, sample_size, **kwargs):        
        all_idx = np.arange(len(self.entries))
        sample_idx = np.random.choice(all_idx, size = sample_size)    
        
        # retrieve required field
        try:
            field = kwargs["field"]
            # retrieve index of required fields
            key_idx = self.fields[field]
        except KeyError:
            # by default output entries
            key_idx = 0
               
            
        # return required fields
        out_data = []
        for i in sample_idx:
            out_data.append(self.data[i][key_idx])
        
        return out_data
            
            
