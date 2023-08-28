
# coding: utf-8

# # generate_cne_ids.py

# ## Usage
# 
# generate_cne_ids.py cnefinder_out_dir
# 
# ## Goal
# 
# Cluster overlapping CNEs from one species across multiple CNEFinder runs (pairwise comparisons)
# 
# ## Input
# 
# directory containing all CNEfinder output files. 
# 
# ## Output
# 
# 1. Dictionary of CNEs and their coordinates for each species, named 'unique_non_overlap_cnes.txt'  
# {"aaur": {"aaur_cne_0": [0, 0], "aaur_cne_1": [47284721, 47284885], "aaur_cne_2": [216228632, 216228834], etc  
#   
# 2. Coordinates of each non-redundant CNE. One file oper species, names 'species_cne_coords.tsv'. 
# 
# 
# ***

# In[1]:


import pandas as pd
from Bio import SeqIO
import sys
import json
from collections import Counter
from os import listdir
from os.path import isfile, join
import itertools


# #### Increase recursion limit
# For very large numbers of CNEs, increase value if run fails
sys.setrecursionlimit(15000)


# #### User input
working_dir = sys.argv[1]


# #### List CNEFinder output files

cne_files = [f for f in listdir(working_dir) if isfile(join(working_dir, f))]
print("cne_files:", cne_files)


# #### Create empty dictionary of cnes for each species
# 
# This will be populated by calling parse_pairwise_comps

all_species_cne_dict= {}
for file in cne_files:
    species = file.split(".out")[0].split("_vs_")
    for species_name in species:
        if species_name not in all_species_cne_dict:
            first_cne_id = species_name + "_cne_0"
            all_species_cne_dict[species_name] = {first_cne_id: [[0,0]]}


# #### Column names for reading CNEFinder output files
#column_names = ['ref_chrom', 'ref_start', 'ref_end', 'ref_strand', 'query_chrom', 'query_start', 'query_end', 'query_strand',
 #                     'ref_length', 'query_length', 'sim']


# #### Function that parses cnefinder output
# 
# Identify cnes with overlapping coordinates from different cnefinder runs
def parse_pairwise_comps(species, comp):
    for index, row in comp.iterrows():
        ref_query = "ref"
        for species_name in species:
            i = 0
            # Create a dictionary for this species
            species_dict = all_species_cne_dict[species_name]
            # Determine if species was query or reference in cnefinder run
            # Retrieve coordinates of cne for each species in pair
            if ref_query == "ref":
                start = row['ref_start']
                end = row['ref_end']
            elif ref_query == "query":
                start = row['query_start']
                end = row['query_end']
            # Loop through existing cnes for this species
            for cne_id, cne_coord in all_species_cne_dict[species_name].copy().items():
                for coord_set in cne_coord:
                    cne_start = coord_set[0]
                    cne_end = coord_set[1]
                    # Search for overlaps between current cne and existing cnes
                    # if overlap, append coordinates to cne
                    if (cne_start <= end <= cne_end) or (start <= cne_end <= end) :
                        all_species_cne_dict[species_name][cne_id].append([start, end])
                        i = 1
                        break
            # If no overlap was found, create new_cne
            if i == 0:
                new_cne_id = species_name + "_cne_" + str(len(all_species_cne_dict[species_name]))
                all_species_cne_dict[species_name][new_cne_id] = []
                all_species_cne_dict[species_name][new_cne_id].append([start, end])
            ref_query = "query"


# #### Function that extends CNE coordinates
# For each cne with set of overlapping coordinates, create the largest range that includes all coordinates
def extend_cne_coords(species_dict):
    new_dict = {}
    for cne_id, cne_coords in species_dict.items():
        new_start = (min([j for i in cne_coords for j in i]))
        new_end = (max([j for i in cne_coords for j in i]))
        new_dict[cne_id] = [new_start, new_end]
    return(new_dict)


# #### Function that recursively merge overlapping coordinates
def recursive_merge(inter, start_index = 0):
    for i in range(start_index, len(inter) - 1):
        if inter[i][1] > inter[i+1][0]:
            new_start = inter[i][0]
            new_end = inter[i+1][1]
            inter[i] = [new_start, new_end]
            del inter[i+1]
            return recursive_merge(inter.copy(), start_index=i)
    return inter


# #### Parse each cnefinder output file
for file in cne_files:
    print("Processing cne file: ", file)
    species = file.split(".out")[0].split("_vs_")
    comp = pd.read_csv(working_dir + file, sep = "\t")#, names = column_names)
    parse_pairwise_comps(species, comp)
    print("file ", file, "processed")


# #### Create dictionary of all cnes with extended coordinates
print("Extending CNE coordinates")
all_species_cne_dict_merged = {}
for species, species_dict in all_species_cne_dict.items():
    merged_species_dict = extend_cne_coords(species_dict)
    all_species_cne_dict_merged[species] = merged_species_dict


# #### Recursively merge CNEs with new overlaps
# Overlapping cnes were generated sequentially.  
# As a result, cnes that were not overlapping may now overlap.  
# cnes need to be merged recursively to avoid this problem.  

print("Recursively merging CNEs with new overlaps")
unique_non_overlap_cnes = {}
for species, species_dict in all_species_cne_dict_merged.items():
    unique_non_overlap_cnes[species] = {}
    coord_list =  sorted(list(species_dict.values()))
    merged_list = recursive_merge(coord_list)
    for i in range(len(merged_list)):
        cne_id = species + "_cne_" + str(i)
        unique_non_overlap_cnes[species][cne_id] = merged_list[i]


# #### Write dictionary of non-overlapping CNEs to file

print("Write dictionary to file")
# Write dictionary of non-overlapping, unique cnes to file
with open('unique_non_overlap_cnes.txt', 'w') as file:
     file.write(json.dumps(unique_non_overlap_cnes)) # use 'json.loads' to do the reverse


# #### Write cne coordinates to file for downstream analysis

print("Write cne coordinates for downstream analysis")
# Write file of cne coordinates for each species
for species, cnes in unique_non_overlap_cnes.items():
    file_name = species + "_cne_coords.tsv"
    df = pd.DataFrame.from_dict(cnes,'index')
    df.to_csv(file_name, sep='\t', header = False)

print("All done, Bye")


