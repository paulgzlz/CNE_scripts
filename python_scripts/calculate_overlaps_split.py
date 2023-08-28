# # calculate_overlaps_split.py

# ## Usage
# 
# calculate_overlaps.py cnefinder_output_file_1 cnefinder_output_file_2 output_dir
# 
# ## Goal
# 
# Identify overlapping CNEs across pairwise comparisons that have a species in common.  
# Extend coordinates of each CNE until it includes all overlapping CNEs.
# 
# ## Input
# 
# Two CNEFinder output files
# 
# ## Output
# 
# overlap file (json dictionary) named common_species_vs_specific_species_1_specific_species_2.txt  
# For example, run with: 'aaur_vs_epal.out' and 'aaur_vs_hsym.out' will output: 'aaur_vs_epal_hsym.txt 
# 
# ***

import pandas as pd
from Bio import SeqIO
import sys
import json
from collections import Counter
from os import listdir
from os.path import isfile, join
import itertools

# #### User input
file1 = sys.argv[1]
file2 = sys.argv[2]
output_dir = sys.argv[3]

# #### Identify common and specific species names

comp1 = file1.split("/")[-1].split(".out")[0].split("_vs_")
comp2 = file2.split("/")[-1].split(".out")[0].split("_vs_")
for species in comp1:
    if species in comp2:
        common_species = species
    else:
        specific_comp1 = species
for species in comp2:
    if species not in comp1:
        specific_comp2 = species


# #### Identify query and reference species in both files (=comparisons)
ref_comp1 = comp1[0]
query_comp1 = comp1[1]
ref_comp2 = comp2[0]
query_comp2 = comp2[1]

# #### Create column names with correct species names for reading CNEFinder files

column_names_comp1 = [ref_comp1 + '_chrom', ref_comp1 +'_start', ref_comp1 + '_end', ref_comp1 + 'strand',
                  query_comp1 + '_chrom', query_comp1 + '_start', query_comp1 + '_end', query_comp1 + 'strand',
                  ref_comp1 + '_length', query_comp1 + '_length', 'sim']
column_names_comp2 = [ref_comp2 + '_chrom', ref_comp2 +'_start', ref_comp2 + '_end', ref_comp2 + 'strand',
                  query_comp2 + '_chrom', query_comp2 + '_start', query_comp2 + '_end', query_comp2 + 'strand',
                  ref_comp2 + '_length', query_comp2 + '_length', 'sim']


# #### Read CNEFinder files

comp1_cnes = pd.read_csv(file1, sep = "\t", names = column_names_comp1, header=0)
comp2_cnes = pd.read_csv(file2, sep = "\t", names = column_names_comp2, header=0)


# #### Create variables for searching CNEFinder file

common_start = common_species + '_start'
common_end = common_species + '_end'
spec_comp1_start = specific_comp1 + '_start'
spec_comp1_end = specific_comp1 + '_end'
spec_comp2_start = specific_comp2 + '_start'
spec_comp2_end = specific_comp2 + '_end'


# #### Create dictionary of common CNEs

def check_overlap(coord_pair1, coord_pair2):
    overlap = False
    start1 = coord_pair1[0]
    end1 = coord_pair1[1]
    start2 = coord_pair2[0]
    end2 = coord_pair2[1]
    if (start1 <= end2 <= end1) or (start1 <= start2 <= end1) or (start2 <= start1 <= end2) :
        overlap = True
    return overlap


common_cnes = {}
i = 1
# Iterate over rows of CNEFinder file 1
for index, row in comp1_cnes.iterrows():
    # Retrieve CNE coords for common species and species specific to file 1
    common_start_comp1 = row[common_start]
    common_end_comp1 = row[common_end]
    spec_comp1_start_coord = row[spec_comp1_start]
    spec_comp1_end_coord = row[spec_comp1_end]
    coord_spec_comp1 = str(spec_comp1_start_coord) + ":" + str(spec_comp1_end_coord)
    # Iterate over all rows of CNEFinder file 2 to look for overlapping CNEs
    for index, row in comp2_cnes.iterrows():
        # Retrieve CNE coords for common species and species specific to file 2
        common_start_comp2 = row[common_start]
        common_end_comp2 = row[common_end]
        spec_comp2_start_coord = row[spec_comp2_start]
        spec_comp2_end_coord = row[spec_comp2_end]
        coord_spec_comp2 = str(spec_comp2_start_coord) + ":" + str(spec_comp2_end_coord)
        # Test if CNE overlaps
        if check_overlap([common_start_comp1, common_end_comp1],[common_start_comp2, common_end_comp2] ):
   # Give a new CNE_ID to CNE and store coordinates in dictionary
            cne_id = common_species + "_cne_" + str(i)
            coord_common = str(min(common_start_comp1, common_start_comp2)) + ":"                 + str(max(common_end_comp1, common_end_comp2))
            common_cnes[cne_id] = {common_species: coord_common,
                                     specific_comp1: coord_spec_comp1,
                                     specific_comp2: coord_spec_comp2}
            i = i + 1


# #### Write dictionary to file

output_file_name = output_dir + file1.split("/")[-1] + "_overlap_" + specific_comp1 + "_" + file2.split("/")[-1] + ".txt"
with open(output_file_name, 'w') as file:
    file.write(json.dumps(common_cnes)) # use `json.loads` to do the reverse


