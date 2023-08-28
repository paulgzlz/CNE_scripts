
import json
import pandas as pd
import  csv
import sys
from collections import defaultdict
import glob


# #### User input

cne_dict = sys.argv[1]
cf_output_dir = sys.argv[2]


# #### Test input paths
#cne_dict = '../../results_for_paper/cnidaria_final/filtering/filtered_cne_dict.txt'
#cf_output_dir = '../../results_for_paper/cnidaria_final/cnefinder_output/stranded/'


# #### List CNEFinder output files
cf_output_files = glob.glob(cf_output_dir + '*.out')

# #### Open dictionary of non_overlapping cnes created with generate_cne_ids.py

print("reading dictionary of CNEs")
with open(cne_dict) as json_file:
    all_species_cne_dict = json.load(json_file)


# #### Function that retrieves CNE_ids from each 2-species CNE
# Returns list of list:  
# [[sp1_cne_X, sp2_cne_X], [sp1_cne_X, sp2_cne_X] ... ]

def retrieve_pairwise_links(cne_file):
    print("cne_file is:", cne_file)
    output_list = []
    # Identiy reference and query species
    #ref_species = cne_file.split("/")[-1].split("_vs_")[0]
    #query_species = cne_file.split("/")[-1].split("_vs_")[1].split(".")[0]
    # Read CNEFinder file 
    #column_names = ['ref_chrom', 'ref_start', 'ref_end',
    #                  'query_chrom', 'query_start', 'query_end',
    #                  'ref_length', 'query_length', 'sim']
    
    cnes = pd.read_csv(cne_file, sep = "\t")#, names = column_names)
    ref_species = cnes.columns[0].split("_")[0]
    query_species = cnes.columns[3].split("_")[0]    

    # Iterate over each row (=CNE) in CNEFidner output file
    for index, row in cnes.iterrows():
        # Create list that will contain the two CNE IDs in pair
        pairwise_list = []
        #ref_start = row['ref_start']
        #ref_end = row['ref_end']
        #query_start = row['query_start']
        #query_end = row['query_end']
        ref_start = row[1]
        ref_end = row[2]
        query_start = row[4]
        query_end = row[5]        
        # Retrieve corresponding CNE_ids by searching overlapping coordinates in unique_non_overlap_cnes.txt
        # unique_non_overlap_cnes.txt was loaded as dictionary (all_species_cne_dict)
        for cne_id , cne_id_coords in all_species_cne_dict[ref_species].items():
            if (cne_id_coords[0] <= ref_end <= cne_id_coords[1]) or  (ref_start <= cne_id_coords[1] <= ref_end):
                pairwise_list.append(cne_id)
                break
        for cne_id , cne_id_coords in all_species_cne_dict[query_species].items():
            if (cne_id_coords[0] <= query_end <= cne_id_coords[1]) or  (query_start <= cne_id_coords[1] <= query_end):
                pairwise_list.append(cne_id)
                break
        # Add list to final list of lists. 
        # If a CNE in pair was filtered out (absent from unique_non_overlap_cnes), discard pair
        if len(pairwise_list) == 2 and pairwise_list not in output_list:
            output_list.append(pairwise_list)
    return(output_list)


# #### Create list of links and add all pairwise links

print("Identifying all pairwise links")
all_pairwise_links = []
for cf_file in cf_output_files:
    pairwise_links = retrieve_pairwise_links(cf_file)
    all_pairwise_links = all_pairwise_links + pairwise_links
print("Done")


with open('pairwise_links.json', 'w') as f:
    json.dump(all_pairwise_links, f)


