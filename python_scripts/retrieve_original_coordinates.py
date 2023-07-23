# # retrieve_original_coordinates.py

# ## Usage
#
# retrieve_original_coordinates.py coord_dir padded_fasta_dir
#
# ## Goal
#
# CNEfinder was run on FASTA file with concatenated scaffolds (with NNN... in between)
# This script converts CNE coordinates back to coordinates that correspond to the original FASTA file.
#
# ## Input
#
# coord_dir: directory containing CNE coordinates files.
# Each species has its own coordinate file, named species_id_cne_coords.tsv (e.g. aaur_cne_coords.tsv).
#
# padded_fasta_dir: directory containing padded single-scaffold FASTA files.
# files must be named species_id_pad.fa
#
#
# ## Output
#
# One coordinate file for each species, named 'species_prefix_orig_coords.tsv'
#
#
# ***

import pandas as pd
from Bio import SeqIO
from collections import OrderedDict
import glob
import sys


# #### User input
coord_dir = sys.argv[1]
padded_fasta_dir = sys.argv[2]
out_dir = sys.argv[3]

# #### Increase recursion limit
sys.setrecursionlimit(10**9)


# #### List coordinate files
coord_files = [f for f in glob.glob(coord_dir + "*.tsv")]
print("Found ", len(coord_files), " files in : ", coord_dir)


# #### Scaffold finder function
# Add length of next scaffolds until it finds the scaffold where the CNE is located
def find_scaffold(start_coord, end_coord, cumul_coord, current_scaffold, scaffold_iter):
    cumul_coord += scaffold_lengths[current_scaffold]
    if end_coord < cumul_coord:
        scaff_coord_end = scaffold_lengths[current_scaffold] - (cumul_coord - end_coord)
        scaff_coord_start = scaffold_lengths[current_scaffold] - (cumul_coord - start_coord)
        scaff_coords = [scaff_coord_start, scaff_coord_end]
        return(current_scaffold, scaff_coords)
    else:
        # print("run next iteration")
        current_scaffold = next(scaffold_iter)
        return(find_scaffold(start_coord, end_coord, cumul_coord, current_scaffold, scaffold_iter))


# #### Scaffold length dictionary generator
def create_scaffold_length_dict(fasta_file):
    scaffold_lengths = OrderedDict([])
    print("Retrieving coordinates from original (padded) fasta file: ", fasta_file)
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        scaffold_id = seq_record.id
        scaffold_length = len(seq_record)
        scaffold_lengths[scaffold_id] = scaffold_length
    return(scaffold_lengths)


# #### Main function
def retrieve_original_coordinates(coord_df, scaffold_length_dict):
    # Create output data frame
    orig_coordinates = pd.DataFrame(columns= ['single_sc_start', 'single_sc_end', 'scaffold', 'orig_start', 'orig_end'])
    for index, row in coord_df.iterrows():
        scaffold_iter = iter(scaffold_length_dict)
        cumul_coord = 0
        current_scaffold = next(scaffold_iter)
        start = row['start']
        end = row['end']
        cne_id = row['cne_id']
        #print(cne_id)
        #print("run find_scaffold")
        scaffold = find_scaffold(start, end, cumul_coord, current_scaffold, scaffold_iter)
        scaffold_id = scaffold[0]
        orig_start = scaffold[1][0]
        orig_end = scaffold[1][1]
        orig_coordinates.loc[cne_id] = [start, end, scaffold_id, orig_start, orig_end]
    orig_coordinates.index.name = 'cne_id'
    orig_coordinates.reset_index(inplace=True)
    return(orig_coordinates)


# #### Process each coordinate file

for file in coord_files:
    print("Processing: ", file)
    # Retrieve species name
    species_prefix = file.split("/")[-1].split("_")[0]
    # Retrieve expected fasta file name
    fasta_file = padded_fasta_dir + species_prefix + "_pad.fa"
    # Create output file name
    output_file_name = out_dir + species_prefix + "_orig_coords.tsv"
    # Read coordinate_file
    coord_df = pd.read_csv(file, sep= "\t", names = ['cne_id', 'start', 'end'])
    #coord_df = coord_df.loc[1:]
    # Create ordered dict to hold the scaffold lengths
    scaffold_lengths = create_scaffold_length_dict(fasta_file)
    orig_coordinates = retrieve_original_coordinates(coord_df, scaffold_lengths)
    print("Writing original coordinates file to: ", output_file_name)
    orig_coordinates.to_csv(output_file_name, sep="\t", index=False)

