
import sys
from pathlib import Path
import math
from collections import Counter

file1 = sys.argv[1]
file2 = sys.argv[2]

print("Creating directories for split CNEFinder files and overlap files")

Path("./split_cf_files").mkdir(parents=True, exist_ok=True)
Path("./overlap_files").mkdir(parents=True, exist_ok=True)


# #### Number of lines per file (recommended=10000)
lines_per_file = 10000
print("Number of lines per file: ", lines_per_file)

# #### Split cnefinder input files
def split_cf_file(cf_file):
    smallfile = None
    filenames = []
    with open(cf_file) as big_file:
        counter = 0
        for lineno, line in enumerate(big_file):
            if lineno % lines_per_file == 0:
                counter+=1
                if smallfile:
                    smallfile.close()
                small_filename = "split_cf_files/" + cf_file.split("/")[-1] + '_{}'.format(counter)
                filenames.append(small_filename)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    return filenames

file_names_1 = split_cf_file(file1)
file_names_2 = split_cf_file(file2)
combs = [(x,y) for x in file_names_1 for y in file_names_2]

comp1 = file1.split("/")[-1].split(".out")[0].split("_vs_")
comp2 = file2.split("/")[-1].split(".out")[0].split("_vs_")
species_list = comp1 + comp2
specific_species = []
for sp, count in Counter(species_list).items():
    if count == 2:
        common_species = sp
    elif count == 1:
        specific_species.append(sp)
swarm_file = common_species + "_overlap_" + specific_species[0] + "_" + specific_species[1] + ".swarm"
print("Writing swarm file:", swarm_file)

with open(swarm_file, 'a') as output_file:
    for comparison in combs:
        output_file.write('python calculate_overlaps_split.py ' +  comparison[0] + " " + comparison[1] + " overlap_files/ ;\n")

print("Done. Make sure calculate_overlaps_split.py is in the current directory.")
print("Run swarm with 2G")

