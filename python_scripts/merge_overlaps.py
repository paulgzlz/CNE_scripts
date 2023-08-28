import glob
import json
import sys

overlap_dir = sys.argv[1]
output_file_name = sys.argv[2]

#overlap_dir = "adig_spis_pdam_spis_split/"
#output_file_name = "spis_overlap_adig_pdam.txt"


# #### List dictionaries

overlap_files = [f for f in glob.glob(overlap_dir + "*.txt")]
print("Found ", len(overlap_files), " files in : ", overlap_dir)

counter = 1
all_overlaps = {}
for overlap_file in overlap_files:
    with open(overlap_file) as json_file:
        overlap_dict = json.load(json_file)
    print(len(overlap_dict))
    for cne_id, coord_dict in overlap_dict.items():
        new_id = cne_id.split("_cne_")[0] + "_cne_" + str(counter)
        all_overlaps[new_id] = coord_dict
        counter += 1

print("writing output file:", output_file_name)
with open(output_file_name, 'w') as outfile:
    json.dump(all_overlaps, outfile)
print("Done, bye.")


