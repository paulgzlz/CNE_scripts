import json
import csv
import glob
from collections import defaultdict

pairwise_file = 'pairwise_links.json'
threeway_dir = "threeway_links/"

threeway_files = glob.glob(threeway_dir + "*json")

# #### Function that merges all sublists having common elements.
def merge_common(lists):
    neigh = defaultdict(set)
    visited = set()
    for each in lists:
        for item in each:
            neigh[item].update(each)
    def comp(node, neigh = neigh, visited = visited, vis = visited.add):
        nodes = set([node])
        next_node = nodes.pop
        while nodes:
            node = next_node()
            vis(node)
            nodes |= neigh[node] - visited
            yield node
    for node in neigh:
        if node not in visited:
            yield sorted(comp(node))


# #### Read pairwise_links
print("reading list of pairwise CNEs")
with open(pairwise_file) as json_file:
    pairwise_links = json.load(json_file)


# #### Read all threeway links
print("reading threeway links. Number of files: ", len(threeway_files) )
all_threeway_links = []
counter = 0
for file in threeway_files:
    if counter % 10 == 0:
        print(counter, " files processed.")  
    with open(file) as json_file:
        threeway_links = json.load(json_file)
        all_threeway_links = all_threeway_links + threeway_links
        counter = counter + 1
#### Combine pairwise and threeway links
print("Combining links")
all_cne_links = pairwise_links + all_threeway_links


#### Merge all links
print("Merging homologous clusters, this may take some time")
merged_list = list(merge_common(all_cne_links))
print("Merging done. Number of clusters: ", len(merged_list))


# #### Write to file for downstream analyses
print("Writing CNE clusters to file: merged_cne_clusters.csv")
with open("merged_cne_clusters.csv","w") as f:
    wr = csv.writer(f)
    wr.writerows(merged_list)


