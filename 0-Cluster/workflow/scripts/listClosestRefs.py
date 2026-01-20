""""
This code creates a list of all the closest reference leaves represented among the query sequences.
"""

import json
import sys

queryDistanceInfo = sys.argv[1]
Karthik_novel_seq_file = sys.argv[2]
output_ref_list_f = sys.argv[3]

with open(queryDistanceInfo, 'r') as f:
    dist_dict = json.load(f)

    closest_ref_list = []

    for query, info in dist_dict.items():
        for ref in info["closest_ref"]:
            if ref not in closest_ref_list:
                closest_ref_list.append(ref)

# write reference leaf list
with open(output_ref_list_f, 'w') as f:
    for item in closest_ref_list:
        f.write("%s\n" % item)


