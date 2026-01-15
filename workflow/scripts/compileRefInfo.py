""""
This code adds closest reference leaf (by branch distance) and sample ID to all hit sequences
from dsrAB HMMER search that passed the score threshold (so, all duplicate sequences).

INPUT:
    queryDistanceInfoCSV :
        CSV file with hit sequence info for unique hit sequences only (no duplicates)
    biosampleJSON :
        key = sample_ID, value "query_IDs" = WGS accession number
    scoreThresholdCSV :
        CSV file with hit sequences for all hits that pass the score threshold

OUTPUT:
    closestRefInfo_allScoreThresholdHits.csv
"""

import csv
import json
import sys
import ast

queryDistanceCSV = sys.argv[1]
scoreThresholdCSV = sys.argv[2]
biosampleJSON = sys.argv[3]
outputCSV = sys.argv[4]

# iterate through scoreThreshold CSV and add reference info and sample info for each row

with open(scoreThresholdCSV, 'r') as f:
    scoreThreshold_reader = csv.reader(f)

    with open(biosampleJSON, 'r') as f:
        biosample_dict = json.load(f)

        with open(outputCSV, 'w') as f:
            writer = csv.writer(f)

            writer.writerow(["hit_id",
                             "query_id",
                             "biosample_id",
                             "closest_ref",
                             "distance",
                             "SRM_present"])

            next(scoreThreshold_reader)  # skip header
            for row in scoreThreshold_reader:
                hit_id = row[0]  # first column in csv file. hit_id includes full sequence info
                query_id = hit_id.split(".")[0]  # query_id is the WGS accession prefix

                # find biosample_id for this query_id from dict
                for biosample_id, metadata in biosample_dict.items():
                    if query_id in metadata["query_ids"]:
                        sample_id = biosample_id
                        break  # exit for loop once you've found the sample_id

                    # find closest_ref and distance from query_distance_dict
                    # closest_ref = []
                    # for hit, dist_info in queryDistance_dict.items():
                    #     if hit_id in hit:
                    #         closest_ref = dist_info["closest_ref"]
                    #         distance = dist_info["dist"]
                    #     break
                with open(queryDistanceCSV, 'r') as f:
                    queryDistance_reader = csv.DictReader(f)

                    closest_ref = []
                    distance = None
                    SRM_present = False
                    for queryDist_row in queryDistance_reader:
                        if hit_id in queryDist_row["id"]:
                            # closest_ref = ast.literal_eval(queryDist_row["closest_ref"])
                            SRM_present = True
                            distance = queryDist_row["dist"]
                            if type(queryDist_row["closest_ref"]) is str:
                                closest_ref = ast.literal_eval(queryDist_row["closest_ref"])
                            else:
                                closest_ref = queryDist_row["closest_ref"]
                            break
                        if hit_id in queryDist_row["edups"] and hit_id not in queryDist_row["id"]:
                            SRM_present = True
                            # closest_ref = ast.literal_eval(queryDist_row["closest_ref"])
                            if type(queryDist_row["closest_ref"]) is str:
                                closest_ref = ast.literal_eval(queryDist_row["closest_ref"])
                            else:
                                closest_ref = queryDist_row["closest_ref"]
                            distance = queryDist_row["dist"]
                            break
                        if hit_id in queryDist_row["sdups"]:
                            SRM_present = True
                            # closest_ref.append(ref)  # make list of closest_ref if shorter dup
                            if type(queryDist_row["closest_ref"]) is str:
                                for ref in ast.literal_eval(queryDist_row["closest_ref"]):
                                    if ref not in closest_ref:
                                        closest_ref.append(ref)  # make list of closest_ref if shorter dup
                            else:
                                for ref in queryDist_row["closest_ref"]:
                                    if ref not in closest_ref:
                                        closest_ref.append(ref)
                            if not distance:
                                distance = queryDist_row["dist"]
                            else:
                                distance = "NA"

                    writer_row = [hit_id, query_id, sample_id, closest_ref, distance, SRM_present]
                    writer.writerow(writer_row)
