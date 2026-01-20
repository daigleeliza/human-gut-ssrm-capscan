"""
creates new FASTA and CSV files with only hits that have passed the score threshold
    INPUT:
        all_hits_fasta: FASTA of all compiled hits
        all_hits_csv: CSV of all compiled HMMER data

    OUTPUT:
        score_threshold_fasta: FASTA of hits with score >= 100
        score_threshold_csv: CSV of hits with score >= 100
"""

import csv
from Bio import SeqIO
import sys

all_hits_fasta = sys.argv[1]
all_hits_csv = sys.argv[2]
score_threshold_fasta = sys.argv[3]
score_threshold_csv = sys.argv[4]
score_threshold = float(sys.argv[5])

# create counters for number of total hits and number of hits that meet score threshold
n_allHits = 0
n_scoreThresholdHits = 0

scoreThreshold_ID_list = []

# write new CSV file for sequences that meet score threshold
with open(all_hits_csv, mode='r') as f_read:
    with open(score_threshold_csv, mode='w') as f_write:
        csv_reader = csv.reader(f_read)
        csv_writer = csv.writer(f_write)
        for row in csv_reader:
            n_allHits += 1
            score = float(row[7])
            if score >= score_threshold:
                n_scoreThresholdHits += 1
                csv_writer.writerow(row)
                hitID = row[0]
                scoreThreshold_ID_list.append(hitID)


# write new FASTA file with hits that meet score threshold
scoreThreshold_records = []
allHit_records = SeqIO.parse(all_hits_fasta, format="fasta")

for record in allHit_records:
    if record.id in scoreThreshold_ID_list:
        scoreThreshold_records.append(record)

SeqIO.write(scoreThreshold_records, score_threshold_fasta, format="fasta")

