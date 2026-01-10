
from Bio import SeqIO
import sys

Anantharaman2018_seqs_file = sys.argv[1]
output_file = sys.argv[2]

query_id_list = []
for record in SeqIO.parse(Anantharaman2018_seqs_file, "fasta"):
    query_id = "QUERY___" + str(record.id)
    query_id_list.append(query_id)

with open(output_file, 'w') as f:
    for item in query_id_list:
        f.write("%s\n" % item)

