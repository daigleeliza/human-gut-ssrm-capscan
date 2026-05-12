from Bio import SeqIO
import sys

queryIDs = sys.argv[1]
input_bact = sys.argv[2]
input_arch = sys.argv[3]
output_bact = sys.argv[4]
output_arch = sys.argv[5]

with open(queryIDs) as f:
    ids = set(line.strip() for line in f)

for input_faa, output_faa in [(input_bact, output_bact), (input_arch, output_arch)]:
    with open(input_faa) as infile, open(output_faa, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            queryID = record.id.split(".")[0]
            if queryID in ids:
                SeqIO.write(record, outfile, "fasta")