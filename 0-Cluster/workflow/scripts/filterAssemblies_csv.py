import sys

queryIDs = sys.argv[1]
input_bact = sys.argv[2]
input_arch = sys.argv[3]
output_bact = sys.argv[4]
output_arch = sys.argv[5]

with open(queryIDs) as f:
    ids = set(line.strip() for line in f)

for input_csv, output_csv in [(input_bact, output_bact), (input_arch, output_arch)]:
    with open(input_csv) as infile, open(output_csv, "w") as outfile:
        for line in infile:
            queryID = line.split(".")[0]
            if queryID in ids:
                outfile.write(line)