""""
Identifies duplicate sequences (exact and shorter duplicates) between Anantharaman2018's dsrAB fasta file and
Mueller2015's dsrAB reference file.
Creates new FASTA file with all duplicates (shorter and exact) removed from Anantharaman2018 sequences
Creates json file with duplicate (exact and shorter) information

INPUT:
    hit_seqs: FASTA file containing duplicate sequences
OUTPUT:
    noEDups_fasta: FASTA file without exact duplicates
    noDups_fasta: FASTA file without exact duplicates or shorter duplicates
    noDups_json: JSON file with duplicate information
"""

import sys
from Bio import SeqIO
import os
import json

Anantharaman2018_seqs_fasta = sys.argv[1]  # Anantharaman2018's sequences in FASTA format with duplicate sequences
Mueller2015_seqs_fasta = sys.argv[2]  # Mueller2015's sequences which we are using as reference sequences
outDir = sys.argv[3]  # output directory

# create output file paths
basename = os.path.splitext(os.path.basename(Anantharaman2018_seqs_fasta))[0]
noEDups_fasta = outDir + "/" + basename + "_noEDups.faa"
noDups_fasta = outDir + "/" + basename + "_noDups.faa"
noDups_json = outDir + "/" + basename + "_noDups.json"

###  REMOVE EXACT DUPLICATES ###

# Create our hash table of Mueller2015's sequences
Mueller2015_sequences = {}

# Make dict for Mueller2015's sequences
for seq_record in SeqIO.parse(Mueller2015_seqs_fasta, "fasta"):
    sequence = str(seq_record.seq).upper()
    if sequence not in Mueller2015_sequences:
        Mueller2015_sequences[sequence] = seq_record.id

# Create hash table to which unique (non-duplicated) sequences from Anantharaman2018 will be added
unique_Anantharaman2018_sequences = {}

# Using the Biopython fasta parse we can read our fasta input
for seq_record in SeqIO.parse(Anantharaman2018_seqs_fasta, "fasta"):
    # Take the current sequence
    sequence = str(seq_record.seq).upper()
    # If the sequence passed in the test "is it clean?" and it isn't in the
    # hash table, the sequence and its id are going to be in the hash
    if sequence not in Mueller2015_sequences:
        unique_Anantharaman2018_sequences[sequence] = seq_record.id
    # If it is already in the hash table, we're just gonna concatenate the ID
    # of the current sequence to another one that is already in the hash table
    else:
        unique_Anantharaman2018_sequences[sequence] += "," + seq_record.id

    # # Write the clean sequences

    # Create a file in the same directory where you ran this script
    with open(noEDups_fasta, "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in unique_Anantharaman2018_sequences:
            output_file.write(">" + unique_Anantharaman2018_sequences[sequence] + "\n" + sequence + "\n")

## REMOVE SHORT DUPLICATES ###

noDups_tempInfo_dict = {}
noEDups_dict = SeqIO.index(noEDups_fasta, "fasta")
Mueller2015_seqs_dict = SeqIO.index(Mueller2015_seqs_fasta, "fasta")

# Create duplicate info dict " noDups_info_dict

# iterate and compare through entire noEDups_dict
for record in noEDups_dict:
    has_dups = False
    seq = noEDups_dict[record].seq

    # check for shorter duplicates in Anantharaman2018's sequences
    for comparison_record in noEDups_dict:
        comparison_seq = noEDups_dict[comparison_record].seq
        if record != comparison_record:
            if seq in comparison_seq:
                has_dups = True
                if comparison_record not in noDups_tempInfo_dict:
                    noDups_tempInfo_dict[comparison_record] = record
                elif noDups_tempInfo_dict[comparison_record] != "":
                    noDups_tempInfo_dict[comparison_record] += ";" + record
                else:
                    noDups_tempInfo_dict[comparison_record] = record

    # check for shorter duplicates in Mueller2015's sequences
    for comparison_record in Mueller2015_seqs_dict:
        comparison_seq = Mueller2015_seqs_dict[comparison_record].seq
        if record != comparison_record:
            if seq in comparison_seq:
                has_dups = True
                if comparison_record not in noDups_tempInfo_dict:
                    noDups_tempInfo_dict[comparison_record] = record
                elif noDups_tempInfo_dict[comparison_record] != "":
                    noDups_tempInfo_dict[comparison_record] += ";" + record
                else:
                    noDups_tempInfo_dict[comparison_record] = record
    if not has_dups:
        if record not in noDups_tempInfo_dict:
            noDups_tempInfo_dict[record] = ""

# remove records that appear as values for other record keys
dups = list(noDups_tempInfo_dict.values())
dups_list = []
for item in dups:
    dups_list += item.split(";")

noDups_info_dict = {}

for key in noDups_tempInfo_dict.keys():
    if key not in dups_list:
        noDups_info_dict[key] = noDups_tempInfo_dict[key]

# write json file for noDups duplicate info
with open(noDups_json, "w+") as output_file:
    json.dump(noDups_info_dict, output_file, indent=4)

# write fasta file for noDups sequences
reference_dict = SeqIO.index(noEDups_fasta, "fasta")

with open(noDups_fasta, "w+") as output_file:
    # Just read the hash table and write on the file as a fasta format
    for record in noDups_info_dict:
        if record in reference_dict:
            sequence = str(reference_dict[record].seq)
            output_file.write(">" + record + "\n" + sequence + "\n")

