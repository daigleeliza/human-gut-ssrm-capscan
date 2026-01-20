"""
removes reference sequences from MSA (FASTA format) of hit and reference sequences
"""

import sys
from Bio import SeqIO

msa_withRef = sys.argv[1]
msa_refOnly = sys.argv[2]
msa_noRef = sys.argv[3]

msa_withRef_records = SeqIO.parse(msa_withRef, format='fasta')
ref_records = SeqIO.parse(msa_refOnly, format='fasta')

ref_IDs = []

for ref_record in ref_records:
    ref_IDs.append(ref_record.id)

hit_records = []
for record in msa_withRef_records:
    if record.id not in ref_IDs:
        hit_records.append(record)

SeqIO.write(hit_records, handle=msa_noRef, format='fasta')
