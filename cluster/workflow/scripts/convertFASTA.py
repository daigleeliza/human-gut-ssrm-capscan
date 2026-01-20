"""
converts HMMER hit output file from Stockholm format to FASTA
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

output_msa = sys.argv[1]
FASTA = sys.argv[2]

records = []
msa = SeqIO.parse(output_msa, "stockholm")
for hit in msa:
    ungapped_seq = str(hit.seq.ungap("-"))
    ungapped_seq_upper = ungapped_seq.upper()
    record = SeqRecord(
        Seq(ungapped_seq_upper),
        id=hit.id,
        description=""
        )
    records.append(record)
SeqIO.write(records, FASTA, "fasta")
