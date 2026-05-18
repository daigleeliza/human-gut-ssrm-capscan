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
# if this is not working, recommended to use what was used for asr portion 
#esl-reformat -u fasta {input.A} |
#		awk -v coassembly="{wildcards.coassembly}" '
#			/^>/ {{
#				sub(/\\s*\\[.*$/, "", $0)
#				contig = substr($0, 2)
#				print ">" coassembly "." contig
#				next
#			}}
#			{{print}}
#		' > {output.A}
