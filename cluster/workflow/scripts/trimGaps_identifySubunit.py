"""
Function: identify query sequences in alignment as dsrA, dsrB, or both.
Steps:
1. Identify if sequence falls within start/end regions of dsrA or dsrB.
2. If sequence spans both regions, identify of % of gap position in A and B regions.
3. If seq has higher % gap than threshold in ONE region, trim sequence to span only opposite region.
4. If seq has higher % gap than threshold in BOTH regions, remove from MSA.

Input:
	QUERY_MSA: MSA of query sequences aligned to ref. Contains only query. Format = fasta

Output:
	PERC_GAPS_CSV: Dictionary of query id's and percentage of gaps.

	QUERY_IDENTITIES: Dictionary of query and their respective dsrA/B identities.
		{
		QUERY_ID1: "dsrA",
		QUERY_ID2: "dsrB",
		QUERY_ID3: "dsrAB",
		QUERY_ID4: "remove"
		}
	CLEAN_MSA: MSA of cleaned query sequences, trimmed to fit region or removed if gappy.

"""

from Bio import SeqIO
import re
import csv
import json
import sys

QUERY_MSA = sys.argv[1]
DUPS_DICT = sys.argv[2]
QUERY_INFO = sys.argv[3]
CLEAN_MSA = sys.argv[4]
GAP_THRESHOLD = int(sys.argv[5])

SUBUNIT_CUTOFF = 685  # amino acid position between dsrA and dsrB


def main():
	records = SeqIO.parse(QUERY_MSA, format="fasta")

	gap_info = get_gap_info(records)

	records_dict = SeqIO.index(QUERY_MSA, format="fasta")

	clean_seqs(records_dict, gap_info)

	rename_queries(gap_info)

	dict_f = open(DUPS_DICT)
	dups_dict = json.load(dict_f)

	get_dups_info(dups_dict, gap_info)

	write_clean_msa(gap_info, CLEAN_MSA)

	write_csv(gap_info)


def get_gap_info(records):
	gaps_dict = {}

	for record in records:
		seq = str(record.seq)

		s_chars = seq_char(seq)
		s_start = s_chars[0]
		s_end = s_chars[1]
		s_len = s_chars[2]

		perc_gaps = calc_gaps(seq, s_start, s_end, s_len)

		gaps_dict[record.id] = {}
		gaps_dict[record.id]["p_gap"] = perc_gaps
		gaps_dict[record.id]["start"] = s_start
		gaps_dict[record.id]["end"] = s_end
		gaps_dict[record.id]["len"] = s_len

		region, trim_stat = id_region(seq, record.id, gaps_dict)

		gaps_dict[record.id]["region"] = region
		gaps_dict[record.id]["trim_stat"] = trim_stat

	return gaps_dict


def seq_char(seq):
	s_match = re.search(r'[a-z]', seq, re.I)
	s_start = s_match.start() + 1

	reverse_seq = seq[::-1]

	e_match = re.search(r'[a-z]', reverse_seq, re.I)
	s_end = (len(seq) - e_match.end()) + 1

	s_len = (s_end - s_start) + 1

	return s_start, s_end, s_len


def calc_gaps(seq, start, end, length):
	n_gaps = 0

	for residue in range(start - 1, end - 1):
		if seq[residue] == "-":
			n_gaps += 1

	return round((n_gaps / length) * 100)


def id_region(seq, query_id, gaps_dict):
	if gaps_dict[query_id]["end"] <= SUBUNIT_CUTOFF:
		gaps_dict[query_id]["a_gap"] = gaps_dict[query_id]["p_gap"]
		gaps_dict[query_id]["b_gap"] = "NA"
		return "A", "keep"

	elif gaps_dict[query_id]["start"] > SUBUNIT_CUTOFF:
		gaps_dict[query_id]["b_gap"] = gaps_dict[query_id]["p_gap"]
		gaps_dict[query_id]["a_gap"] = "NA"
		return "B", "keep"

	else:
		a_gap = calc_gaps(seq, gaps_dict[query_id]["start"], SUBUNIT_CUTOFF, (SUBUNIT_CUTOFF - gaps_dict[query_id]["start"]))
		b_gap = calc_gaps(seq, SUBUNIT_CUTOFF, gaps_dict[query_id]["end"], (gaps_dict[query_id]["end"] - SUBUNIT_CUTOFF))

		gaps_dict[query_id]["a_gap"] = a_gap
		gaps_dict[query_id]["b_gap"] = b_gap

		return check_gap_threshold(a_gap, b_gap)


def check_gap_threshold(a_gap, b_gap):
	if a_gap > GAP_THRESHOLD and b_gap > GAP_THRESHOLD:
		return "gap threshold not met", "remove"

	elif b_gap > GAP_THRESHOLD and a_gap < GAP_THRESHOLD:
		return "A", "trim B"

	elif a_gap > GAP_THRESHOLD and b_gap < GAP_THRESHOLD:
		return "B", "trim A"

	elif a_gap < GAP_THRESHOLD and b_gap < GAP_THRESHOLD:
		return "AB", "keep"


def clean_seqs(records_dict, gap_info):
	for query in gap_info.keys():
		original_seq = str(records_dict[query].seq)
		if gap_info[query]["trim_stat"] == "keep":
			new_seq = original_seq

		elif gap_info[query]["trim_stat"] == "trim A":
			new_seq = trimmer("A", original_seq)

		elif gap_info[query]["trim_stat"] == "trim B":
			new_seq = trimmer("B", original_seq)

		elif gap_info[query]["trim_stat"] == "remove":
			new_seq = "NA"

		gap_info[query]["seq"] = new_seq


def trimmer(region, seq):
	if region == "B":
		return seq[0:SUBUNIT_CUTOFF] + ("-" * (len(seq) - SUBUNIT_CUTOFF))

	if region == "A":
		return ("-" * SUBUNIT_CUTOFF) + seq[SUBUNIT_CUTOFF:len(seq)]


def write_clean_msa(gap_info, fasta_file):
	with open(fasta_file, "w+") as output_file:
		for query in gap_info:
			if gap_info[query]["trim_stat"] != "remove":
				sequence = gap_info[query]["seq"]
				seq_id = gap_info[query]["sample_name"]
				output_file.write(">" + seq_id + "\n" + sequence + "\n")


def rename_queries(gap_info):
	n_sample = 1
	for query in gap_info.keys():
		gap_info[query]["sample_name"] = "S" + str(n_sample)

		n_sample += 1


def get_dups_info(dups_dict, query_dict):
	for query in query_dict.keys():
		for dup_query in dups_dict.keys():
			if query in dup_query:
				exact_dups = dup_query.split(",")
				shorter_dups = dups_dict[dup_query].split(";")

				n_exact_dups = len(exact_dups) - 1
				n_shorter_dups = len(shorter_dups)

				query_dict[query]["n_exactdups"] = n_exact_dups
				query_dict[query]["n_shorterdups"] = n_shorter_dups
				query_dict[query]["exactdups"] = exact_dups
				query_dict[query]["shorterdups"] = shorter_dups


def write_csv(gap_dict):
	id_l = []
	pgap_l = []
	start_l = []
	end_l = []
	len_l = []
	agap_l = []
	bgap_l = []
	reg_l = []
	trim_l = []
	seq_l = []
	sample_l = []
	n_edups_l = []
	n_sdups_l = []
	edups_l = []
	sdups_l = []

	# with open(QUERY_INFO, "w+") as csvfile:
	csvfile = open(QUERY_INFO, "w+")
	writer = csv.writer(csvfile)

	for key, value in gap_dict.items():
		id_l.append(key)
		pgap_l.append(gap_dict[key]["p_gap"])
		start_l.append(gap_dict[key]["start"])
		end_l.append(gap_dict[key]["end"])
		len_l.append(gap_dict[key]["len"])
		agap_l.append(gap_dict[key]["a_gap"])
		bgap_l.append(gap_dict[key]["b_gap"])
		reg_l.append(gap_dict[key]["region"])
		trim_l.append(gap_dict[key]["trim_stat"])
		seq_l.append(gap_dict[key]["seq"])
		sample_l.append(gap_dict[key]["sample_name"])
		n_edups_l.append(gap_dict[key]["n_exactdups"])
		n_sdups_l.append(gap_dict[key]["n_shorterdups"])
		edups_l.append(gap_dict[key]["exactdups"])
		sdups_l.append(gap_dict[key]["shorterdups"])

	writer.writerow(["id", "pgap", "start", "end", "len", "agap", "bgap", "region", "trim_stat", "seq",
					"sample_name", "n_edups", "n_sdups", "edups", "sdups"])
	for i in range(len(id_l)):
		writer.writerow([id_l[i],
						pgap_l[i],
						start_l[i],
						end_l[i],
						len_l[i],
						agap_l[i],
						bgap_l[i],
						reg_l[i],
						trim_l[i],
						seq_l[i],
						sample_l[i],
						n_edups_l[i],
						n_sdups_l[i],
						edups_l[i],
						sdups_l[i]])


if __name__ == '__main__':
	main()