## Python script to parse HMMER domtbl output into csv file for further analysis
## Based on code provided by RC Christensen on 13May2021

import sys
import csv
from Bio import SearchIO


output_domtbl = sys.argv[1]
COMPILED_TABLE_CSV = sys.argv[2]

domtbl_records = SearchIO.parse(output_domtbl, "hmmsearch3-domtab")
hit_id_list = []
targ_id_list = []
target_len_list = []
hit_score_list = []
hit_evalue_list = []
hit_bias_list = []
hsp_len_list = []
hsp_score_list = []
hsp_evalue_cond_list = []
hsp_evalue_ind_list = []
hsp_bias_list = []

f=open(COMPILED_TABLE_CSV, mode='a')
f.close()

for qresult in domtbl_records:
    for hit in qresult:
        if (hit.evalue <= 0.01):
            for hsp in hit:
                if (hsp.evalue_cond <= 0.01):
                    for hspfrag in hsp:

                        hit_id = hit.id + "/" + str(hspfrag.hit_start + 1) + "-" + str(hspfrag.hit_end)
                        hsp_len = (hsp.hit_end - hspfrag.hit_start)
                        hit_id_list.append(hit_id)
                        targ_id_list.append(hit.id)
                        target_len_list.append(hit.seq_len)
                        hit_score_list.append(hit.bitscore)
                        hit_evalue_list.append(hit.evalue)
                        hit_bias_list.append(hit.bias)
                        hsp_len_list.append(hsp_len)
                        hsp_evalue_cond_list.append(hsp.evalue_cond)
                        hsp_evalue_ind_list.append(hsp.evalue)
                        hsp_bias_list.append(hsp.bias)
                        hsp_score_list.append(hsp.bitscore)

    # add to csv file
    with open(COMPILED_TABLE_CSV, mode='a') as file:
        freq_writer = csv.writer(file, delimiter=",")
        for i in range(len(hit_id_list)):
            freq_writer.writerow([hit_id_list[i],
                                  targ_id_list[i],
                                  target_len_list[i],
                                  hit_score_list[i],
                                  hit_evalue_list[i],
                                  hit_bias_list[i],
                                  hsp_len_list[i],
                                  hsp_score_list[i],
                                  hsp_evalue_cond_list[i],
                                  hsp_evalue_ind_list[i],
                                  hsp_bias_list[i],
                                  output_domtbl])

