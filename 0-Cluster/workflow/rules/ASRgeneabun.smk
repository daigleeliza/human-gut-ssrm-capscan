# fix formatting of csv
rule fixCSV:
	input: 
		hmm=join(config["asrHMMERDir"],"compiled_asr_{ABC}_hits_bitFilter.csv")
	output: 
		hmm=join(config["asrHMMERDir"],"compiled_asr_{ABC}_hits_bitFilter_formatted.csv")
	threads: 1
	resources:
		mem_mb=4000,
		runtime=30
	shell:
		"""
		awk -F',' 'BEGIN {{ OFS="," }}
		{{
			if (match($NF, /(Stool_[0-9]+|Capsule_[0-9]+)/, a)) {{
				$1 = a[1] "." $1
				$2 = a[1] "." $2
			}}
			print
		}}
		' {input.hmm} > {output.hmm}
		""" 

#extract geneReads that are found in filtered csv
rule countASR:
	input:
		hmm=expand(join(config["asrHMMERDir"],"compiled_asr_{ABC}_hits_bitFilter_formatted.csv"), ABC=["A","B","C"]),
		geneReads=join(config["geneReadsDir"],"{coassembly}/{sample}_microbeCensus_RPKG.txt")
	output: "workflow/out/asrABC_CapsuleStool_Abundances/{coassembly}/{sample}_asr_ABC_RPKG_redo.txt"
	threads: 1
	resources:
		mem_mb=4000,
		runtime=60
	shell:
		"""
		set -euo pipefail
		temp_file=$(mktemp)
		awk 'NR==2{{print $0"\\tgeneLoc\\tgeneNum\\tScore"}}' {input.geneReads} > $temp_file
		for hmm_file in {input.hmm}; do
			awk -F',' '{{print $1","$2","$8}}' "$hmm_file" | while IFS=',' read -r col1 col2 col8; do
				sample_subj=$(echo "$col2" | awk -F'.' '{{print $1}}')
		
				if [ "$sample_subj" == "{wildcards.coassembly}" ]; then
					contig=$(echo "$col2" | awk -F'_' '{{print $2"_"$3}}' | sed 's/.*\\.//')
					gene_number=$(echo "$col2" | awk -F'_' '{{print $NF}}')
					tail -n +3 {input.geneReads} | awk -F'\\t' -v Chr="$contig" -v gene_number="$gene_number" -v c1="$col1" -v c2="$col2" -v c8="$col8" '
						$2 == Chr && substr($1, index($1,"_")+1) == gene_number {{print $0"\\t"c1"\\t"c2"\\t"c8}}
					' >> $temp_file
				fi
			done
		done
		mv $temp_file {output}
		"""

#concatenate all coassemblies' hits
rule concatGeneHits_asr:
	input:
		expand(
			join("workflow/out/asrABC_CapsuleStool_Abundances/{coassembly}/{sample}_asr_ABC_RPKG_redo.txt"),
			zip,
			coassembly=[c for c in coassemblies for _ in get_subject_sample_list_dropped(c)],
			sample=[s for c in coassemblies for s in get_subject_sample_list_dropped(c)])
	output:
		"workflow/out/asrABC_CapsuleStool_Abundances/concat_asrABC_Abundances_RPKG_redo.txt"
	threads: 1
	resources:
		mem_mb=2000,
		time=120
	#make one file including all samples: keep only tab separated portion of input, add column for identifying sample, rename columns
	shell:
		"""
		echo -e "coassembly\\tsample\\tgeneID\\tcontig\\tstart\\tstop\\tstrand\\tlength\\tmapped\\tRPKG\\tgeneLoc\\tgeneNum\\tScore" > {output}
		for x in {input}; do
			coassembly=$(basename "$(dirname "$x")")
			sample=$(basename "$x" | awk -F'_' '{{print $1"_"$2}}')
			tail -n +2 "$x" | awk -F'\\t' -v c="$coassembly" -v s="$sample" 'BEGIN{{OFS="\\t"}} {{
				split($2, arr2, "_")
				split($1, arr1, "_")
				geneID = arr2[2] "_" arr1[2]
				print c,s,geneID,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
			}}'
		done >> {output}
		"""
