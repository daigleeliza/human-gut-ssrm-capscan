#set microorganism=microorganism
#this is for filtering by score in column 8 (max. individual score instead of overall score)
rule filterHMMCapsule:
	input:
		join(config["dsrHMMERDir"],"coassemblies/CapsuleData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits.csv")
	output:
		join(config["dsrHMMERDir"],"coassemblies/CapsuleData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits_FILTERED_col8.csv")
	threads: 1
	resources:
		mem_mb=1000,
		time=30
	shell:
		"""
		awk -F',' '$8 >= 100' {input} > {output}
		"""

rule filterHMMStool:
	input:
		join(config["dsrHMMERDir"],"coassemblies/StoolData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits.csv")
	output:
		join(config["dsrHMMERDir"],"coassemblies/StoolData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits_FILTERED_col8.csv")
	threads: 1
	resources:
		mem_mb=1000,
		time=30
	shell:
		"""
		awk -F',' '$8 >= 100' {input} > {output}
		"""

rule countHMMinSamplesCapsule:
	input:
		hmm=expand(join(config["dsrHMMERDir"],"coassemblies/CapsuleData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits_FILTERED_col8.csv"), microorganism=microorganism),
		geneReads=join(config["geneReadsDir"],"{coassembly}/{sample}_RPKM.txt")
	output: "workflow/out/dsrAB_CapsuleStool_Abundances/Capsule/{coassembly}/{sample}_dsrAB_col8score_RPKM.txt"
			#"workflow/out/dsrAB_CapsuleStool_Abundances/{coassembly}/{sample}_dsrAB_RPKM.txt"
	threads: 1
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		temp_file=$(mktemp)
		awk 'NR==2{{print $0"\\tgeneLoc\\tgeneNum\\tScore"}}' {input.geneReads} > $temp_file
		for hmm_file in {input.hmm}; do
			awk -F',' '{{print $1","$2","$8}}' "$hmm_file" | while IFS=',' read -r col1 col2 col8; do
				sample_subj=$(echo "$col2" | awk -F'.' '{{print $1}}')
		
				if [ "$sample_subj" == "{wildcards.coassembly}" ]; then
					contig=$(echo "$col2" | awk -F'_' '{{print $2"_"$3}}' | sed 's/.*\\.//')
					gene_number=$(echo "$col2" | awk -F'_' '{{print $NF}}')
					tail -n +2 {input.geneReads} | awk -F'\\t' -v Chr="$contig" -v gene_number="$gene_number" -v c1="$col1" -v c2="$col2" -v c8="$col8" '
						$2 == Chr && $1 ~ "_"gene_number {{print $0"\\t"c1"\\t"c2"\\t"c8}}
					' >> $temp_file
				fi
			done
		done
		mv $temp_file {output}
		"""

rule concatGeneHitsCapsule:
	input:
		expand(
			join("workflow/out/dsrAB_CapsuleStool_Abundances/Capsule/{coassembly}/{sample}_dsrAB_col8score_RPKM.txt"),
			zip,
			coassembly=[c for c in capsule for _ in get_subject_sample_list_dropped(c)],
			sample=[s for c in capsule for s in get_subject_sample_list_dropped(c)])
	output:
		"workflow/out/dsrAB_CapsuleStool_Abundances/Capsule/concat_dsrAB_col8score_Capsule_Abundances_RPKM.txt"
		#"workflow/out/dsrAB_CapsuleStool_Abundances/concat_dsrAB_Capsule_Abundances_RPKM.txt"
	threads: 1
	resources:
		mem_mb=2000,
		time=30
	#make one file including all capsule samples: keep only tab separated portion of input, add column for identifying sample, rename columns
	shell:
		"""
		echo -e "coassembly\\tsample\\tgeneID\\tcontig\\tstart\\tstop\\tstrand\\tlength\\tmapped\\tRPKM\\tgeneLoc\\tgeneNum\\tScore" > {output}
		for x in {input}; do
			coassembly=$(basename "$(dirname "$x")")
			sample=$(basename "$x" _dsrAB_RPKM.txt)
			tail -n +3 "$x" | awk -F'\\t' -v c="$coassembly" -v s="$sample" 'BEGIN{{OFS="\\t"}} {{print c,s,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}'
		done >> {output}
		"""

		#for column 4 filtering echo -e "coassembly\\tsample\\tgeneID\\tcontig\\tstart\\tstop\\tstrand\\tlength\\tmapped\\tRPKM\\tgeneLoc\\tgeneNum\\tbitScore" > {output}

rule countHMMinSamplesStool:
	input:
		hmm=expand(join(config["dsrHMMERDir"],"coassemblies/StoolData_parsedHMMDomainTableSummary_dsrAB_{microorganism}_hits_FILTERED_col8.csv"), microorganism=microorganism),
		geneReads=join(config["geneReadsDir"],"{coassembly}/{sample}_RPKM.txt")
	output: "workflow/out/dsrAB_CapsuleStool_Abundances/Stool/{coassembly}/{sample}_dsrAB_col8score_RPKM.txt"
			#"workflow/out/dsrAB_CapsuleStool_Abundances/{coassembly}/{sample}_dsrAB_RPKM.txt"
	threads: 1
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		temp_file=$(mktemp)
		awk 'NR==2{{print $0"\\tgeneLoc\\tgeneNum\\tScore"}}' {input.geneReads} > $temp_file
		for hmm_file in {input.hmm}; do
			awk -F',' '{{print $1","$2","$8}}' "$hmm_file" | while IFS=',' read -r col1 col2 col8; do
				sample_subj=$(echo "$col2" | awk -F'.' '{{print $1}}')
		
				if [ "$sample_subj" == "{wildcards.coassembly}" ]; then
					contig=$(echo "$col2" | awk -F'_' '{{print $2"_"$3}}' | sed 's/.*\\.//')
					gene_number=$(echo "$col2" | awk -F'_' '{{print $NF}}')
					tail -n +2 {input.geneReads} | awk -F'\\t' -v Chr="$contig" -v gene_number="$gene_number" -v c1="$col1" -v c2="$col2" -v c8="$col8" '
						$2 == Chr && $1 ~ "_"gene_number {{print $0"\\t"c1"\\t"c2"\\t"c8}}
					' >> $temp_file
				fi
			done
		done
		mv $temp_file {output}
		"""

rule concatGeneHitsStool:
	input:
		expand(
			join("workflow/out/dsrAB_CapsuleStool_Abundances/Stool/{coassembly}/{sample}_dsrAB_col8score_RPKM.txt"),
			zip,
			coassembly=[c for c in stool for _ in get_subject_sample_list_dropped(c)],
			sample=[s for c in stool for s in get_subject_sample_list_dropped(c)])
	output:
		"workflow/out/dsrAB_CapsuleStool_Abundances/Stool/concat_dsrAB_col8score_Stool_Abundances_RPKM.txt"
		#"workflow/out/dsrAB_CapsuleStool_Abundances/concat_dsrAB_Stool_Abundances_RPKM.txt"
	threads: 1
	resources:
		mem_mb=2000,
		time=30
	#make one file including all stool samples: keep only tab separated portion of input, add column for identifying sample, rename columns
	shell:
		"""
		echo -e "coassembly\\tsample\\tgeneID\\tcontig\\tstart\\tstop\\tstrand\\tlength\\tmapped\\tRPKM\\tgeneLoc\\tgeneNum\\tScore" > {output}
		for x in {input}; do
			coassembly=$(basename "$(dirname "$x")")
			sample=$(basename "$x" _dsrAB_RPKM.txt)
			tail -n +3 "$x" | awk -F'\\t' -v c="$coassembly" -v s="$sample" 'BEGIN{{OFS="\\t"}} {{print c,s,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}}'
		done >> {output}
		"""
