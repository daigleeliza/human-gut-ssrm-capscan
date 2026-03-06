rule runHMMER:
	input:
		join(config["prodigalDir"],"coassemblies/{coassembly}_contigs_out.faa")
	output:
		hmmOut_asr_A = join(config["asrHMMERDir"],"out/{coassembly}_asr_A.hmm.out"),
		domOut_asr_A = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_A.domtblout"),
		hmmOut_asr_B = join(config["asrHMMERDir"],"out/{coassembly}_asr_B.hmm.out"),
		domOut_asr_B = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_B.domtblout"),
		hmmOut_asr_C = join(config["asrHMMERDir"],"out/{coassembly}_asr_C.hmm.out"),
		domOut_asr_C = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_C.domtblout"),
		msa_asr_A = join(config["asrHMMERDir"],"sto/{coassembly}_asr_A.sto"),
		msa_asr_B = join(config["asrHMMERDir"],"sto/{coassembly}_asr_B.sto"),
		msa_asr_C = join(config["asrHMMERDir"],"sto/{coassembly}_asr_C.sto")
	threads: 2
	resources:
		mem_mb = 6000,
		time = 120
	params:
		hmm_profile_asr_A="config/asrA.HMM",
		hmm_profile_asr_B="config/asrB.HMM",
		hmm_profile_asr_C="config/asrC.HMM"
	conda:
		"../../workflow/envs/hmmer.yml"
	shell:
		"""
		hmmsearch -o {output.hmmOut_asr_A} --domtblout {output.domOut_asr_A} -A {output.msa_asr_A} {params.hmm_profile_asr_A} {input}
		hmmsearch -o {output.hmmOut_asr_B} --domtblout {output.domOut_asr_B} -A {output.msa_asr_B} {params.hmm_profile_asr_B} {input}
		hmmsearch -o {output.hmmOut_asr_C} --domtblout {output.domOut_asr_C} -A {output.msa_asr_C} {params.hmm_profile_asr_C} {input}
		"""

rule parseHMMER:
	input:
		domOut_asr_A = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_A.domtblout"),
		domOut_asr_B = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_B.domtblout"),
		domOut_asr_C = join(config["asrHMMERDir"],"domtbl/{coassembly}_asr_C.domtblout")
	output:
		asr_A=join(config["asrHMMERDir"], "summary/{coassembly}_asr_A_hits.csv"),
		asr_B=join(config["asrHMMERDir"], "summary/{coassembly}_asr_B_hits.csv"),
		asr_C=join(config["asrHMMERDir"], "summary/{coassembly}_asr_C_hits.csv")
	conda:
		"../../workflow/envs/biopython.yml"
	shell:
		"""
		python3 workflow/scripts/parse_hmmer_domtable.py {input.domOut_asr_A} {output.asr_A}
		python3 workflow/scripts/parse_hmmer_domtable.py {input.domOut_asr_B} {output.asr_B}
		python3 workflow/scripts/parse_hmmer_domtable.py {input.domOut_asr_C} {output.asr_C}
		"""

rule combineCSV_ASR:
	input:
		asr_A=expand(join(config["asrHMMERDir"],"summary/{coassembly}_asr_A_hits.csv"), coassembly=coassemblies),
		asr_B=expand(join(config["asrHMMERDir"],"summary/{coassembly}_asr_B_hits.csv"), coassembly=coassemblies),
		asr_C=expand(join(config["asrHMMERDir"],"summary/{coassembly}_asr_C_hits.csv"), coassembly=coassemblies)
	output:
		asr_A=join(config["asrHMMERDir"],"compiled_asr_A_hits.csv"),
		asr_B=join(config["asrHMMERDir"],"compiled_asr_B_hits.csv"),
		asr_C=join(config["asrHMMERDir"],"compiled_asr_C_hits.csv")
	shell:
		"""
		cat {input.asr_A} > {output.asr_A}
		cat {input.asr_B} > {output.asr_B}
		cat {input.asr_C} > {output.asr_C}
		"""

rule filterHMM_ASR:
	input:
		asr_A=join(config["asrHMMERDir"],"compiled_asr_A_hits.csv"),
		asr_B=join(config["asrHMMERDir"],"compiled_asr_B_hits.csv"),
		asr_C=join(config["asrHMMERDir"],"compiled_asr_C_hits.csv")
	output:
		asr_A=join(config["asrHMMERDir"],"compiled_asr_A_hits_bitFilter.csv"),
		asr_B=join(config["asrHMMERDir"],"compiled_asr_B_hits_bitFilter.csv"),
		asr_C=join(config["asrHMMERDir"],"compiled_asr_C_hits_bitFilter.csv")
	threads: 1
	resources:
		mem_mb=1000,
		time=30
	params:
		asr_A_cutoff=config["asr_A_tc"],
		asr_B_cutoff=config["asr_B_tc"],
		asr_C_cutoff=config["asr_C_tc"]
	shell:
		"""
		awk -F',' '$8 >= {params.asr_A_cutoff}' {input.asr_A} > {output.asr_A}
		awk -F',' '$8 >= {params.asr_B_cutoff}' {input.asr_B} > {output.asr_B}
		awk -F',' '$8 >= {params.asr_C_cutoff}' {input.asr_C} > {output.asr_C}
		"""


# awk statement to remove info on hmmer search and to include the coassembly at the beginning of the line
rule testConvert:
	input: 
		A=join(config["asrHMMERDir"],"sto/{coassembly}_asr_A.sto"),
		B=join(config["asrHMMERDir"],"sto/{coassembly}_asr_B.sto"),
		C=join(config["asrHMMERDir"],"sto/{coassembly}_asr_C.sto")
	output: 
		A=join(config["asrHMMERDir"],"faa/{coassembly}_asr_A.faa"),
		B=join(config["asrHMMERDir"],"faa/{coassembly}_asr_B.faa"),
		C=join(config["asrHMMERDir"],"faa/{coassembly}_asr_C.faa")
	conda: "../../workflow/envs/hmmer.yml"
	threads: 1
	resources:
		mem_mb=1000,
		runtime=30
	shell:
		"""
		esl-reformat -u fasta {input.A} |
		awk -v coassembly="{wildcards.coassembly}" '
			/^>/ {{
				sub(/\\s*\\[.*$/, "", $0)
				contig = substr($0, 2)
				print ">" coassembly "." contig
				next
			}}
			{{print}}
		' > {output.A}
		esl-reformat -u fasta {input.B} |
		awk -v coassembly="{wildcards.coassembly}" '
			/^>/ {{
				sub(/\\s*\\[.*$/, "", $0)
				contig = substr($0, 2)
				print ">" coassembly "." contig
				next
			}}
			{{print}}
		' > {output.B}
		esl-reformat -u fasta {input.C} |
		awk -v coassembly="{wildcards.coassembly}" '
			/^>/ {{
				sub(/\\s*\\[.*$/, "", $0)
				contig = substr($0, 2)
				print ">" coassembly "." contig
				next
			}}
			{{print}}
		' > {output.C}
		"""

rule combineFAA_ASR:
	input:
		asr_A=expand(join(config["asrHMMERDir"],"faa/{coassembly}_asr_A.faa"), coassembly=coassemblies),
		asr_B=expand(join(config["asrHMMERDir"],"faa/{coassembly}_asr_B.faa"), coassembly=coassemblies),
		asr_C=expand(join(config["asrHMMERDir"],"faa/{coassembly}_asr_C.faa"), coassembly=coassemblies)
	output:
		asr_A=join(config["asrHMMERDir"],"compiled_asr_A_hits.faa"),
		asr_B=join(config["asrHMMERDir"],"compiled_asr_B_hits.faa"),
		asr_C=join(config["asrHMMERDir"],"compiled_asr_C_hits.faa")
	conda: "../envs/biopython.yml"
	threads: 1
	resources:
		mem_mb=1000,
		runtime=30
	shell:
		"""
		cat {input.asr_A} > {output.asr_A}
		cat {input.asr_B} > {output.asr_B}
		cat {input.asr_C} > {output.asr_C}
		"""
