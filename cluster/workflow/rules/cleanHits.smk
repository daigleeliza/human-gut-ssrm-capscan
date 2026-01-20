rule scoreFilter_NCBI:
	input:
		bactCSV = join(config["summaryDir"],"compiled_dsrAB_bact_hits.csv"),
		archCSV = join(config["summaryDir"],"compiled_dsrAB_arch_hits.csv"),
		bactFASTA=join(config["summaryDir"],"compiled_dsrAB_bact_hits.faa"),
		archFASTA=join(config["summaryDir"],"compiled_dsrAB_arch_hits.faa")
	output:
		bactCSV = join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.csv"),
		archCSV = join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.csv"),
		bactFASTA = join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.faa"),
		archFASTA = join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.faa")
	conda:
		"../envs/biopython.yml"
	resources:
		mem_mb=4000,
		time=60
	params:
		scripts_dir=config["scriptsDir"],
		score_threshold=config["scoreThreshold"]
	shell:
		"""
		python3 {params.scripts_dir}/scoreFilter.py {input.bactFASTA} {input.bactCSV} {output.bactFASTA} {output.bactCSV} {params.score_threshold}
		python3 {params.scripts_dir}/scoreFilter.py {input.archFASTA} {input.archCSV} {output.archFASTA} {output.archCSV} {params.score_threshold}
		"""

rule combineArchBact_NCBI:
	input:
		bactCSV = join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.csv"),
		archCSV = join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.csv"),
		bactFASTA = join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.faa"),
		archFASTA = join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.faa")
	output:
		combinedFASTA=join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold.faa"),
		combinedCSV=join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold.csv")
	conda:
		"../envs/biopython.yml"
	resources:
		time=30,
		mem_mb=2000
	shell:
		"""
		cat {input.bactFASTA} {input.archFASTA} > {output.combinedFASTA}
		cat {input.bactCSV} {input.archCSV} > {output.combinedCSV}
		"""

rule removeDuplicates_NCBI:
	input:
		join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold.faa")
	output:
		join(config["cleanHitsDir"], "compiled_dsrAB_hits_scoreThreshold_noDups.faa"),
		join(config["cleanHitsDir"], "compiled_dsrAB_hits_scoreThreshold_noDups.json"),
		join(config["cleanHitsDir"], "compiled_dsrAB_hits_scoreThreshold_noEDups.faa")
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir = config["scriptsDir"],
		output_dir = config["cleanHitsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/removeDuplicates.py {input} {params.output_dir}
		"""

rule runMAFFT_NCBI:
	input:
		hitsFASTA = join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_noDups.faa"),
		refMSA=config["refMSA"]
	output:
		withRefMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef.faa"),
		noRefMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_noRef.faa")
	conda:
		"../envs/mafft.yml"
	threads: 4
	resources:
		time=30,
		mem_mb=2000
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		mafft --add {input.hitsFASTA} --keeplength {input.refMSA} > {output.withRefMSA}
		python3 {params.scripts_dir}/removeRefFromMSA.py {output.withRefMSA} {input.refMSA} {output.noRefMSA}
		"""

rule trimGaps_identifySubunit_NCBI:
	input:
		MSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_noRef.faa"),
		dupsInfo=join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_noDups.json"),
	output:
		gapPercCSV=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
		trimmedGapsMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_noRef_trimmedGaps.faa")
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir = config["scriptsDir"],
		gapThreshold = config["gapThreshold"]
	shell:
		"""
		python3 {params.scripts_dir}/trimGaps_identifySubunit.py {input.MSA} {input.dupsInfo} {output.gapPercCSV} {output.trimmedGapsMSA} {params.gapThreshold}
		"""

rule compileFinalMSA_NCBI:
	input:
		hitsMSA=join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_noDups_msa_noRef_trimmedGaps.faa"),
		refMSA=config["refMSA"]
	output:
		join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa")
	shell:
		"""
		cat {input.hitsMSA} {input.refMSA} > {output}
		"""
