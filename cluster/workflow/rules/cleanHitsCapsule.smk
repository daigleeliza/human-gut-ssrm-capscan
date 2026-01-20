# bit filter csv and faa files
rule scoreFilter:
	input:
		bactCSV = join(config["dsrHMMERDir"],"summary/StoolCapsuleData_parsedHMMDomainTableSummary_dsrAB_bacterial_hits.csv"),
		archCSV = join(config["dsrHMMERDir"],"summary/StoolCapsuleData_parsedHMMDomainTableSummary_dsrAB_archaeal_hits.csv"),
		bactFASTA=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_compiled_dsrAB_bact_hits.faa"),
		archFASTA=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_compiled_dsrAB_arch_hits.faa")
	output:
		bactCSV = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_bact_hits_scoreThreshold.csv"),
		archCSV = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_arch_hits_scoreThreshold.csv"),
		bactFASTA = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_bact_hits_scoreThreshold.faa"),
		archFASTA = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_arch_hits_scoreThreshold.faa")
	conda:
		"../envs/biopython.yml"
	log:
		bact="Eliza/logs/bitfilter/bact.log",
		arch="Eliza/logs/bitfilter/arch.log"
	params:
		score_threshold=100
	shell:
		"""
		python3 workflow/scripts/scoreFilter.py {input.bactFASTA} {input.bactCSV} {output.bactFASTA} {output.bactCSV} {params.score_threshold} >> {log.bact} 2>&1
		python3 workflow/scripts/scoreFilter.py {input.archFASTA} {input.archCSV} {output.archFASTA} {output.archCSV} {params.score_threshold} >> {log.arch} 2>&1
		"""

rule combineArchBact:
	input:
		bactCSV = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_bact_hits_scoreThreshold.csv"),
		archCSV = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_arch_hits_scoreThreshold.csv"),
		bactFASTA = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_bact_hits_scoreThreshold.faa"),
		archFASTA = join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_arch_hits_scoreThreshold.faa")
	output:
		combinedFASTA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold.faa"),
		combinedCSV= join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold.csv")
	shell:
		"""
		cat {input.bactFASTA} {input.archFASTA} > {output.combinedFASTA}
		cat {input.bactCSV} {input.archCSV} > {output.combinedCSV}
		"""

rule removeDuplicates:
	input:
		join(config["cleanHitsDir"], "StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold.faa")
	output:
		join(config["cleanHitsDir"], "StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold_noDups.faa"),
		join(config["cleanHitsDir"], "StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold_noDups.json"),
		join(config["cleanHitsDir"], "StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold_noEDups.faa")
	conda:
		"../envs/biopython.yml"
	params:
		output_dir = config["cleanHitsDir"]
	log: "Eliza/logs/nodups.log"
	shell:
		"""
		python3 workflow/scripts/removeDuplicates.py {input} {params.output_dir} >> {log} 2>&1
		"""

rule runMAFFT:
	input:
		hitsFASTA = join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold_noDups.faa"),
		refMSA=config["refMSA"]
	output:
		withRefMSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_withRef.faa"),
		noRefMSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_noRef.faa")
	conda:
		"../envs/mafft.yml"
	threads: 4
	resources:
		time=30,
		mem_mb=2000
	shell:
		"""
		mafft --add {input.hitsFASTA} --keeplength {input.refMSA} > {output.withRefMSA}
		python3 workflow/scripts/removeRefFromMSA.py {output.withRefMSA} {input.refMSA} {output.noRefMSA}
		"""

rule trimGaps_identifySubunit:
	input:
		MSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_noRef.faa"),
		dupsInfo=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold_noDups.json")
	output:
		gapPercCSV=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
		trimmedGapsMSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_noRef_trimmedGaps.faa")
	conda:
		"../envs/biopython.yml"
	params:
		gapThreshold = config["gapThreshold"]
	shell:
		"""
		python3 workflow/scripts/trimGaps_identifySubunit.py {input.MSA} {input.dupsInfo} {output.gapPercCSV} {output.trimmedGapsMSA} {params.gapThreshold}
		"""

rule compileFinalMSA:
	input:
		hitsMSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_noRef_trimmedGaps.faa"),
		refMSA=config["refMSA"]
	output:
		join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa")
	shell:
		"""
		cat {input.hitsMSA} {input.refMSA} > {output}
		"""
