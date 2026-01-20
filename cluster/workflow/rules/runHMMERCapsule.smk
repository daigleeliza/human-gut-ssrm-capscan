# set coassembly = coassemblies
rule runProdigal:
	input: join(config["coAssemblyDir"],"{coassembly}/megahit/filtered_contigs_1kb.fasta")
	output:
		geneCoords=join(config["prodigalGeneCoordDir"],"{coassembly}_geneCoord.out"),
		proteinSeqs=join(config["prodigalProteinSeqDir"], "{coassembly}.faa")
	conda:
		"../envs/prodigal.yml"
	shell:
		"""
		prodigal -i {input} -o {output.geneCoords} -a {output.proteinSeqs}
		"""

rule runHMMER:
	input:
		join(config["prodigalProteinSeqDir"],"{coassembly}.faa")
	output:
		hmmOut_bact = join(config["hmmOutDir"],"{coassembly}_dsrAB_bact.hmm.out"),
		domOut_bact = join(config["domtblDir"],"{coassembly}_dsrAB_bact.domtblout"),
		hmmOut_arch = join(config["hmmOutDir"],"{coassembly}_dsrAB_arch.hmm.out"),
		domOut_arch = join(config["domtblDir"],"{coassembly}_dsrAB_arch.domtblout"),
		msa_bact= join(config["msaDir"],"sto/{coassembly}_dsrAB_bact.sto"),
		msa_arch= join(config["msaDir"],"sto/{coassembly}_dsrAB_arch.sto")
	params:
		hmm_profile_bact=config["dsrAB_bact_HMMProfile"],
		hmm_profile_arch=config["dsrAB_arch_HMMProfile"]
	conda:
		"../envs/hmmer.yml"
	shell:
		"""
		hmmsearch -o {output.hmmOut_bact} --domtblout {output.domOut_bact} -A {output.msa_bact} {params.hmm_profile_bact} {input}
		hmmsearch -o {output.hmmOut_arch} --domtblout {output.domOut_arch} -A {output.msa_arch} {params.hmm_profile_arch} {input}
		"""

rule parseHMMER:
	input:
		bact=join(config["domtblDir"],"{coassembly}_dsrAB_bact.domtblout"),
		arch=join(config["domtblDir"],"{coassembly}_dsrAB_arch.domtblout")
	output:
		bact=join(config["summaryDir"], "{coassembly}_dsrAB_bact_hits.csv"),
		arch=join(config["summaryDir"], "{coassembly}_dsrAB_arch_hits.csv")
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/parse_hmmer_domtable.py {input.bact} {output.bact}
		python3 {params.scripts_dir}/parse_hmmer_domtable.py {input.arch} {output.arch}
		"""

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA:
	input:
		bact=join(config["msaDir"],"sto/{coassembly}_dsrAB_bact.sto"),
		arch=join(config["msaDir"],"sto/{coassembly}_dsrAB_arch.sto")
	output:
		bact=join(config["msaDir"],"faa/{coassembly}_dsrAB_bact.faa"),
		arch=join(config["msaDir"],"faa/{coassembly}_dsrAB_arch.faa")
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/convertFASTA.py {input.bact} {output.bact}
		python3 {params.scripts_dir}/convertFASTA.py {input.arch} {output.arch}
		"""

# files included in the repo starting below
# combine stool and capsule faa files
rule combineMSA:
	input: 
		bact=expand(join(config["dsrHMMERDir"], "faa/{coassembly}_dsrAB_bact.faa"), coassembly=coassemblies),
		arch=expand(join(config["dsrHMMERDir"], "faa/{coassembly}_dsrAB_arch.faa"), coassembly=coassemblies)
	output:
		bact=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_compiled_dsrAB_bact_hits.faa"),
		arch=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_compiled_dsrAB_arch_hits.faa")
	shell:
		"""
		cat {input.bact} > {output.bact}
		cat {input.arch} > {output.arch}
		"""
# combine stool and capsule csv files
rule combineCSV:
	input: 
		bact=expand(join(config["dsrHMMERDir"], "coassemblies/{sampletype}Data_parsedHMMDomainTableSummary_dsrAB_bacterial_hits.csv"), sampletype=sampletype),
		arch=expand(join(config["dsrHMMERDir"], "coassemblies/{sampletype}Data_parsedHMMDomainTableSummary_dsrAB_archaeal_hits.csv"), sampletype=sampletype)
	output:
		bact=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_parsedHMMDomainTableSummary_dsrAB_bacterial_hits.csv"),
		arch=join(config["dsrHMMERDir"],"summary/StoolCapsuleData_parsedHMMDomainTableSummary_dsrAB_archaeal_hits.csv")
	shell:
		"""
		cat {input.bact} > {output.bact}
		cat {input.arch} > {output.arch}
		"""
