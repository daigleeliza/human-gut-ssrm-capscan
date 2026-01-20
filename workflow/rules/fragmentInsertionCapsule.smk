# make MSA from our result msa (with hits and ref seqs) with Anantharaman2018's novel seqs
rule runMAFFT_2:
	input:
		hitsWithRefMSA=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
		Anantharaman2018Seqs="config/Anantharaman2018_dsrA_dsrB_noDups.faa"
	output:
		resultMSA=join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa")
	conda:
		"../envs/mafft.yml"
	threads: 4
	resources:
		time=30,
		mem_mb=2000
	shell:
		"""
		mafft --add {input.Anantharaman2018Seqs} --keeplength {input.hitsWithRefMSA} > {output.resultMSA}
		"""

#can reuse this from the NCBI rule
#rule estimateModelParams:
#	input:
#		refMSA=config["refMSA"],
#		refTree="config/dsrAB_consensus_phylogeny.newick"
#	output:
#		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtension']}"),
#		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtension']}")
#	conda:
#		"../envs/raxml.yml"
#	threads: 8
#	resources:
#		time=30,
#		mem_mb=16000
#	params:
#		outputDir=config["modelParamsDir"],
#		treeFileExtension=config["treeFileExtension"]
#	shell:
#		"""
#		raxmlHPC -f e -m PROTGAMMADAYHOFF -T 8 -s {input.refMSA} -t {input.refTree} -n {params.treeFileExtension} -w $(pwd)/workflow/out/modelParams
#		"""

# FRAGMENT INSERTION
rule runRAXML:
	input:
		seqAlignment=join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa"),
		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtension']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtension']}")
	params:
		outputDir=config["raxmlOutputDir"],
		#Daigle2026_frag_capsule
		fileExtension=config["treeFileExtension_frag_capsule"]
	conda:
		"../envs/raxml.yml"
	threads: 4
	resources:
		time=360,
		mem_mb=8000
	output:
		join(config["raxmlOutputDir"], f"RAxML_info.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classification.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classificationLikelihoodWeights.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_entropy.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_portableTree.{config['treeFileExtension_frag_capsule']}.jplace")
	shell:
		"""
		raxmlHPC -f v -T {threads} -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w $(pwd)/workflow/out/raxmlOutput
		"""

# CLEAN TREE
rule removeBootstrapValues:
	input:
		join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_frag_capsule']}"),
	output:
		join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_capsule.newick"),
	resources:
		time=5,
		mem_mb=500
	shell:
		"""
		python3 workflow/scripts/removeBootstrapValues.py {input} {output}
		"""
