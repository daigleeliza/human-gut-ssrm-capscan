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

rule runRAXML_capsule_no_2018:
	input:
		seqAlignment=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtension']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtension']}")
	params:
		outputDir=config["raxmlOutputDir"],
		fileExtension=config["treeFileExtension_no_2018_capsule"]
	conda:
		"../envs/raxml.yml"
	threads: 4
	resources:
		time=360,
		mem_mb=8000
	output:
		join(config["raxmlOutputDir"], f"RAxML_info.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classification.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classificationLikelihoodWeights.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_entropy.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.{config['treeFileExtension_no_2018_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_portableTree.{config['treeFileExtension_no_2018_capsule']}.jplace")
	shell:
		"""
		raxmlHPC -f v -T {threads} -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w $(pwd)/workflow/out/raxmlOutput
		"""

# CLEAN TREE
rule removeBootstrapValues:
	input:
		original=join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.{config['treeFileExtension_no_2018_capsule']}"),
		labelled=join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_no_2018_capsule']}")
	output:
		original=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_no_2018_capsule.newick"),
		labelled=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_no_2018_capsule.newick")
	resources:
		time=10,
		mem_mb=500
	shell:
		"""
		python3 workflow/scripts/removeBootstrapValues.py {input.original} {output.original}
		python3 workflow/scripts/removeBootstrapValues.py {input.labelled} {output.labelled}
		"""
