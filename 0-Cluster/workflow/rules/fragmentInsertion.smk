# estimating model parameters
rule estimateModelParams:
	input:
		refMSA=config["refMSA"],
		refTree="config/dsrAB_consensus_phylogeny.newick"
	output:
		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtension']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtension']}")
	conda:
		"../envs/raxml.yml"
	threads: 8
	resources:
		time=30,
		mem_mb=16000
	params:
		outputDir=config["modelParamsDir"],
		treeFileExtension=config["treeFileExtension"]
	shell:
		"""
		raxmlHPC -f e -m PROTGAMMADAYHOFF -T 8 -s {input.refMSA} -t {input.refTree} -n {params.treeFileExtension} -w $(pwd)/workflow/out/modelParams
		"""

# FRAGMENT INSERTION
rule runRAXML_NCBI_no_2018_assemblies:
	input:
		seqAlignment=join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_assemblies_noDups_msa_withRef_trimmedGaps.faa"),
		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtension']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtension']}")
	params:
		outputDir=config["raxmlOutputDir"],
		fileExtension=config["treeFileExtension_no_2018_assemblies"]
	conda:
		"../envs/raxml.yml"
	threads: 4
	resources:
		time=700,
		mem_mb=10000
	output:
		join(config["raxmlOutputDir"], f"RAxML_info.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_classification.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_classificationLikelihoodWeights.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_entropy.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.{config['treeFileExtension_no_2018_assemblies']}"),
		join(config["raxmlOutputDir"], f"RAxML_portableTree.{config['treeFileExtension_no_2018_assemblies']}.jplace")
	shell:
		"""
		raxmlHPC -f v -T {threads} -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w $(pwd)/workflow/out/raxmlOutput
		"""

# CLEAN TREE
rule removeBootstrapValues:
	input:
		original=join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.{config['treeFileExtension_no_2018_assemblies']}"),
		labelled=join(config["raxmlOutputDir"], f"RAxML_labelledTree.{config['treeFileExtension_no_2018_assemblies']}")
	output:
		original=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_no_2018_assemblies.newick"),
		labelled=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_no_2018_assemblies.newick")
	resources:
		time=10,
		mem_mb=500
	shell:
		"""
		python3 workflow/scripts/removeBootstrapValues.py {input.original} {output.original}
		python3 workflow/scripts/removeBootstrapValues.py {input.labelled} {output.labelled}
		"""
