# root output tree from fragment insertion without Karthik seqs
# this will be the refTree for estimating model params
rule rootTree:
	input:
		join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap.newick")
	output:
		join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_rooted.newick")
	shell:
		"""
		python3 workflow/scripts/rootTree.py {input} 
		"""

# GET MODEL PARAMS
rule estimateModelParams:
	input:
		refMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
		refTree=config["hitAndRefResultTree"]
	output:
		modelParams=join(config["modelParamsDir"], f"RAxML_binaryModelParameters.{config['treeFileExtensionNovel']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtensionNovel']}")
	conda:
		"../envs/raxml.yml"
	threads: 8
	resources:
		time=75,
		mem_mb=16000
	params:
		outputDir=config["modelParamsDir"]
		treeFileExtension=config["treeFileExtensionNovel"]
	shell:
		"""
		raxmlHPC -f e -m PROTGAMMADAYHOFF -T 8 -s {input.refMSA} -t {input.refTree} -n {params.treeFileExtension} -w $(pwd)/workflow/out/modelParams
		"""

# FRAGMENT INSERTION
# just use capsules MSA without 2018 sequences since it is faster and we just want the original labelled tree at the end of this
rule runRAXML:
	input:
		seqAlignment=join(config["cleanHitsDir"],"StoolCapsule_compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa"),
		modelParams=join(config["modelParamsDir"], f"R AxML_binaryModelParameters.{config['treeFileExtensionNovel']}"),
		refTree=join(config["modelParamsDir"], f"RAxML_result.{config['treeFileExtensionNovel']}")
	params:
		outputDir=config["raxmlOutputDir"],
		fileExtension=config["treeFileExtension_frag_capsule"]
	output:
		join(config["raxmlOutputDir"], f"RAxML_info.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classification.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_classificationLikelihoodWeights.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_entropy.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_labelledTree.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_originalLabelledTree.raxmlnovel.{config['treeFileExtension_frag_capsule']}"),
		join(config["raxmlOutputDir"], f"RAxML_portableTree.raxmlnovel.{config['treeFileExtension_frag_capsule']}.jplace")
	conda:
		"../envs/raxml.yml"
	threads: 4
	resources:
		time=150,
		mem_mb=8000
	shell:
		"""
		raxmlHPC -f v -T {threads} -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w $(pwd)/workflow/out/raxmlOutput
		"""

#rule removeBootstrapValues:
#	input:
#		join(config["raxmlOutputDir"],"RAxML_labelledTree.raxmlnovel")
#	output:
#		join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_wKarthikSeqs.newick")
#	shell:
#		"""
#		python3 workflow/scripts/removeBootstrapValues.py {input} {output}
#		"""

#rule rootTree:
#	input:
#		join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_wKarthikSeqs.newick")
#	output:
#		join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_rooted_wKarthikSeqs.newick")
#	shell:
#		"""
#		python3 workflow/scripts/rootTree.py {input} {output}
#		"""
