# GET MODEL PARAMS
rule estimateModelParams:
	input:
		refMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
		refTree=config["hitAndRefResultTree"]
	params:
		outputDir=config["modelParamsDir"]
	shell:
		"""
		raxmlHPC -f e -m PROTGAMMADAYHOFF -s {input.refMSA} -t {input.refTree} -n PARAMS01092023 -w {params.outputDir}
		"""

# FRAGMENT INSERTION
rule runRAXML:
	input:
		seqAlignment=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_plusKarthikSeqs.faa"),
		modelParams=config["binaryModelParams"],
		refTree=config["refTree"]
	params:
		outputDir=config["raxmlOutputDir"],
		fileExtension=config["raxmlOutputFileExtension"]
	output:
		join(config["raxmlOutputDir"], "RAxML_info.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_classification.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_classificationLikelihoodWeights.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_entropy.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_labelledTree.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_originalLabelledTree.raxmlnovel"),
		join(config["raxmlOutputDir"], "RAxML_portableTree.raxmlnovel.jplace")
	shell:
		"""
		raxmlHPC -f v -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w {params.outputDir}
		"""

rule removeBootstrapValues:
	input:
		join(config["raxmlOutputDir"],"RAxML_labelledTree.raxmlnovel")
	output:
		join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_wKarthikSeqs.newick")
	shell:
		"""
		python3 workflow/scripts/removeBootstrapValues.py {input} {output}
		"""

rule getBranchDistances:
	input:
		tree=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_rooted_wKarthikSeqs.newick"),
		query_info=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
		karthik_seq_list="workflow/out/treeInfo/KarthikSeqs_queryIDs_inTree.txt"
	output:
		csv="workflow/out/treeInfo/queryDistanceInfo_wKarthikSeqsAsRef.csv",
		json="workflow/out/treeInfo/queryDistanceInfo_wKarthikSeqsAsRef.json"
	shell:
		"""
		python3 workflow/scripts/getBranchDistances.py {input.tree} {input.query_info} {output.csv} {output.json} {input.karthik_seq_list}
		"""

rule pruneTree:
	input:
		result_tree=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_wKarthikSeqs.newick"),
		ref_tree=join(config["raxmlOutputDir"], "RAxML_labelledTree_noBootstrap_wKarthikSeqs.newick"),
		distance_info="workflow/out/treeInfo/queryDistanceInfo_wKarthikSeqs.csv",
	output:
		pruned_result_tree="workflow/out/treeInfo/prunedResultTree_wKarthikSeqs.nw",
		pruned_ref_tree="workflow/out/treeInfo/prunedRefTree_wKarthikSeqs.nw"
	shell:
		"""
		python3 workflow/scripts/pruneTree.py {input.result_tree} {input.ref_tree} {input.distance_info} {output.pruned_result_tree} {output.pruned_ref_tree}
		"""

rule rootTree:
	input:
		join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_wKarthikSeqs.newick")
	output:
		join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap_rooted_wKarthikSeqs.newick")
	shell:
		"""
		python3 workflow/scripts/rootTree.py {input} {output}
		"""
