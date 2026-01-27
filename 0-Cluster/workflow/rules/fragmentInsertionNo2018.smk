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

## get tree treeInfo

rule getBranchDistances:
	input:
		tree=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_no_2018_capsule.newick"),
		query_info=join(config["cleanHitsDir"],"StoolCapsule_arch_bact_compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
		Anantharaman2018_seq_list="workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
	output:
		csv="workflow/out/treeInfo/queryDistanceInfo_no_2018_capsule.csv",
		json="workflow/out/treeInfo/queryDistanceInfo_no_2018_capsule.json"
	conda:
		"../envs/biopython.yml"
	threads: 1
	resources:
		time=400,
		mem_mb=2000
	shell:
		"""
		python3 workflow/scripts/getBranchDistances.py {input.tree} {input.query_info} {input.Anantharaman2018_seq_list} {output.csv} {output.json}
		"""

rule compileRefInfo:
	input:
		queryDistanceCSV="workflow/out/treeInfo/queryDistanceInfo_no_2018_capsule.csv",
		scoreThresholdCSV=join(config["cleanHitsDir"], "StoolCapsule_arch_bact_compiled_dsrAB_hits_scoreThreshold.csv"),
	output:
		"workflow/out/treeInfo/closestRefInfo_allScoreThresholdHits_no_2018_capsule.csv"
	conda:
		"../envs/biopython.yml"
	threads: 1
	resources:
		time=30,
		mem_mb=1000
	shell:
		"""
		python3 workflow/scripts/compileRefInfoCapsule.py {input.queryDistanceCSV} {input.scoreThresholdCSV} {output}
		"""

rule listClosestRefLeaves:
	input:
		queryDistanceJSON="workflow/out/treeInfo/queryDistanceInfo_no_2018_capsule.json",
		AnantharamanNoDups="config/Anantharaman2018_dsrA_dsrB_noDups.faa"
	output: "workflow/out/treeInfo/closestRefList_no_2018_capsule.txt"
	conda:
		"../envs/biopython.yml"
	threads: 1
	resources:
		time=5,
		mem_mb=500
	shell:
		"""
		python3 workflow/scripts/listClosestRefs.py {input.queryDistanceJSON} {input.AnantharamanNoDups} {output}
		"""
