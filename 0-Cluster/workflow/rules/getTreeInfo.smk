rule getAnantharamanQueriesIDs:
	input: "config/Anantharaman2018_dsrA_dsrB_noDups.faa"
	output: "workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
	conda:
		"../envs/biopython.yml"
	resources:
		time=5,
		mem_mb=500
	shell:
		"""
		python3 workflow/scripts/Anantharaman2018_querySeqsInTree.py {input} {output}
		"""

# find distance to nearest reference leaf
rule getBranchDistances:
	input:
		tree=join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap_no_2018_assemblies.newick"),
		query_info=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_assemblies_noDups_gapPercentageInfo.csv"),
		Anantharaman2018_seq_list="workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
	output:
		csv="workflow/out/treeInfo/queryDistanceInfo_no_2018_assemblies.csv",
		json="workflow/out/treeInfo/queryDistanceInfo_no_2018_assemblies.json"
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

# adds closest reference leaf (by branch distance) and sample ID to all hit sequences from dsrAB HMMER search that passed the score threshold (so, all duplicate sequences)
rule compileRefInfo:
	input:
		queryDistanceCSV="workflow/out/treeInfo/queryDistanceInfo_no_2018_assemblies.csv",
		scoreThresholdCSV=join(config["cleanHitsDir"], "compiled_dsrAB_hits_scoreThreshold_assemblies.csv"),
		biosampleJSON=config["biosampleMetadata"]
	output:
		"workflow/out/treeInfo/closestRefInfo_allScoreThresholdHits_no_2018_assemblies.csv"
	conda:
		"../envs/biopython.yml"
	threads: 1
	resources:
		time=30,
		mem_mb=1000
	shell:
		"""
		python3 workflow/scripts/compileRefInfo.py {input.queryDistanceCSV} {input.scoreThresholdCSV} {input.biosampleJSON} {output}
		"""

# make list of closest reference leaves found in this dataset
rule listClosestRefLeaves:
	input:
		queryDistanceJSON="workflow/out/treeInfo/queryDistanceInfo_no_2018_assemblies.json",
		AnantharamanNoDups="config/Anantharaman2018_dsrA_dsrB_noDups.faa"
	output: "workflow/out/treeInfo/closestRefList_no_2018_assemblies.txt"
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
