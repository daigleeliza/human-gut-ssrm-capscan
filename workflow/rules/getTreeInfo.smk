rule getAnantharamanQueriesIDs:
    input: join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa")
    output: "workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
    shell:
        "python3 Ananthtaraman2018_querySeqsInTree.py {input} {output}"


rule getBranchDistances:
    input:
        tree=join(config["raxmlOutputDir"],"RAxML_originalLabelledTree_noBootstrap.newick"),
        query_info=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
        Anantharaman2018_seq_list="workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt"
    output:
        csv="workflow/out/treeInfo/queryDistanceInfo.csv",
        json="workflow/out/treeInfo/queryDistanceInfo.json"
    shell:
        """
        python3 workflow/scripts/getBranchDistances.py {input.tree} {input.query_info} {input.Anantharaman2018_seq_list} {output.csv} {output.json}
        """

rule compileRefInfo:
    input:
        queryDistanceCSV="workflow/out/treeInfo/queryDistanceInfo.csv",
        scoreThresholdCSV=join(config["cleanHitsDir"], "compiled_dsrAB_hits_scoreThreshold.csv"),
        biosampleJSON=config["biosampleMetadata"]
    output:
        "workflow/out/treeInfo/closestRefInfo_allScoreThresholdHits.csv"
    shell:
        """
        python3 workflow/scripts/compileRefInfo.py {input.queryDistanceCSV} {input.scoreThresholdCSV} {input.biosampleJSON} {output}
        """

rule listClosestRefLeaves:
    input:
        queryDistanceJSON="workflow/out/treeInfo/queryDistanceInfo.csv",
        AnantharamanNoDups=join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa"),
    output: "workflow/out/treeInfo/closestRefList.txt"
    shell:
        """
        python3 workflow/scripts/listClosestRefs.py {input.queryDistanceJSON} {input.AnantharamanNoDups} {output}
        """
