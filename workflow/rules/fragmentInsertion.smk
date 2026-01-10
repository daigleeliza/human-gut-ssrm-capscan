# identify novel sequences in Anantharaman2018 (not already present in Mueller2015 reference sequences)
rule IdentifyNovelSeqs:
    input:
        Anantharaman2018Seqs="config/Anantharaman2018_dsrA_dsrB.faa",
        Mueller2015Seqs="config/Mueller2015_dsrrAB.faa"
    params:
        outputDir=config["cleanHitsDir"]
    output:
        no_exact_dups_f=join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noEDups.faa"),
        no_short_dups_f=join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa"),
        dups_info=join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.json")
    shell:
        """
        python3 workflow/scripts/Anantharaman2018_novel_seqs.py {input.Anantharaman2018Seqs} {input.Mueller2015Seqs} {params.outputDir}
        """

# make MSA from our result msa (with hits and ref seqs) with Anantharaman2018's novel seqs
rule runMAFFT_2:
    input:
        hitsWithRefMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
        Anantharaman2018Seqs=join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa")
    output:
        resultMSA=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa")
    shell:
        """
        /usr/local/bin/mafft --add {input.Anantharaman2018Seqs} --keeplength {input.hitsWithRefMSA} > {output.resultMSA}
        """

rule estimateModelParams:
    input:
        refMSA=config["refMSA"],
        refTree="config/dsrAB_consensus_phylogeny.newick"
    params:
        outputDir=config["modelParamsDir"],
        treeFileExtension=config["treeFileExtension"]
    shell:
        """
        raxmlHPC -f e -m PROTGAMMADAYHOFF -s {input.refMSA} -t {input.refTree} -n {params.treeFileExtension} -w {params.outputDir}
        """

# FRAGMENT INSERTION
rule runRAXML:
    input:
        seqAlignment=join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa"),
        modelParams=join(join(config["modelParamsDir"], "RAxML_binaryModelParameters."), config["treeFileExtension"]),
        refTree=join(join(config["modelParamsDir"], "RAxML_result."), config["treeFileExtension"])
    params:
        outputDir=config["raxmlOutputDir"],
        fileExtension=config["treeFileExtension"]
    output:
        join(join(config["raxmlOutputDir"], "RAxML_info."), config["treeFileExtension"]),
        join(join(config["raxmlOutputDir"], "RAxML_classification."), config["treeFileExtension"]),
        join(join(config["raxmlOutputDir"], "RAxML_classificationLikelihoodWeights."), config["treeFileExtension"]),
        join(join(config["raxmlOutputDir"], "RAxML_entropy."), config["treeFileExtension"]),
        join(join(config["raxmlOutputDir"], "RAxML_labelledTree."), config["treeFileExtension"]),
        join(join(config["raxmlOutputDir"], "RAxML_originalLabelledTree."), config["treeFileExtension"]),
        join(join(join(config["raxmlOutputDir"], "RAxML_portableTree.", config["treeFileExtension"])), ".jplace")
    shell:
        """
        raxmlHPC -f v -R {input.modelParams} -r {input.refTree} -s {input.seqAlignment} -m PROTGAMMADAYHOFF -G 0.1 -n {params.fileExtension} -w {params.outputDir}
        """

# CLEAN TREE
rule removeBootstrapValues:
    input:
        join(join(config["raxmlOutputDir"],"RAxML_labelledTree."), config["treeFileExtension"])
    output:
        join(config["raxmlOutputDir"],"RAxML_labelledTree_noBootstrap.newick")
    shell:
        """
        python3 workflow/scripts/removeBootstrapValues.py {input} {output}
        """
