from os.path import join

configfile: "config/config.yaml"

# get list of accession numbers from accession number file
with open(config["accessionNumFile"], 'r') as f:
    accession_num_file = f.read()
    ACCESSION_NUMBERS = accession_num_file.split()
CAPSULE_SAMPLES = {'Capsule_1', 'Capsule_2', 'Capsule_3', 'Capsule_4', 'Capsule_5', 'Capsule_6', 'Capsule_7', 'Capsule_8', 'Capsule_9', 'Capsule_10', 'Capsule_11', 'Capsule_12', 'Capsule_13', 'Capsule_14', 'Capsule_15', 'Stool_1', 'Stool_2', 'Stool_3', 'Stool_4', 'Stool_5', 'Stool_6', 'Stool_7', 'Stool_8', 'Stool_9', 'Stool_10', 'Stool_11', 'Stool_12', 'Stool_13', 'Stool_14', 'Stool_15'}

rule all:
    input:
        # download human gut metagenomic data from NCBI
#        expand(join(config["repositoryDir"],"{accession_num}"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["dumpDir"],"{accession_num}.fa"),accession_num=ACCESSION_NUMBERS),
        # run HMMER search of dsrAB on human gut metagenomic data
#        expand(join(config["prodigalProteinSeqDir"], "{accession_num}.faa"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["hmmOutDir"],"{accession_num}_dsrAB_bact.hmm.out"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["domtblDir"],"{accession_num}_dsrAB_bact.domtblout"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["hmmOutDir"],"{accession_num}_dsrAB_arch.hmm.out"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["domtblDir"],"{accession_num}_dsrAB_arch.domtblout"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["summaryDir"],"{accession_num}_dsrAB_bact_hits.csv"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["summaryDir"],"{accession_num}_dsrAB_arch_hits.csv"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["msaDir"],"sto/{accession_num}_dsrAB_bact.sto"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["msaDir"],"sto/{accession_num}_dsrAB_arch.sto"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["msaDir"],"faa/{accession_num}_dsrAB_bact.faa"), accession_num=ACCESSION_NUMBERS),
#        expand(join(config["msaDir"],"faa/{accession_num}_dsrAB_arch.faa"), accession_num=ACCESSION_NUMBERS),
        # clean hit sequences returned by HMMER search
        expand(join(config["hmmOutDir"],"{capsule_sample}_dsrAB_bact.hmm.out"), capsule_sample = CAPSULE_SAMPLES),
#        expand(join(config["summaryDir"],"{accession_num}_dsrAB_bact_hits.csv"), accession_num = CAPSULE_SAMPLES),
        # clean hit sequences returned by HMMER search
#        join(config["summaryDir"],"compiled_dsrAB_bact_hits.csv"),
#        join(config["summaryDir"],"compiled_dsrAB_arch_hits.csv"),
#        join(config["msaDir"],"faa/compiled_dsrAB_bact_hits.faa"),
#        join(config["msaDir"],"faa/compiled_dsrAB_arch_hits.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.csv"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.csv"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_bact_hits_scoreThreshold.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_arch_hits_scoreThreshold.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold.csv"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_noDups.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_noDups.json"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_hits_scoreThreshold_noEDups.faa"),
#        join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_noDups_msa_withRef.faa"),
#        join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_noDups_msa_noRef.faa"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_gapPercentageInfo.csv"),
#        join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_noRef_trimmedGaps.faa"),
#        join(config["cleanHitsDir"], "compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps.faa"),
#        # fragment insertion of hit sequences passing quality filter on to reference dsrAB tree
#        join(config["cleanHitsDir"],"compiled_dsrAB_scoreThreshold_noDups_msa_withRef_trimmedGaps_withAnantharaman2018Seqs.faa"),
#        join(join(config["modelParamsDir"],"RAxML_binaryModelParameters."),config["treeFileExtension"]),
#        join(join(config["modelParamsDir"],"RAxML_result."),config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_info."), config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_classification."), config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_classificationLikelihoodWeights."), config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_entropy."), config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_labelledTree."), config["treeFileExtension"]),
#        join(join(config["raxmlOutputDir"],"RAxML_originalLabelledTree."), config["treeFileExtension"]),
#        join(join(join(config["raxmlOutputDir"],"RAxML_portableTree.",config["treeFileExtension"])), ".jplace"),
#        # calculate branch distances between hit sequences and reference sequences and other tree info
#        "workflow/out/treeInfo/Anantharaman2018_queryIDs_inTree.txt",
#        "workflow/out/treeInfo/queryDistanceInfo.csv",
#        "workflow/out/treeInfo/queryDistanceInfo.json",
#        "workflow/out/treeInfo/closestRefInfo_allScoreThresholdHits.csv",
#        "workflow/out/treeInfo/closestRefList.txt"

#include:  "workflow/rules/downloadMetagenomicData.smk"
include:  "workflow/rules/runHMMER.smk"
#include:  "workflow/rules/cleanHits.smk"
#include:  "workflow/rules/fragmentInsertion.smk"
#include:  "workflow/rules/getTreeInfo.smk"
