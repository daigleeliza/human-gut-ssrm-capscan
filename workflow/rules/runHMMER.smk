rule runProdigal:
    input: join(config["dumpDir"],"{accession_num}.fa")
    output:
        geneCoords=join(config["prodigalGeneCoordDir"],"{accession_num}_geneCoord.out"),
        proteinSeqs=join(config["prodigalProteinSeqDir"], "{accession_num}.faa")
    conda:
        config["prodigalEnv"]
    shell:
        """
        prodigal -i {input} -o {output.geneCoords} -a {output.proteinSeqs}
        """

rule runHMMER:
    input:
        join(config["prodigalProteinSeqDir"],"{accession_num}.faa")
    output:
        hmmOut_bact = join(config["hmmOutDir"],"{accession_num}_dsrAB_bact.hmm.out"),
        domOut_bact = join(config["domtblDir"],"{accession_num}_dsrAB_bact.domtblout"),
        hmmOut_arch = join(config["hmmOutDir"],"{accession_num}_dsrAB_arch.hmm.out"),
        domOut_arch = join(config["domtblDir"],"{accession_num}_dsrAB_arch.domtblout"),
        msa_bact= join(config["msaDir"],"sto/{accession_num}_dsrAB_bact.sto"),
        msa_arch= join(config["msaDir"],"sto/{accession_num}_dsrAB_arch.sto")
    params:
        hmm_profile_bact=config["dsrAB_bact_HMMProfile"],
        hmm_profile_arch=config["dsrAB_arch_HMMProfile"]
    conda:
        config["hmmerEnv"]
    shell:
        """
        hmmsearch -o {output.hmmOut_bact} --domtblout {output.domOut_bact} -A {output.msa_bact} {params.hmm_profile_bact} {input}
        hmmsearch -o {output.hmmOut_arch} --domtblout {output.domOut_arch} -A {output.msa_arch} {params.hmm_profile_arch} {input}
        """

rule parseHMMER:
    input:
        bact=join(config["domtblDir"],"{accession_num}_dsrAB_bact.domtblout"),
        arch=join(config["domtblDir"],"{accession_num}_dsrAB_arch.domtblout")
    output:
        bact=join(config["summaryDir"], "{accession_num}_dsrAB_bact_hits.csv"),
        arch=join(config["summaryDir"], "{accession_num}_dsrAB_arch_hits.csv")
    params:
        scripts_dir=config["scriptsDir"]
    shell:
        """
        python3 {params.scripts_dir}/parse_hmmer_domtable.py {input.bact} {output.bact}
        python3 {params.scripts_dir}/parse_hmmer_domtable.py {input.arch} {output.arch}
        """

# convert the HMMER output MSA from stockholm format (.sto) to fasta (.faa)
rule convertMSA:
    input:
        bact=join(config["msaDir"],"sto/{accession_num}_dsrAB_bact.sto"),
        arch=join(config["msaDir"],"sto/{accession_num}_dsrAB_arch.sto")
    output:
        bact=join(config["msaDir"],"faa/{accession_num}_dsrAB_bact.faa"),
        arch=join(config["msaDir"],"faa/{accession_num}_dsrAB_arch.faa")
    params:
        scripts_dir=config["scriptsDir"]
    shell:
        """
        python3 {params.scripts_dir}/convertFASTA.py {input.bact} {output.bact}
        python3 {params.scripts_dir}/convertFASTA.py {input.arch} {output.arch}
        """

rule combineCSV:
    input:
        bact=expand(join(config["summaryDir"],"{accession_num}_dsrAB_bact_hits.csv"),accession_num=ACCESSION_NUMBERS),
        arch=expand(join(config["summaryDir"],"{accession_num}_dsrAB_arch_hits.csv"),accession_num=ACCESSION_NUMBERS)
    output:
        bact=join(config["summaryDir"],"compiled_dsrAB_bact_hits.csv"),
        arch=join(config["summaryDir"],"compiled_dsrAB_arch_hits.csv")
    params:
        input_dir=config["summaryDir"]
    shell:
        """
        xargs -0 cat {params.input_dir}/*bact_hits.csv  > {output.bact}
        xargs -0 cat {params.input_dir}/*arch_hits.csv  > {output.arch}
        """

rule combineMSA:
    input:
        bact=expand(join(config["msaDir"],"faa/{accession_num}_dsrAB_bact.faa"),accession_num=ACCESSION_NUMBERS),
        arch=expand(join(config["msaDir"],"faa/{accession_num}_dsrAB_arch.faa"), accession_num=ACCESSION_NUMBERS)
    output:
        bact=join(config["msaDir"],"faa/compiled_dsrAB_bact_hits.faa"),
        arch=join(config["msaDir"],"faa/compiled_dsrAB_arch_hits.faa")
    params:
        input_dir=join(config["msaDir"], "faa")
    shell:
        """
        xargs -0 cat {params.input_dir}/*dsrAB_bact.faa  > {output.bact}
        xargs -0 cat {params.input_dir}/*dsrAB_arch.faa  > {output.arch}
        """
