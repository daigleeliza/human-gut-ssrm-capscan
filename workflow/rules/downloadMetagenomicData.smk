rule prefetch:
    output:
        join(config["repositoryDir"],"{accession_num}")
    params:
        acc_num=lambda w: {w.accession_num}
    conda:
        config["sraEnv"]
    shell:
        "prefetch {params.acc_num} -o {output}"


# dump accession numbers one at a time
rule dump:
    input:
        join(config["repositoryDir"],"{accession_num}")
    output:
        join(config["dumpDir"],"{accession_num}.fa")
    conda:
        config["sraEnv"]
    shell:
        "vdb-dump -f fasta {input} --output-file {output}"
