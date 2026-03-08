rule prefetch:
    output:
        join(config["repositoryDir"],"{accession_num}")
    params:
        acc_num=lambda w: {w.accession_num}
    conda:
        "../../workflow/envs/sra-tools.yml"
    shell:
        "prefetch {params.acc_num} -o {output}"


# dump accession numbers one at a time
rule dump:
    input:
        join(config["repositoryDir"],"{accession_num}")
    output:
        join(config["dumpDir"],"{accession_num}.fa")
    conda:
        "../../workflow/envs/sra-tools.yml"
    shell:
        "vdb-dump -f fasta {input} --output-file {output}"
