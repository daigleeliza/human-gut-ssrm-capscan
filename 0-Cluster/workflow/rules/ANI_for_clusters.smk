rule makeClusterFasta:
	input:
		refFASTA = config["refMSA_nuc"],
		otu_in_cluster_list=join(config["clusterIdsDir"], "{cluster}_otus_in_cluster_015.txt")
		#otu_in_cluster_list=join(config["clusterIdsDir"], "{cluster}_otus_in_cluster.txt")
	output:
		join(config["clusterSeqsDir"], "{cluster}_cluster_seqs_015.fasta")
	conda:
		"../envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		python3 {params.scripts_dir}/selectOTUSeqsInCluster.py {input.refFASTA} {input.otu_in_cluster_list} {output}
		"""

rule makeDB:
	input:
		join(config["clusterSeqsDir"], "{cluster}_cluster_seqs_015.fasta")
	output:
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.ndb"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nhr"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nin"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.njs"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.not"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nsq"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.ntf"),
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nto")
	conda:
		"../envs/blast.yml"
	params:
		join(config["clusterBlastDB"],"{cluster}_db_nuc_015")
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		makeblastdb -in {input} -out {params} -dbtype nucl
		"""

rule runBLAST:
	input:
		query = join(config["clusterSeqsDir"], "{cluster}_cluster_seqs_015.fasta"),
		ndb = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.ndb"),
		nhr = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nhr"),
		nin = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nin"),
		njs = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.njs"),
		not_ = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.not"),
		nsq = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nsq"),
		ntf = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.ntf"),
		nto = join(config["clusterBlastDB"],"{cluster}_db_nuc_015.nto")
		# query=config["primers"]
	output:
		join(config["clusterBlastOutput"],"{cluster}_cluster_blast_015.out")
	conda:
		"../envs/blast.yml"
	params:
		format="10 qseqid sseqid pident",
		db=join(config["clusterBlastDB"],"{cluster}_db_nuc_015")
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		blastn -db {params.db} -query {input.query} -out {output} -max_hsps 1 -max_target_seqs 20 -outfmt '{params.format}'
		"""

rule calcANIperCluster:
	input:
		join(config["clusterBlastOutput"],"{cluster}_cluster_blast_015.out")
	output:
		join(config["clusterANI"],"{cluster}_clusterANI_015.csv")
	conda:
		"../envs/biopython.yml"
	resources:
		mem_mb=4000,
		time=60
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/calcualteANIperCluster.py {input} {wildcards.cluster} {output}
		"""

#rule compile:
#	output:
#		join(config["clusterANI"],"{group}_clusterANI_compiled.csv")
#	params:
#		dir=join(config["clusterANI"],"{group}/")
#	shell:
#		"""
#		cat {params.dir}/*.csv > {output}
#		"""