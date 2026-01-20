rule makeGTF:
	input:
		gff=join(config["prodigalDir"],"coassemblies/{coassembly}_contigs_prodigal.out")
	output:
		gtf=join(config["prodigalDir"],"coassemblies/{coassembly}_contigs_prodigal_GTF.gtf")
	threads: 1
	resources:
		mem_mb=1000,
		time=150
	conda:
		"../envs/featureCounts.yml"
	shell:
		#convert to GTF for use with featureCounts, using gffread, and get reads
		"""
		gffread {input.gff} -T -o {output.gtf}
		sed -i 's/transcript_id/gene_id/g' {output.gtf}
		"""

# The output files from the following rules are in the GitHub. Prodigal output, GTF, and bam files can be accessed by contacting authors.
rule perGeneCounts:
	input:
		gtf=join(config["prodigalDir"],"coassemblies/{coassembly}_contigs_prodigal_GTF.gtf"),
		bam=join(config["mapDir"],"bam/CoAssembly/{coassembly}/{sample}.bam")
	output:
		geneReads=join(config["mapDir"],"bam/CoAssembly_gene_reads/{coassembly}/{sample}_reads_per_gene.txt")
	threads: 8
	resources:
		mem_mb=8000,
		time=150
	conda:
		"../envs/featureCounts.yml"
	shell:
		#convert to GTF for use with featureCounts, using gffread, and get reads
		"""
		featureCounts -f -p --countReadPairs -t CDS -B -a {input.gtf} -o {output.geneReads} {input.bam} --largestOverlap -T {threads}
		"""

rule RPKM:
	input:
		geneReads=join(config["mapDir"],"bam/CoAssembly_gene_reads/{coassembly}/{sample}_reads_per_gene.txt")
	output:
		geneReadsRPKM=join(config["mapDir"],"bam/CoAssembly_gene_reads/{coassembly}/{sample}_RPKM.txt")
	threads: 1
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		totalMapped=$(awk 'NR > 2 {{sum += $7}} END {{print sum}}' {input.geneReads})
		awk -v total="$totalMapped" 'BEGIN{{OFS="\\t"}}
		NR == 1 {{print; next}}
		NR == 2 {{print $0, "RPKM"; next}}
		{{
			rpkm = ($7 * 1e9) / ($6 * total)
			print $0, rpkm
		}}' {input.geneReads} > {output.geneReadsRPKM}
		"""
