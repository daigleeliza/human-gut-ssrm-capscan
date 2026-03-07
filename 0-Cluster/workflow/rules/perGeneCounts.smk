# convert GFF to GTF
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
# get reads per gene
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

# get average genome size and genome equivalents for each sample
#set sample=samples
rule runMicrobeCensus:
	input: 
		filt1=join(config["filterDir"],"{sample}-filtered.1.fastq.gz"),
		filt2=join(config["filterDir"],"{sample}-filtered.2.fastq.gz")
	output: "workflow/out/microbeCensus/{sample}.txt"
	conda: "../envs/microbeCensus.yml"
	params:
		tmpDir = config['tmpDir']
	threads: 4
	resources: 
		mem_mb=150000,
		time=240
	shell:
		"""
		export TMPDIR={params.tmpDir}
		run_microbe_census.py -t {threads} {input.filt1},{input.filt2} {output}
		"""

# get list of all the genome equivalents for the samples 
rule getGenomeEquivalents:
	input: expand(join("workflow/out/microbeCensus/{sample}.txt"), sample=[s for c in coassemblies for s in get_subject_sample_list_dropped(c)])
	output: "workflow/out/microbeCensus/genome_equivalents.txt"
	resources:
		mem_mb = 16000,
		time = 30
	shell:
		"""
		echo -e "sample\\tgenome_equivalents" > {output}
		for f in {input}; do
			sample=$(basename $f .txt)
			genomeEquivalents=$(grep "genome_equivalents" $f | awk '{{print $2}}')
			echo -e "$sample\t$genomeEquivalents" >> {output}
		done
		"""

# calculate RPKG
rule RPKG:
	input:
		geneReads=join(config["mapDir"],"bam/CoAssembly_gene_reads/{coassembly}/{sample}_reads_per_gene.txt"),
		genomeEquivalents="workflow/out/microbeCensus/genome_equivalents.txt"
	output:
		geneReadsRPKG=join(config["mapDir"],"bam/CoAssembly_gene_reads/{coassembly}/{sample}_microbeCensus_RPKG.txt")
	threads: 1
	resources:
		mem_mb=4000,
		time=60
	shell:
		"""
		genomeEquivalents=$(grep "{wildcards.sample}" {input.genomeEquivalents} | awk '{{print $2}}')
		awk -v genEq="$genomeEquivalents" 'BEGIN{{OFS="\\t"}}
		NR == 1 {{print; next}}
		NR == 2 {{print $0, "RPKG"; next}}
		{{
			rpkg = ($7 * 1e3) / ($6 * genEq)
			print $0, rpkg
		}}' {input.geneReads} > {output.geneReadsRPKG}
		"""
