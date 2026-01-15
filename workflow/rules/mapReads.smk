# build bowtie2 index files from co-assemblies
# set sample = samples, coassembly = coassemblies
rule buildIndexFiles_coassemblies:
	input: join(config["coAssemblyDir"],"{coassembly}/megahit/filtered_contigs_1kb.fasta")
	output: join(config["coAssemblyDir"],"{coassembly}/bowtie2_index.1.bt2")
	threads: 8
	params:
		base=join(config["coAssemblyDir"],"{coassembly}/bowtie2_index")
	conda:
		"../envs/bowtie2-no-builds.yml"
	shell:
		"""
		bowtie2-build {input} {params.base} --seed 2525 --threads {threads} --offrate 3
		"""

# Use bowtie2 to map raw reads to coassemblies.
# set sample = samples, set coassembly = coassemblies
rule mapReads:
	input:
		filt1=join(config["filterDir"],"{sample}-filtered.1.fastq.gz"),
		filt2=join(config["filterDir"],"{sample}-filtered.2.fastq.gz"),
		index=join(config["coAssemblyDir"], "{coassembly}/bowtie2_index.1.bt2")
	output:
		bt2log=join(config["SAMDir"],"CoAssembly/{coassembly}/{sample}.bt2.log"),
		samfile=join(config["SAMDir"],"CoAssembly/{coassembly}/{sample}.sam")
	threads: 8
	params:
		#coassembly=lambda wildcards: df.query('samplename == @wildcards.sample')['co_assembly'].item(),
		base=join(config["coAssemblyDir"],"{coassembly}/bowtie2_index"),
		SAMDir=config['SAMDir']
	conda:
		"../envs/bowtie2-no-builds.yml"
	shell:
		"""
		mkdir -p {params.SAMDir}/CoAssembly/{wildcards.coassembly}
		bowtie2 --local -x {params.base} -1 {input.filt1} -2 {input.filt2} -S {output.samfile} --threads {threads} --sensitive-local --reorder -t 2> {output.bt2log}
		"""

# convert SAM to BAM, remove unmapped reads; get bam index
# set sample = samples, set coassembly = coassemblies
rule mappedReadsBAM:
	input:
		sam=join(config["SAMDir"],"CoAssembly/{coassembly}/{sample}.sam")
	output:
		bam=join(config["mapDir"],"bam/CoAssembly/{coassembly}/{sample}.bam"),
		bam_index=join(config["mapDir"],"bam/CoAssembly_index/{coassembly}/{sample}.bam.bai")
	threads: 4
	resources:
		mem_mb=16000,
		time=240
	params:
		mapDir=config['mapDir'],
		tmpDir=config['tmpDir']
	conda:
		"../envs/bbmap.yml"
	shell:
		"""
		samtools view -b -F 4 {input.sam} | samtools sort -@ {threads} -T {params.tmpDir} -o {output.bam}
		samtools index -@ {threads} {output.bam} {output.bam_index}
		"""
		#mkdir -p {params.mapDir}/bam/CoAssembly/{wildcards.coassembly}

# get read counts per contig for each sample
# set sample = samples, set coassembly = coassemblies
rule contigReads:
	input:
		bam=join(config["mapDir"],"bam/CoAssembly/{coassembly}/{sample}.bam"),
		bam_index=join(config["mapDir"],"bam/CoAssembly_index/{coassembly}/{sample}.bam.bai")
	output:
		join(config["mapDir"],"bam/CoAssembly_contig_reads/{coassembly}/{sample}_reads_per_contig.txt")
	threads: 1
	resources:
		mem_mb=16000,
		time=240
	conda:
		"../envs/bbmap.yml"
	shell:
		"""
		samtools idxstats {input.bam} > {output}
		"""
