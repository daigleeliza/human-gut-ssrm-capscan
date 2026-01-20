rule prepAnantharaman:
	input: 
		Anantharaman2018="config/Anantharaman2018_dsrA_dsrB.faa",
		Mueller2015="config/Mueller2015_dsrAB.faa"
	output: 
		noDups_fasta="config/Anantharaman2018_dsrA_dsrB_noDups.faa",
		noEDups_fasta="config/Anantharaman2018_dsrA_dsrB_noEDups.faa",
		noDups_json="config/Anantharaman2018_dsrA_dsrB_noDups.json"
	conda: "../envs/biopython.yml"
	threads: 1
	resources:
		time=30,
		mem_mb=4000
	shell:
		"""
		python3 workflow/scripts/Anantharaman2018_novel_seqs.py {input.Anantharaman2018} {input.Mueller2015} {output.noDups_fasta} {output.noEDups_fasta} {output.noDups_json}
		"""
