rule prepAnantharaman:
	input: 
		Anantharaman2018="config/Anantharaman2018_dsrA_dsrB.faa",
		Mueller2015="config/Mueller2015_dsrAB.faa"
	output: join(config["cleanHitsDir"], "Anantharaman2018_dsrA_dsrB_noDups.faa")
	threads: 1
	resources:
		time=60,
		mem_mb=4000
	shell:
		"""
		python3 workflow/scripts/Anantharaman2018_novel_seqs.py {input.Anantharaman2018} {input.Mueller2015} {output}
		"""
