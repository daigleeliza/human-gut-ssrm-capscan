rule makeTaxaTable:
	input: "config/Mueller2015_dsrABref_wTaxa.faa"
	output: "workflow/out/clustering_otus/taxa_table_0420.csv"
	conda:
		"../envs/biopython.yml"
	threads: 1
	resources: 
		time=30,
		mem_mb=4000
	shell:
		"""
		python3 workflow/scripts/make_taxonTable.py {input} {output}
		"""
