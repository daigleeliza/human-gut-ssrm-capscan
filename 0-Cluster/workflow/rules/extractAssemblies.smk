# extract query IDs from csv file of biosamples included
rule extractQueryIDs:
	input: assemblies = config["assemblies"]
	output: queryIDs = config["assembliesQueryIDs"]
	conda:
		"../../workflow/envs/R_minimal.yml"
	script:
		"../scripts/extractQueryIDs.R"



# filter compiled csv files for queryIDs included in our list of assemblies
rule extractAssembliesCSV:
	input:
		queryIDs = config["assembliesQueryIDs"],
		bact=join(config["summaryDir"],"compiled_dsrAB_bact_hits.csv"),
		arch=join(config["summaryDir"],"compiled_dsrAB_arch_hits.csv")
	output: 
		filterBact=join(config["summaryDir"],"compiled_dsrAB_bact_hits_assemblies.csv"),
		filterArch=join(config["summaryDir"],"compiled_dsrAB_arch_hits_assemblies.csv")
	resources:
		time=30,
		mem_mb=4000
	conda:
		"../../workflow/envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/filterAssemblies_csv.py {input.queryIDs} {input.bact} {input.arch} {output.filterBact} {output.filterArch}
		"""

# filter compiled faa files for biosamples included in our list of assemblies
rule extractAssembliesMSA:
	input:
		queryIDs = config["assembliesQueryIDs"],
		bact=join(config["summaryDir"],"compiled_dsrAB_bact_hits.faa"),
		arch=join(config["summaryDir"],"compiled_dsrAB_arch_hits.faa")
	output:
		filterBact=join(config["summaryDir"],"compiled_dsrAB_bact_hits_assemblies.faa"),
		filterArch=join(config["summaryDir"],"compiled_dsrAB_arch_hits_assemblies.faa")
	resources:
		time=30,
		mem_mb=4000
	conda:
		"../../workflow/envs/biopython.yml"
	params:
		scripts_dir=config["scriptsDir"]
	shell:
		"""
		python3 {params.scripts_dir}/filterAssemblies_faa.py {input.queryIDs} {input.bact} {input.arch} {output.filterBact} {output.filterArch}
		"""