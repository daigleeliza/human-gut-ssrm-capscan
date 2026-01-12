# SSRM_2022
Snakemake pipeline:
  1. prepAnantharaman.smk: remove sequences from Anantharaman that are already found in Mueller
  2. downloadMetagenomicData.smk: pull data from NCBI using accession numbers in NCBI_WGS_accessionNumbers.txt and convert to .faa
  3. runHMMER.smk
  4. cleanHits.smk
  5. fragmentInsertion.smk
  6. getTreeInfo.smk
