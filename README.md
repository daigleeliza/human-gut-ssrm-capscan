# SSRM_2022
Snakemake pipeline:
- for NCBI data (Snakefile_NCBI -> copy/paste into Snakefile so Snakemake can recognize the file name)
  1. prepAnantharaman.smk: remove sequences from Anantharaman that are already found in Mueller
  2. downloadMetagenomicData.smk: pull data from NCBI using accession numbers in NCBI_WGS_accessionNumbers.txt and convert to .faa
  3. runHMMER.smk
  4. cleanHits.smk
  5. fragmentInsertion.smk
  6. getTreeInfo.smk
- for Capscan/Stool data (Snakefile_Capsule)
  1. runHMMERCapsule.smk
  2. cleanHitsCapsule.smk
  3. fragmentInsertionCapsule.smk
  4. getTreeInfoCapsule.smk
- for gene abundances
  1. mapReads.smk
  2. perGeneCounts.smk
  3. getGeneAbunCapsule.smk
