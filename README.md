# SSRM_2022
Snakemake pipeline (cluster folder -> to be downloaded onto cluster and run via slurm submission):
- for NCBI data (Snakefile_NCBI -> copy/paste into Snakefile so Snakemake can recognize the file name)
  1. prepAnantharaman.smk: Remove sequences from Anantharaman that are already found in Mueller.
  2. downloadMetagenomicData.smk: Pull data from NCBI using accession numbers in NCBI_WGS_accessionNumbers.txt and convert to .faa.
  3. runHMMER.smk: Run prodigal to predict genes then run hmmer search to find dsrAB genes. convert output to fasta and csv formats.
  4. cleanHits.smk: Filter fasta and csv files by bit score (100), remove duplicate hmmer hits, and run multiple sequence alignment between the fasta and the reference sequences. 
  5. fragmentInsertion.smk: Use RAxML to estimate model parameters and insert query sequences into the reference tree.
  6. getTreeInfo.smk: Get identity of closest reference leaves to our query sequences.
- for Capscan/Stool data (Snakefile_Capsule): same steps as rules above, with file names changed.
  1. runHMMERCapsule.smk
  2. cleanHitsCapsule.smk
  3. fragmentInsertionCapsule.smk
  4. getTreeInfoCapsule.smk
- for gene abundances
  1. mapReads.smk: map raw reads to capsule/stool coassemblies and get reads per contig
  2. perGeneCounts.smk: get reads per gene and calculate RPKMs
  3. getGeneAbunCapsule.smk: extract RPKMs of dsrAB genes to make dataset for analysis
 
Move into RStudio for data visualization
  1. Ensure packages listed in config are installed
  2. Data management
  3. Analysis
