# SSRM_2026
Snakemake pipeline (cluster folder -> to be downloaded onto cluster and run via slurm submission; change name of "Snakefile_NCBI/Capsule" to "Snakefile" so it is recognizable):
- for NCBI data (Snakefile_NCBI)
  1. prepAnantharaman.smk: Remove sequences from Anantharaman that are already found in Mueller.
  2. downloadMetagenomicData.smk: Pull data from NCBI using accession numbers in NCBI_WGS_accessionNumbers.txt and convert to .faa.
  3. extractAssemblies.smk: Selects NCBI data that are gut metagenomic assemblies.
  4. runHMMER.smk: Run prodigal to predict genes then run hmmer search to find dsrAB genes. convert output to fasta and csv formats.
  5. cleanHits.smk: Filter fasta and csv files by bit score (100), find/remove duplicate HMMER hits, and run multiple sequence alignment between the fasta and the reference sequences. 
  6. fragmentInsertion.smk: Use RAxML to estimate model parameters and insert query sequences into the reference tree.
  7. getTreeInfo.smk: Find distances from queries to closest reference leaves.
  8. makeTaxaTable.smk: Makes a taxa table for making a phyloseq object and identifying reference sequences.
- for Capscan/Stool data (Snakefile_Capsule): same steps as rules above, with file names changed.
  1. runHMMERCapsule.smk
  2. cleanHitsCapsule.smk
  3. fragmentInsertionCapsule.smk
  4. getTreeInfoCapsule.smk
- for gene abundances
  1. mapReads.smk: map raw reads to capsule/stool coassemblies and get reads per contig
  2. perGeneCounts.smk: get reads per gene and calculate RPKGs
  3. getGeneAbunCapsule.smk: extract RPKGs of dsrAB genes to make dataset for analysis
  4. ASRhmmer.smk: run hmmer search for asrABC genes
  5. ASRgeneabun.smk: extract RPKGs of asrABC genes
 
Move into RStudio for data visualization
  1. Ensure packages listed in config are installed
  2. Data management
  3. Analysis
     -> use ANI_for_clusters.smk to find ANI within clusters made in clustering_cooccur_prevalences.Rmd
