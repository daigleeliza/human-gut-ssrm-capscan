srb_lite_notes_ages_countries <- read.csv(snakemake@input[["assemblies"]], stringsAsFactors = FALSE)

all_ids <- unlist(strsplit(srb_lite_notes_ages_countries$query_ids, ","))
all_ids <- trimws(all_ids)  # remove whitespace
all_ids <- gsub("\\[", "", all_ids)
all_ids <- gsub("\\]", "", all_ids)
all_ids <- gsub("'", "", all_ids)  # remove brackets and quotes

write(all_ids, snakemake@output[["queryIDs"]])