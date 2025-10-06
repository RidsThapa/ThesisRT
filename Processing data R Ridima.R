# Libraries
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# Save TnSeq data
wig_file <- "combined_wig.txt"
wig_df <- readr::read_tsv(wig_file, comment = "#", col_names = TRUE)

# Column naming for only relevant information
colnames(wig_df)[1] <- "position"
colnames(wig_df)[ncol(wig_df)] <- "gene_id" 
rep_cols <- paste0("rep", 1:(ncol(wig_df) - 2))
colnames(wig_df) <- c("position", rep_cols, "gene_id")

# Turn replicate columns into a total insertions
tnseq_counts <- wig_df %>%
  dplyr::mutate(insertions = rowSums(dplyr::across(starts_with("rep")))) %>%
  dplyr::select(gene_id, position, insertions)

# Clean gene identifiers to get gene symbol
tnseq_counts <- tnseq_counts %>%
  dplyr::mutate(
    gene_id = stringr::str_extract(gene_id, "(?<=\\().+?(?=\\))"),
    gene_id = ifelse(is.na(gene_id), gene_id, gene_id),
    gene_id = tolower(gene_id),
    essentiality = "unknown"
  )

readr::write_csv(tnseq_counts, "tnseq_clean.csv")

# Save domain data
domain_file <- "uniprot-proteome_UP000001584.tsv"
dom_df <- readr::read_tsv(domain_file)

# Column naming for only relevant information
dom_df <- dom_df %>%
  dplyr::rename(
    gene_id     = `Gene Names`,
    domain_name = Pfam,
    domain_info = `Domain [FT]`
  )

# Change aliases + Pfam IDs into rows
dom_df <- dom_df %>%
  dplyr::mutate(gene_id = tolower(gene_id)) %>%
  tidyr::separate_rows(gene_id, sep = " ") %>%
  dplyr::filter(gene_id != "")
dom_df <- dom_df %>%
  tidyr::separate_rows(domain_name, sep = ";") %>%
  dplyr::mutate(domain_name = trimws(domain_name)) %>%
  dplyr::filter(domain_name != "")

# Extract start/end coords from Domain ft
dom_df <- dom_df %>%
  dplyr::mutate(
    domain_start = as.numeric(stringr::str_extract(domain_info, "\\d+(?=\\.\\.)")),
    domain_end   = as.numeric(stringr::str_extract(domain_info, "(?<=\\.\\.)\\d+"))
  ) %>%
  dplyr::select(gene_id, domain_name, domain_start, domain_end) %>%
  dplyr::distinct()

readr::write_csv(dom_df, "domains_final.csv")


# Check if data overlaps
tnseq <- readr::read_csv("tnseq_clean.csv", show_col_types = FALSE)
domains <- readr::read_csv("domains_final.csv", show_col_types = FALSE)

common_genes <- intersect(na.omit(unique(tnseq$gene_id)),
                          na.omit(unique(domains$gene_id)))
