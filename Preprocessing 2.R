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

# Replicate columns into total insertions
tnseq_counts <- wig_df %>%
  dplyr::mutate(insertions = rowSums(dplyr::across(starts_with("rep")))
  ) %>%
  dplyr::select(gene_id, position, insertions)

# Clean gene identifiers
tnseq_counts <- tnseq_counts %>%
  dplyr::mutate(
    gene_id = ifelse(
      stringr::str_detect(gene_id, "\\("),
      stringr::str_extract(gene_id, "(?<=\\().+?(?=\\))"),
      gene_id),
    
    gene_id = tolower(trimws(gene_id)),
    position = as.numeric(position),
    insertions = as.numeric(insertions)
    
  ) %>%
  dplyr::filter(!is.na(gene_id), gene_id != "", !is.na(position), !is.na(insertions))

# Save domain data
domain_file <- "uniprot-proteome_UP000001584.tsv"

dom_df <- readr::read_tsv(domain_file, show_col_types = FALSE)

# Column naming for only relevant information
dom_df <- dom_df %>%
  dplyr::rename(gene_names = `Gene Names`, domain_info = `Domain [FT]`, protein_length = Length)

# Keep UniProt protein information
dom_df <- dom_df %>%
  dplyr::select(Entry, gene_names, protein_length, domain_info)

# Load H37Rv genome annotation
gff_file <- "H37Rv.gff"

gff <- readr::read_tsv(gff_file, comment = "#", col_names = FALSE, show_col_types = FALSE)

# Name GFF columns
colnames(gff) <- c(
  "seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Keep protein coding genes
gff_genes <- gff %>%
  dplyr::filter(type == "CDS", seqid == "NC_000962.3")


# Extract gene names and Rv locus tags
gff_genes <- gff_genes %>%
  dplyr::mutate(
    gene_name = stringr::str_extract(
      attributes,
      "(?<=gene=)[^;]+"),
    
    locus_tag = stringr::str_extract(
      attributes,
      "(?<=locus_tag=)[^;]+"),
    
    gene_name = tolower(gene_name), locus_tag = tolower(locus_tag)
  )


# Keep relevant genomic information
gff_genes <- gff_genes %>%
  dplyr::select( gene_name, locus_tag,  gene_start = start,  gene_end = end, strand) %>%
  dplyr::filter(!is.na(locus_tag)) %>%
  dplyr::distinct()


# Keep genes with one CDS record
gff_genes <- gff_genes %>%
  dplyr::group_by(locus_tag) %>%
  dplyr::filter(dplyr::n() == 1) %>%
  dplyr::ungroup()

# Create gene matching table, matching gene names
gene_matches <- gff_genes %>%
  dplyr::select(gene_id = gene_name, locus_tag, gene_start, gene_end, strand)

# Match using Rv locus tags
locus_matches <- gff_genes %>%
  dplyr::select(gene_id = locus_tag, locus_tag, gene_start, gene_end, strand)

# Combine gene names and locus tags
matching_table <- dplyr::bind_rows(gene_matches,  locus_matches) %>%
  dplyr::filter(!is.na(gene_id),  gene_id != "") %>%
  dplyr::distinct()

# Check unclear gene names
unsure_genes <- matching_table %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise(
    matching_genes = dplyr::n_distinct(locus_tag),
    .groups = "drop") %>%
  dplyr::filter(matching_genes > 1)


# Save unclear genes and remove unclear gene names
readr::write_csv(unsure_genes,  "unsure_genes.csv")

matching_table <- matching_table %>%
  dplyr::filter(!gene_id %in% unsure_genes$gene_id)

# Save gene matching table
readr::write_csv(matching_table,  "gene_matching_table.csv")

# Match TnSeq genes to GFF
tnseq_matched <- tnseq_counts %>%
  dplyr::left_join(matching_table,  by = "gene_id")

# Save unmatched TnSeq genes
unmatched_tnseq <- tnseq_matched %>%
  dplyr::filter(is.na(locus_tag)) %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct()

readr::write_csv(unmatched_tnseq,  "unmatched_tnseq.csv")

# Make TnSeq gene identifier standard
tnseq_final <- tnseq_matched %>%
  dplyr::filter(!is.na(locus_tag)) %>%
  dplyr::mutate(gene_id = locus_tag) %>%
  dplyr::select(gene_id,  position,  insertions) %>%
  dplyr::group_by(gene_id,  position) %>%
  dplyr::summarise(insertions = max(insertions),  .groups = "drop")

# Change UniProt gene aliases into rows
dom_aliases <- dom_df %>%
  dplyr::mutate(gene_names = tolower(gene_names)) %>%
  tidyr::separate_rows(gene_names,  sep = "\\s+") %>%
  dplyr::rename(gene_id = gene_names) %>%
  dplyr::filter(!is.na(gene_id),  gene_id != "")

# Match UniProt genes to GFF
dom_matched <- dom_aliases %>%
  dplyr::left_join(matching_table,  by = "gene_id")

# Check if one UniProt protein matches multiple genes
unsure_proteins <- dom_matched %>%
  dplyr::filter(!is.na(locus_tag)) %>%
  dplyr::group_by(Entry) %>%
  dplyr::summarise(matching_genes = dplyr::n_distinct(locus_tag),
    .groups = "drop") %>%
  dplyr::filter(matching_genes > 1)

# Save unclear proteins
readr::write_csv(unsure_proteins, "unsure_proteins.csv")

# Remove unclear proteins
dom_matched <- dom_matched %>%
  dplyr::filter(!Entry %in% unsure_proteins$Entry)

# Save unmatched UniProt proteins
unmatched_domains <- dom_matched %>%
  dplyr::group_by(Entry) %>%
  dplyr::filter(all(is.na(locus_tag))) %>%
  dplyr::ungroup() %>%
  dplyr::select(Entry,  gene_id) %>%
  dplyr::distinct()

readr::write_csv(unmatched_domains,  "unmatched_domains.csv")


# Keep matched UniProt proteins
dom_matched <- dom_matched %>%
  dplyr::filter(!is.na(locus_tag)) %>%
  dplyr::distinct(Entry,  locus_tag,  .keep_all = TRUE)


# Check protein and CDS lengths
dom_matched <- dom_matched %>%
  dplyr::mutate(cds_length = gene_end - gene_start + 1,
    length_difference = abs(cds_length - protein_length * 3)
  )

# Save proteins with different lengths
length_mismatch <- dom_matched %>%
  dplyr::filter(is.na(length_difference) |
      length_difference > 3)

readr::write_csv(length_mismatch,  "length_mismatch.csv")

# Keep proteins with compatible lengths
dom_matched <- dom_matched %>%
  dplyr::filter(!is.na(length_difference),  length_difference <= 3)

# Extract individual UniProt domains
dom_matched <- dom_matched %>%
  dplyr::mutate(domain_info = stringr::str_extract_all(
      domain_info,  'DOMAIN\\s+\\d+\\.\\.\\d+;\\s*/note="[^"]*"')
  ) %>%
  tidyr::unnest(domain_info)


# Extract amino acid coordinates
dom_matched <- dom_matched %>%
  dplyr::mutate(
    aa_start = as.numeric(stringr::str_extract(domain_info,"(?<=DOMAIN )\\d+"
      )),
    aa_end = as.numeric(stringr::str_extract(  domain_info,"(?<=\\.\\.)\\d+"
      )),
    domain_name = stringr::str_match(domain_info,'/note="([^"]+)"')[, 2]
  )


# Check amino acid coordinates
dom_matched <- dom_matched %>%
  dplyr::filter(!is.na(aa_start),  !is.na(aa_end),
    aa_start >= 1,
    aa_start <= aa_end,
    aa_end <= protein_length
  )

# Convert amino acid to genomic coordinates
dom_matched <- dom_matched %>%
  dplyr::mutate(domain_start = dplyr::case_when(
      strand == "+" ~
        gene_start + (aa_start - 1) * 3,
      strand == "-" ~
        gene_end - aa_end * 3 + 1,
      TRUE ~ NA_real_
    ),
    
    domain_end = dplyr::case_when(
      strand == "+" ~
        gene_start + aa_end * 3 - 1,
      strand == "-" ~
        gene_end - (aa_start - 1) * 3,
      TRUE ~ NA_real_
    )
  )


# Check genomic coordinates
dom_matched <- dom_matched %>%
  dplyr::filter(!is.na(domain_start),  !is.na(domain_end),
    domain_start >= gene_start,
    domain_end <= gene_end,
    domain_start <= domain_end)



# Save final domain data
domains_final <- dom_matched %>%
  dplyr::mutate(gene_id = locus_tag) %>%
  dplyr::select(gene_id,  domain_name,  domain_start,  domain_end) %>%
  dplyr::distinct()

# Save TnSeq data
readr::write_csv(tnseq_final,  "tnseq_clean.csv")

# Save domain data
readr::write_csv(domains_final,  "domains_final.csv")

# Check if data overlaps
tnseq <- readr::read_csv("tnseq_clean.csv",  show_col_types = FALSE)

domains <- readr::read_csv("domains_final.csv",  show_col_types = FALSE)

# Find common genes
common_genes <- intersect(na.omit(unique(tnseq$gene_id)),  na.omit(unique(domains$gene_id)))

# Display results
cat("Total TnSeq genes:",  length(unique(tnseq$gene_id)),  "\n")

cat("Total domain genes:",  length(unique(domains$gene_id)),  "\n")

cat("Overlapping genes:",  length(common_genes),  "\n")

cat("Total genomic domains:",  nrow(domains),  "\n")

# Example overlapping genes
head(common_genes, 20)