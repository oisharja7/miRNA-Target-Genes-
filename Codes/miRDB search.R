setwd("D:/DEG analysis_oisharja/Significant csvs")

# read mirna deg file
mirna_deg <- read.csv("sig_prim_norm_miRNA.csv", stringsAsFactors = FALSE)

# extract miRNA names 
mirna_list <- unique(mirna_deg$SystematicName)

# View the list of miRNAs
print(mirna_list)

library(tidyverse)

# rank top 50
top_50 <- mirna_deg %>%
  mutate(RankScore = -log10(adj.P.Val) * abs(logFC)) %>% # calculating composite score
  arrange(desc(RankScore)) %>% # sort in descending order
  head(50)

# extract top miRNA names 
mirna_top_list <- unique(top_50$SystematicName)

# View the list of miRNAs
print(mirna_top_list) 

# export 
write.csv(top_50, "top_50_miRNA.csv")

# **removing duplicate probes**



# Group by miRNA and select the row with the highest absolute logFC
unique_mirnas <- mirna_deg %>%
  group_by(SystematicName) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1) %>%
  ungroup()

# Now sort the unique miRNAs by the composite ranking score in descending order
top50_unique <- unique_mirnas %>%
  arrange(desc(RankScore)) %>%
  head(50)

# export 
write.csv(unique_mirnas, "unique_mirnas.csv")

print(unique_mirnas$SystematicName)



# Second contrast (Rec - Norm)

# read 
mirna_rec_norm <- read.csv("sig_rec_norm_miRNA.csv", stringsAsFactors = FALSE)

# calculate composote score 
mirna_rec_norm <- mirna_rec_norm %>%
  mutate(RankScore = -log10(adj.P.Val) * abs(logFC)) 


# Group by miRNA and select the row with the highest absolute logFC
unique_mirnas_rec_norm <- mirna_rec_norm %>%
  group_by(SystematicName) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1) %>%
  ungroup()

top50_rec_norm <- unique_mirnas_rec_norm %>%
  arrange(desc(RankScore)) %>%
  head(50)

# export
write.csv(top50_rec_norm, "top50_rec_norm.csv")

print(top50_rec_norm$SystematicName)

head(top50_rec_norm, 50)


# clean annotated gene list 

# read 
deg_genes <- read.csv("Annotated_Significant_Genes.csv", stringsAsFactors = FALSE)

# clean N/A values 
deg_genes_clean <- deg_genes[!is.na(deg_genes$GENE_SYMBOL), ]

# Extract a vector of unique genes
deg_gene_list <- unique(deg_genes_clean$GENE_SYMBOL)
print(deg_gene_list)


# load mirna target genes 
mirna_target <- read.csv("Combine Sheets result_ miRDB target prediction data.csv", stringsAsFactors = FALSE)

# Extract a vector of unique target gene symbols 
target_gene_list <- unique(mirna_target$Gene.Symbol)

print(target_gene_list)

# overlap between miRNA targets and deg
common_genes <- intersect(target_gene_list, deg_gene_list)

print(common_genes)

# convert as dataframe 
common_genes_df <- data.frame(Gene = common_genes)
head(common_genes_df)

# export
write.csv(common_genes_df, "Common_Genes.csv", row.names = FALSE)


## load new target batch ##
mirna_target_rec_norm <- read.csv("Combine Sheets result_rec_norm  - Combined data.csv", stringsAsFactors = FALSE)

# Extract a vector of unique target gene symbols 
mirna_target_list_rec_norm <- unique(mirna_target_rec_norm$Gene.Symbol)

print(mirna_target_list_rec_norm)


# overlap between miRNA targets and deg
common_genes_rec_norm <- intersect(mirna_target_list_rec_norm, deg_gene_list)

print(common_genes_rec_norm)


# convert as dataframe 
common_genes_rec_norm_df <- data.frame(Gene = common_genes_rec_norm)
head(common_genes_rec_norm_df)

# export
write.csv(common_genes_rec_norm_df, "common_genes_rec_norm_df.csv", row.names = FALSE)
