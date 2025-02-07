# load libraries

library(here)
library(readxl)
library(phyloseq)
library(janitor)
library(tidyverse)
library(Wrench)

#------------------------------ Load Datasets and initial cleaning--------------------

metadata <- read_xlsx(here("data/data-copy_4-11/For Ellen/metadata.xlsx"), skip = 1,
                      col_types = c(rep("text", 9), "numeric")) %>% 
  clean_names() %>% 
  rename(dec_isolate_id = dec_isolate_id_s,
         sequenced_16s = x16s_r_rna_amplicon_sequenced_from_whole_stool) %>% 
  mutate(alternate_metagenome_id = ifelse(alternate_metagenome_id == "NA", NA, alternate_metagenome_id)) %>% 
  mutate(dec_isolate_id = ifelse(dec_isolate_id == "NA", NA, dec_isolate_id)) %>% 
  mutate(sequenced_16s = sequenced_16s == "yes",
         dec_isolate_whole_genome_sequenced = dec_isolate_whole_genome_sequenced == "yes",
         shotgun_metagenome_sequenced_from_whole_stool = shotgun_metagenome_sequenced_from_whole_stool == "yes") %>% 
  mutate(diarrhea = diarrhea == "Case") # make 1/0, 1=Case, 0 = control

write_csv(metadata, here("data/formatted_data/cleaned_metadata.csv"))



alpha_diversity_16s <- read_csv(here("data/data-copy_4-11/For Ellen/16S rRNA amplicon data/alpha diversity/16S_alpha_diversity.csv")) %>% 
  select(-`...1`) %>% 
  clean_names()

phyloseq_obj_16s <- read_rds(here("data/data-copy_4-11/For Ellen/16S rRNA amplicon data/taxonomy/phyloseq_obj_16S.rds"))

metagenome_func_genes_annotated_functions <- read_csv(here("data/data-copy_4-11/For Ellen/Metagenome data/functional genes/metagenome_func_genes_annotated_functions.csv"))  %>% 
  rename(ko = `...1` )
pathway_annot <- read_csv(here("data/data-copy_4-11/For Ellen/Metagenome data/functional genes/KEGG mapping/pathway_annot.csv")) %>% 
  rename(pathway = `...1`)

path_list <- read_tsv(here("data/data-copy_4-11/For Ellen/Metagenome data/functional genes/KEGG mapping/path_list.txt"))
path_ko <- read_tsv(here("data/data-copy_4-11/For Ellen/Metagenome data/functional genes/KEGG mapping/path_ko.txt")) %>% 
  filter(!str_detect(path, "ko"))

VFDB_tab_trans <- read_csv(here("data/data-copy_4-11/For Ellen/Metagenome data/virulence factors/VFDB_tab_trans.csv"))

phyloseq_obj_metagenome <- read_rds(here("data/data-copy_4-11/For Ellen/Metagenome data/taxonomy/phyloseq_obj_metagenome.rds"))

nonpareil_alpha_diversity <- read_csv("data/data-copy_4-11/For Ellen/Metagenome data/alpha diversity/nonpareil_alpha_diversity.csv")

#------------------------------- Format datasets for joining ---------------------

# 16s

# we will use the estimated diversity 
alpha_diversity_16s_formatted <- alpha_diversity_16s %>% 
  rename(sample_id = sample) %>% 
  select(estimator, sample_id, diversity) %>% 
  pivot_wider(names_from = diversity, 
              values_from = estimator) %>% 
  clean_names()

# 16S Phyloseq Taxonomy Data

# find minumum count
phyloseq_16s_counts_vector <- as.vector(otu_table(phyloseq_obj_16s))
phyloseq_16s_counts_vector[phyloseq_16s_counts_vector == 0] <- NA
min(phyloseq_16s_counts_vector, na.rm = TRUE)

phyloseq_16s_formatted <- otu_table(phyloseq_obj_16s) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename_all(.funs = ~str_c("phyloseq_taxa_", .x)) %>% 
  mutate(sample_id = rownames(.)) %>% 
  mutate(across(starts_with("phyloseq_taxa"), ~ifelse(.x == 0, 2, .x))) %>%  # replace 0's with minimum count (2)
  rowwise() %>% 
  mutate(geometric_mean = exp(mean(log(c_across(starts_with("phyloseq_taxa")))))) %>% # calculate per-subject geometric mean
  ungroup() %>% 
  mutate(across(starts_with("phyloseq_taxa"), ~.x/geometric_mean)) # perform the Centered Log Ratio Transform

# DEC Presence data
dec_pathotypes <- metadata %>% 
  filter(dec_isolate_whole_genome_sequenced)  %>% 
  pull(dec_pathotype) %>% 
  str_split(pattern = "( and )|(, )") %>% 
  unlist() %>% 
  str_remove("(and )") %>% 
  unique()

dec_pathotype_presence <- metadata %>% 
  filter(dec_isolate_whole_genome_sequenced, sequenced_16s) %>%
  select(sample_id, dec_pathotype) %>% 
  mutate(dec_pathotype_DAEC = str_detect(dec_pathotype, "DAEC"),
         dec_pathotype_ETEC = str_detect(dec_pathotype, "ETEC"),
         dec_pathotype_EAEC = str_detect(dec_pathotype, "EAEC"),
         dec_pathotype_EIEC = str_detect(dec_pathotype, "EIEC"),
         dec_pathotype_aEPEC = str_detect(dec_pathotype, "aEPEC"),
         dec_pathotype_tEPEC = str_detect(dec_pathotype, "tEPEC")) %>% 
  select(-dec_pathotype)

# metagenome

# functional genes
ko_to_pathway_mapping <- left_join(path_ko, path_list, by = "path") %>% 
  right_join(pathway_annot, by = "pathway")


metagenome_func_genes_annotated_functions_reduced <- metagenome_func_genes_annotated_functions %>% 
  filter(ko %in% ko_to_pathway_mapping$ko)

metagenome_func_genes_annotated_functions_reduced_flipped <- select(metagenome_func_genes_annotated_functions_reduced, -ko) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(sample_id = str_remove(rownames(.), "_prodigal_MG_contigs_blastn")) %>% 
  relocate(sample_id)  %>% 
  mutate(sample_id = str_remove(sample_id, "_"))
  
names(metagenome_func_genes_annotated_functions_reduced_flipped) <- c("sample_id", metagenome_func_genes_annotated_functions_reduced$ko)

# virulence genes

VFDB_tab_trans_flipped <- select(VFDB_tab_trans, -`...1`) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(sample_id = str_extract(rownames(.), "[A-Z]+_?[0-9]*")) %>% 
  relocate(sample_id)  %>% 
  mutate(sample_id = str_remove(sample_id, "_"))

names(metagenome_func_genes_annotated_functions_reduced_flipped) <- c("sample_id", metagenome_func_genes_annotated_functions_reduced$ko)


# taxa counts

phyloseq_metagenome_counts_vector <- as.vector(otu_table(phyloseq_obj_metagenome))
phyloseq_metagenome_counts_vector[phyloseq_metagenome_counts_vector == 0] <- NA
min_phyloseq_metagenome_counts_vector <- min(phyloseq_metagenome_counts_vector, na.rm = TRUE)
phyloseq_metagenome_formatted <- otu_table(phyloseq_obj_metagenome) %>% 
  t() %>% 
  as.data.frame() %>% 
  rename_all(.funs = ~str_c("phyloseq_taxa_", .x)) %>% 
  mutate(sample_id = rownames(.)) %>% 
  mutate(across(starts_with("phyloseq_taxa"), ~ifelse(.x == 0, min_phyloseq_metagenome_counts_vector, .x))) %>%  # replace 0's with minimum count (2)
  rowwise() %>% 
  mutate(geometric_mean = exp(mean(log(c_across(starts_with("phyloseq_taxa")))))) %>% # calculate per-subject geometric mean
  ungroup() %>% 
  mutate(across(starts_with("phyloseq_taxa"), ~.x/geometric_mean)) %>% # perform the Centered Log Ratio Transform
  select(-geometric_mean) %>% 
  mutate(sample_id = str_remove(sample_id, "_"))

# alpha diversity

nonpareil_alpha_diversity_clean <- nonpareil_alpha_diversity %>% 
  mutate(sample_id = str_remove(str_extract(str_extract(`...1`, "//.*\\.npo"), "[A-Z]+_?[0-9]*"), "_")) %>% 
  select(sample_id, diversity)


# ----------------------------- Build full dataset ---------------------------

# 16s
combined_data_16s <- metadata %>% 
  left_join(alpha_diversity_16s_formatted, by = "sample_id") %>% 
  left_join(phyloseq_16s_formatted, by = "sample_id") %>% 
  filter(sequenced_16s == TRUE)

write_csv(combined_data_16s, here("data/formatted_data/combined_data_16s.csv"))

# 16s with DEC pathotype

combined_data_16s_dec <- metadata %>% 
  left_join(alpha_diversity_16s_formatted, by = "sample_id") %>% 
  left_join(phyloseq_16s_formatted, by = "sample_id") %>% 
  left_join(dec_pathotype_presence, by = "sample_id") %>% 
  filter(sequenced_16s == TRUE, dec_isolate_whole_genome_sequenced)

write_csv(combined_data_16s_dec, here("data/formatted_data/combined_data_16s_dec.csv"))


# Bacterial Families to Group Indices by

x_16s <- select(combined_data_16s, -diarrhea, -sample_id,
                -c(alternate_metagenome_id, 
                   sequenced_16s, 
                   dec_isolate_whole_genome_sequenced, 
                   dec_isolate_id, 
                   shotgun_metagenome_sequenced_from_whole_stool, 
                   dec_infection, 
                   dec_pathotype, 
                   participant_age_years, 
                   geometric_mean))

taxa_table <- tax_table(phyloseq_obj_16s) %>% 
  as.data.frame() %>% 
  mutate(taxa_id = str_c("phyloseq_taxa_", rownames(.)))

variable_positions <- data.frame(variable_name = names(x_16s), position = 1:length(names(x_16s)))

family_variables <- taxa_table %>%
  left_join(variable_positions, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Family)) %>% 
  group_by(Family) %>% 
  filter(n() > 0) %>% 
  group_split() %>% 
  map(~pull(.x, position))

names(family_variables) <- taxa_table %>%
  left_join(variable_positions, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Family)) %>% 
  group_by(Family) %>% 
  filter(n() > 0) %>% 
  group_keys() %>% 
  pull(Family) # previously used unique, this messed with order

phylum_variables <- taxa_table %>%
  left_join(variable_positions, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Phylum)) %>% 
  group_by(Phylum) %>% 
  filter(n() > 0) %>% 
  group_split() %>% 
  map(~pull(.x, position))

names(phylum_variables) <- taxa_table %>%
  left_join(variable_positions, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Phylum)) %>% 
  group_by(Phylum) %>% 
  filter(n() > 0) %>% 
  group_keys() %>% 
  pull(Phylum) # previously used unique, this messed with order

taxa_variables <- list("taxa_variables" = variable_positions %>% 
  filter(str_detect(variable_name, "phyloseq_taxa")) %>% 
  pull(position))

diversity_variables <- list("diversity_variables" = c(1, 2, 3))

variable_sets_16s <- c(family_variables, 
                       phylum_variables, 
                       taxa_variables, 
                       diversity_variables)

write_rds(variable_sets_16s, "data/variable_sets_16s.rds")

# add the dec pathotype sets for the dec table

x_16s_dec <- select(combined_data_16s_dec, -diarrhea, -sample_id,
                -c(alternate_metagenome_id, 
                   sequenced_16s, 
                   dec_isolate_whole_genome_sequenced, 
                   dec_isolate_id, 
                   shotgun_metagenome_sequenced_from_whole_stool, 
                   dec_infection, 
                   dec_pathotype, 
                   participant_age_years, 
                   geometric_mean))

variable_positions_dec <- data.frame(variable_name = names(x_16s_dec), position = 1:length(names(x_16s_dec)))

dec_pathotype_variables <- list("dec_pathotype_variables" = filter(variable_positions_dec, str_detect(variable_name, "dec_pathotype_")) %>% 
  pull(position))
variable_sets_16s_dec <- c(family_variables, 
                       phylum_variables, 
                       taxa_variables, 
                       diversity_variables,
                       dec_pathotype_variables)

write_rds(variable_sets_16s_dec, "data/variable_sets_16s_dec.rds")


# metagenome
combined_data_metagenome <- metadata %>% 
  mutate(alternate_metagenome_id = ifelse(alternate_metagenome_id == "R54", NA, alternate_metagenome_id),
         alternate_metagenome_id = ifelse(sample_id == "R055", "R55", alternate_metagenome_id)) %>% # fix data error
  filter(!is.na(alternate_metagenome_id)) %>% 
  mutate(sample_id_phyloseq = ifelse(str_detect(sample_id, "[A-Z]0+[1-9]"), str_remove(sample_id, "0+"), sample_id)) %>%
  left_join(nonpareil_alpha_diversity_clean, by = c("alternate_metagenome_id" = "sample_id")) %>% 
  left_join(phyloseq_metagenome_formatted, by = c("sample_id_phyloseq" = "sample_id")) %>%  
  left_join(metagenome_func_genes_annotated_functions_reduced_flipped, by = c("alternate_metagenome_id" = "sample_id")) %>% 
  left_join(VFDB_tab_trans_flipped, by = c("alternate_metagenome_id" = "sample_id")) %>% 
  select(-c(sequenced_16s, 
            dec_isolate_whole_genome_sequenced, 
            dec_isolate_id, 
            shotgun_metagenome_sequenced_from_whole_stool, 
            dec_pathotype, 
            participant_age_years, 
            sample_id_phyloseq))



write_csv(combined_data_metagenome, here("data/formatted_data/combined_data_metagenome.csv"))

x_metagenome <- select(combined_data_metagenome, -diarrhea, -sample_id, -alternate_metagenome_id, -dec_infection)
variable_positions_metagenome <- data.frame(variable_name = names(x_metagenome), position = 1:length(names(x_metagenome)))

# creating variable sets to run VIMP on
taxa_table_metagenome <- tax_table(phyloseq_obj_metagenome) %>% 
  as.data.frame() %>% 
  mutate(taxa_id = str_c("phyloseq_taxa_", rownames(.)))

taxa_phylums_metagenome <- taxa_table_metagenome %>% 
  left_join(variable_positions_metagenome, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Phylum)) %>% 
  group_by(Phylum) %>% 
  filter(n() > 0) %>% 
  group_split() %>% 
  map(~pull(.x, position))

names(taxa_phylums_metagenome) <- taxa_table_metagenome %>% 
  left_join(variable_positions_metagenome, by = c("taxa_id" = "variable_name")) %>% 
  filter(!is.na(Phylum)) %>% 
  group_by(Phylum) %>% 
  filter(n() > 0) %>% 
  group_keys() %>% 
  pull(Phylum)

functional_gene_classes_metagenome <- variable_positions_metagenome %>% 
  inner_join(ko_to_pathway_mapping, by = c("variable_name" = "ko")) %>%
  group_by(Secondary) %>% 
  filter(n() > 0) %>% 
  group_split() %>% 
  map(~pull(.x, position)) %>% 
  map(unique)

names(functional_gene_classes_metagenome) <- variable_positions_metagenome %>% 
  inner_join(ko_to_pathway_mapping, by = c("variable_name" = "ko")) %>% 
  group_by(Secondary) %>% 
  group_keys(Secondary) %>% 
  pull(Secondary)

virulence_genes_metagenome <- list("virulence_genes" = variable_positions_metagenome %>% filter(str_detect(variable_name, "^V")) %>% pull(position))
functional_genes_metagenome <- list("functional_genes" = variable_positions_metagenome %>% filter(str_detect(variable_name, "^K")) %>% pull(position))
taxa_counts_metagenome <- list("taxa_counts" = variable_positions_metagenome %>% filter(str_detect(variable_name, "^phyloseq")) %>% pull(position))
diversity_metagenome <- list("diversity" = variable_positions_metagenome %>% filter(str_detect(variable_name, "^diversity")) %>% pull(position))


variable_sets_metagenome <- c(taxa_phylums_metagenome, 
                              functional_gene_classes_metagenome, 
                              virulence_genes_metagenome, 
                              functional_genes_metagenome, 
                              taxa_counts_metagenome,
                              diversity_metagenome)

write_rds(variable_sets_metagenome, "data/variable_sets_metagenome.rds")
