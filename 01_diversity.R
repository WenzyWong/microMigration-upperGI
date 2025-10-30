setwd("/data/yzwang/project/migration/")

DIR_RES <- "/data/yzwang/project/migration/result/"

library(dplyr)
library(tibble)
library(paletteer)

feature.g <- read.delim("./qiime2_analysis/merged_results/feature_table_genus.tsv", 
                        comment.char="#") %>%
  filter(grepl("g_", OTU.ID))

mtx.g <- feature.g %>%
  mutate(OTU.ID = gsub("..*g__", "", OTU.ID)) %>%
  filter(OTU.ID != "uncultured")

row.names(mtx.g) <- mtx.g$OTU.ID
mtx.g <- as.matrix(mtx.g[ , 2:ncol(mtx.g)])

abund.g <- apply(mtx.g, MARGIN = 2, FUN = function(x) x/sum(x) * 100)
avg_abund.g <- apply(abund.g, MARGIN = 1, FUN = mean) %>% sort(., decreasing = T)

meta.barrets <- read.csv("./barretts_with_oral/SraRunTable.csv")
table(meta.barrets$env_medium)
table(meta.barrets$pathology)

meta.escc <- read.csv("./escc_with_control/SraRunTable.csv")
table(meta.escc$Host_disease)
table(meta.escc$host_tissue_sampled)
# G: capsule-sponge
# W: oral cavity swabs
# TK: gastric biopsy

# Classification
meta.combined <- bind_rows(
  meta.barrets %>%
    select(Run, pathology, env_medium) %>%
    mutate(
      disease_status = case_when(
        pathology == "Adenocarcinoma" ~ "tumour",
        pathology == "Control" ~ "normal",
        pathology %in% c("High grade dysplasia", "Indefinite for dysplasia",
                         "Low grade dysplasia") ~ "dysplasia",
        pathology == "Intestinal metaplasia, no dysplasia" ~ "metaplasia",
        TRUE ~ NA_character_
      ),
      tissue_type = case_when(
        env_medium %in% c("Oral wash", "Saliva") ~ "oral",
        env_medium == "Squamous" ~ "squamous",
        env_medium == "Cardia and Barretts" ~ "junction",
        TRUE ~ NA_character_
      )
    ),
  meta.escc %>%
    select(Run, Host_disease, host_tissue_sampled) %>%
    mutate(
      disease_status = case_when(
        Host_disease == "neoplasia" ~ "tumour",
        Host_disease == "control" ~ "normal",
        Host_disease == "high-risk" ~ "highrisk",
        TRUE ~ NA_character_
      ),
      tissue_type = case_when(
        host_tissue_sampled == "W" ~ "oral",
        host_tissue_sampled == "G" ~ "whole_esophagus",
        host_tissue_sampled == "TK" ~ "gastric",
        TRUE ~ NA_character_
      )
    )
) %>%
  select(Run, disease_status, tissue_type) %>%
  filter(!is.na(disease_status), !is.na(tissue_type))

# Create the disease groups list
disease.groups <- meta.combined %>%
  mutate(group_name = paste(disease_status, tissue_type, sep = ".")) %>%
  group_by(group_name) %>%
  summarise(runs = list(Run), .groups = "drop") %>%
  deframe()

source("utils/plot_beta_diversity.R")
# Check batch effects
abund.g1 <- abund.g[,colnames(abund.g) %in% meta.barrets$Run]
abund.g2 <- abund.g[,colnames(abund.g) %in% meta.escc$Run]
beta.batch <- plot_beta_diversity(
  abundance_list = list(Cohort1 = abund.g1, Cohort2 = abund.g2),
  colors = c("#485682", "#5C8447"),
  output_file = file.path(DIR_RES, "Beta_diversity_cohorts.pdf")
)

# Extract abundance matrices for each disease group
abund.groups <- lapply(disease.groups, function(runs) {
  cols <- colnames(abund.g)[colnames(abund.g) %in% runs]
  abund.g[, cols, drop = FALSE]
})
sapply(abund.groups, ncol)

beta.nc <- plot_beta_diversity(
  abundance_list = abund.groups[grepl("normal|tumour", names(abund.groups))],
  color = paletteer_d("ggsci::category20_d3")[1:10],
  output_file = file.path(DIR_RES, "Beta_diversity_nc.pdf")
)
abund.nc <- abund.groups[grepl("normal|tumour", names(abund.groups))]

beta.oral <- plot_beta_diversity(
  abundance_list = abund.nc[grepl("oral", names(abund.nc))],
  color = paletteer_d("ggsci::flattastic_flatui")[c(4, 1)],
  output_file = file.path(DIR_RES, "Beta_diversity_oral.pdf")
)
beta.esph <- plot_beta_diversity(
  abundance_list = abund.nc[grepl("whole_esophagus", names(abund.nc))],
  color = paletteer_d("ggsci::flattastic_flatui")[c(5, 2)],
  output_file = file.path(DIR_RES, "Beta_diversity_esophagus.pdf")
)
beta.gas <- plot_beta_diversity(
  abundance_list = abund.nc[grepl("gastric", names(abund.nc))],
  color = paletteer_d("ggsci::flattastic_flatui")[c(6, 3)],
  output_file = file.path(DIR_RES, "Beta_diversity_gastric.pdf")
)

beta.tumour <- plot_beta_diversity(
  abundance_list = abund.groups[grepl("tumour", names(abund.groups))],
  color = paletteer_d("ggsci::flattastic_flatui")[1:5],
  output_file = file.path(DIR_RES, "Beta_diversity_tumour.pdf")
)
beta.normal <- plot_beta_diversity(
  abundance_list = abund.groups[grepl("normal", names(abund.groups))],
  color = paletteer_d("ggsci::flattastic_flatui")[1:5],
  output_file = file.path(DIR_RES, "Beta_diversity_normal.pdf")
)

beta.highrisk <- plot_beta_diversity(
  abundance_list = abund.groups[4:6],
  color = paletteer_d("ggsci::flattastic_flatui")[4:6],
  output_file = file.path(DIR_RES, "Beta_diversity_highrisk.pdf")
)
beta.mildrisk <- plot_beta_diversity(
  abundance_list = abund.groups[7:9],
  color = paletteer_d("ggsci::flattastic_flatui")[7:9],
  output_file = file.path(DIR_RES, "Beta_diversity_mildrisk.pdf")
)
beta.normal <- plot_beta_diversity(
  abundance_list = abund.groups[10:12],
  color = paletteer_d("ggsci::flattastic_flatui")[10:12],
  output_file = file.path(DIR_RES, "Beta_diversity_normal.pdf")
)
