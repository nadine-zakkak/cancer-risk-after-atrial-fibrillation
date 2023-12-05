library(dplyr)
path_lookup <- "./look_ups/"

cancer_codelist <- read.csv(paste0(path_lookup, "cancer_classification_new_final_v1.csv")) |>
  mutate(icd10_3dig = trimws(icd10_3dig),
         icd10_4dig = trimws(icd10_4dig))
body_region <- read.delim(paste0(path_lookup, "body_region.txt"))
organ_system <- read.delim(paste0(path_lookup, "organ_system_v1.txt"))

rm(path_lookup)
