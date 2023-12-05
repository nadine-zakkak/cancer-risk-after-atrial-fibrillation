library(dplyr)

output_dir <- "./appendix_codelists/"

# Cancer history ICD-10 & ICD-9 codes ----
#icd-10 codes
cancer_pheno <- read.delim("./phenotypes/cancersite_v1.2_codelist_4digicd10.txt")
cancer_pheno_hx <- cancer_pheno |> 
  filter(cancer_site_number != 0) |>
  select(icd10_3dig, icd10_4dig, cancer_site_desc, cancer_group_desc) |> as_tibble()

# MERGE ICD-10 4 DIGITS INTO ICD-10 3 DIGITS WHEN POSSIBLE

#icd-9 codes
cancer_icd9 <- c(setdiff(140:208, 173), "2252", "2254")

# Cancer outcomes ICD-10 codes ----
source("./scripts/00_read_codelists.R")
cancer_codelist_outcome <- cancer_codelist |> filter(outcome == 1) |> as_tibble()

cancer_codelist_outcome <- cancer_codelist_outcome |>
  select(icd10_3dig, icd10_4dig, `Body Region`=body_region_desc, `Organ System`=organ_system_desc) |>
  arrange(icd10_4dig)

write.csv(cancer_codelist_outcome, paste0(output_dir, "cancer_outcomes.csv"), row.names = F)

# Covariates ----
## Diabetes ----
lookup_dm_gprd <- read.csv("./look_ups/dm_cprd_codelist_2.csv", header=TRUE, sep=",")
lookup_dm_gprd <- lookup_dm_gprd |> filter(category %in% c(3,4, 6)) |> select(category, medcode, readterm , readcode) 
lookup_dm_hes <- read.csv("./look_ups/dm_hes_codelist_1.csv", header=TRUE, sep=",")
lookup_dm_hes <- lookup_dm_hes |> filter(category %in% c(3,4,6)) |> select(category, icd_code, icd_term) 
lookup_dm_cprd <- read.csv("./look_ups/diabdiag_cprd_codelist_1.csv", header=TRUE, sep=",")
lookup_dm_cprd <- lookup_dm_cprd |> filter(category %in% c(3,4)) |> select(category, medcode, readterm, readcode)

dm_readcodes <- lookup_dm_gprd |> select(readterm, readcode) |> 
  bind_rows(lookup_dm_cprd |> select(readterm, readcode)) |>
  distinct()
dm_icd10 <- lookup_dm_hes |> select(icd_term, icd_code) |> distinct()

write.csv(dm_readcodes, paste0(output_dir, "diabetes_readcodes.csv"), row.names = F)
write.csv(dm_icd10, paste0(output_dir, "diabetes_icd10.csv"), row.names = F)

rm(list=ls(pattern=".*dm.*"))

## Hypertension -------
lookup_ht_cprd <-  read.csv("./look_ups/ht_cprd_codelist_1.csv", header=TRUE, sep=",")
lookup_ht_cprd <- lookup_ht_cprd |> filter(category %in% c(3,4)) |> select(category, medcode, readterm , readcode)
lookup_ht_hes <- read.csv("./look_ups/ht_hes_codelist_1.csv", header=TRUE, sep=",")
lookup_ht_hes <- lookup_ht_hes |> filter(category %in% c(3,4)) |> select(category, icd_code, icd_term)

ht_readcodes <- lookup_ht_cprd |> select(readterm, readcode) |> distinct()
ht_icd10 <- lookup_ht_hes |> select(icd_term, icd_code) |> distinct()

write.csv(ht_readcodes, paste0(output_dir, "hypertension_readcodes.csv"), row.names = F)
write.csv(ht_icd10, paste0(output_dir, "hypertension_icd10.csv"), row.names = F)

rm(list=ls(pattern=".*ht.*"))

## Smoking ----
# assuming non-smokers if no record of smoking
#-- cprd: 1-non 2-ex 3-ever 4-current
#-- hes: 4-current
#-- additional: 1-yes, 2-no, 3-ex
lookup_smk_cprd <- read.csv("./look_ups/smoking_status_cprd_codelist_1.csv")
lookup_smk_cprd[13,"category"] <- 1 #changed "current non-smoker" to category 1 manually - check if keep that way or no
lookup_smk_hes <- read.csv("./look_ups/smoking_status_hes_codelist_2.csv")
lookup_smk_hes <- lookup_smk_hes |> filter(!is.na(category)) |> select(icd_code, icd_term, category)
regex_smk_pattern <- paste(paste0(lookup_smk_hes$icd_code, ".*"), collapse="|")
lookup_bnf <- read.delim("./look_ups/bnfcodes.txt", colClasses = c("double", "character"))
bnf_smk_regex <- "^041002.*" #nicotine dependence 

#record of smoking
lookup_smk_cprd_ever <- lookup_smk_cprd |> filter(category %in% c(2,3,4))
lookup_smk_hes_ever <- lookup_smk_hes |> filter(category %in% c(2,3,4))

smk_readcodes <- lookup_smk_cprd_ever |> select(readterm, readcode) |> distinct()
smk_icd10 <- lookup_smk_hes_ever |> select(icd_term, icd_code) |> distinct()
smk_bnf <- lookup_bnf |> filter(grepl(bnf_smk_regex, bnf)) |> select(bnf) |> pull()

write.csv(smk_readcodes, paste0(output_dir, "smk_readcodes.csv"), row.names = F)
write.csv(smk_icd10, paste0(output_dir, "smk_icd10.csv"), row.names = F)
write.csv(smk_bnf, paste0(output_dir, "smk_bnf.csv"), row.names = F)

rm(list=ls(pattern=".*smk.*"))
rm(lookup_bnf)

## Alcohol ----
# 1:non, 2: ex, 3:occ, 4:current, 5:excess, 6:binge
# 1,2,3: drinker/current; 4,5: excess
# additional: 1: yes, 2:no, 3:ex

# New categories:
# 0: Missing, 1: non-drinker, 2: Ex drinker, 3: Current drinker
lookup_alc_drink <- read.csv("./look_ups/alcohol_drinker_cprd_codelist_1.csv")
lookup_alc_unit <- read.csv("./look_ups/alcohol_units_cat_cprd_codelist_1.csv") |>
  mutate(readcode = as.character(readcode))

lookup_alc_dict <- data.frame(
  category_name = c("Missing", "Non", "Ex", "Current"),
  category_num = c(0,1,2,3))

lookup_alc_drink <-
  lookup_alc_drink |>
  mutate(category_new = case_when(
    category == 1 ~ 1,
    category == 2 ~ 2,
    category %in% c(3:6) ~ 3
  ))

lookup_alc_unit <- lookup_alc_unit |>
  mutate(category_new = case_when(
    category %in% c(1:5) ~ 3
  ))

alc_readcodes <- lookup_alc_drink |> select(readterm, readcode, category_new) |>
  bind_rows(lookup_alc_unit |> select(readterm, readcode, category_new)) |>
  filter(!is.na(category_new)) |>
  left_join(lookup_alc_dict, by = c("category_new"="category_num")) |>
  select(category_name, readterm, readcode, -category_new) |>
  rename("Category" = category_name,
         "Read term" = readterm,
         "Read code" = readcode)

alc_additional_ent <- data.frame(ent = as.character(c(1,2,3)), 
                             desc = c("yes", "no", "ex"),
                             category_new = c(3,1,2)) |>
  left_join(lookup_alc_dict, by = c("category_new"="category_num")) |>
  select(category_name, ent, -category_new, -desc) |>
  rename("Category" = category_name,
         "Data Entry"=ent)

write.csv(alc_readcodes, paste0(output_dir, "alc_readcodes.csv"), row.names = F)
write.csv(alc_additional_ent, paste0(output_dir, "alc_additional_ent.csv"), row.names = F)

rm(list=ls(pattern=".*alc.*"))
