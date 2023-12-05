rm(list=ls()); gc()

library(dplyr) #used
library(lubridate)
library(tidyr)

# Follow-up is saved by month

### Cases end date = min(death, outcome, 2015)
### Controls end date = min(death, AF, outcome, 2015)

# AF = 1 only if it is the earliest "outcome" to occur

# 15/03:
# which "outcome" takes priority 
## for now:
### 1) outcome of interest
### 2) af
### 3) death
### 4) end of study
progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 8_outcomes.R script"), file = progress_file, append = T)

path_lookup <- "./look_ups/"
path_rw <- "./data/final_data/"

final <- readRDS(paste0(path_rw, "confounders_final.rds"))
cancer <- readRDS(paste0(path_rw, "cancer_registration_final_cohort.rds"))  |> rename(epatid = e_patid)
death <- readRDS(paste0(path_rw, "death_patient_final_cohort.rds")) |> rename(epatid = e_patid)

cancer <- cancer |> mutate(site_icd10_o2 = trimws(site_icd10_o2))
death <- death |> mutate(cause = trimws(cause))

#codelists
source("./scripts/00_read_codelists.R")
cancer_codelist_outcome <- cancer_codelist |> filter(outcome == 1)

# Add "earliest end date"
## Cases end date = min(death, 2015)
## Controls end date = min(AF, death, 2015)
end_date <- as.Date('2015-12-31')
final <-
  final |>
  mutate(earliest_end = if_else(case,
                                pmin(death_date, end_date, na.rm = T),
                                pmin(af_date, death_date, end_date, na.rm = T)))

#9/March
# ------ death outcome ------
write(paste(Sys.time(), "Getting death outcome"), file = progress_file, append = T)

# flag if died due to cancer
death_cancer <-
  death |> select(epatid, dor, cause) |>
  mutate(cause_mod = gsub("\\.", "", cause),
         cause_mod = trimws(cause_mod)) |>
  mutate(cancer_died = cause_mod %in% cancer_codelist_outcome$icd10_4dig,
         died = T)

# flag eligible death (any cause and cancer)
death_cancer <-
  death_cancer |>
  left_join(final |> select(epatid, case, af, death_date, af_date, earliest_end), by = 'epatid') |>
  mutate(elig_cancerdied = cancer_died & death_date <= earliest_end,
         elig_died = died & death_date <= earliest_end)

table(death_cancer$cancer_died, death_cancer$elig_cancerdied, dnn = c("cancer death", "elig"), useNA = "ifany")
table(death_cancer$died, death_cancer$elig_died, dnn = c("death", "elig"), useNA = "ifany")
table(death_cancer$case, death_cancer$elig_cancerdied, dnn = c("case", "elig cancer death"), useNA = "ifany")
table(death_cancer$elig_died, death_cancer$elig_cancerdied, dnn = c("death", "cancer death"), useNA = "ifany")

#exploring died but not eligible
death_cancer |>
  filter(died, !elig_died) |>
  select(epatid, case, af, death_date, af_date, earliest_end, died, elig_died)

#exploring cancer died but not eligible
death_cancer |>
  filter(cancer_died, !elig_cancerdied) |>
  select(epatid, case, af, death_date, af_date, earliest_end, cancer_died, elig_cancerdied)

#join cancer death info with all patients info
final_death <-
  final |>
  left_join(death_cancer |> select(epatid, case, elig_cancerdied, elig_died), by = c('epatid', 'case')) |>
  mutate(elig_cancerdied = ifelse(is.na(elig_cancerdied), F, elig_cancerdied),
         elig_died = ifelse(is.na(elig_died), F, elig_died),
         af_outcome = !case & af & af_date <= earliest_end,
         outcome = case_when(
           elig_cancerdied ~ "cancer death",
           af_outcome ~ "af",
           elig_died ~ "other cause death", #if cancer death, already recorded as so
           TRUE ~ "end"
         ),
         follow = time_length(start %--% earliest_end, unit = "month")) |>
  rename(cancer_died_outcome = elig_cancerdied,
         died_outcome = elig_died,
         outcome_date = earliest_end) |>
  distinct()

table(final_death$case, final_death$cancer_died_outcome, dnn=c("case", "elig cancer death"), useNA = "ifany")
100*prop.table(table(final_death$case, final_death$cancer_died_outcome, dnn=c("case", "elig cancer death"), useNA = "ifany"), margin = 1)
100*prop.table(table(final_death$case, final_death$died_outcome, dnn=c("case", "death"), useNA = "ifany"), margin = 1)

saveRDS(final_death, paste0(path_rw, "cancer_death.RDS"))

#clean up
rm(death_cancer, final_death)
gc()

# ------ Cancer diagnosis -----
write(paste(Sys.time(), "Getting cancer incidence"), file = progress_file, append = T)

cancer_filter <-
  cancer |>
  mutate(cancer = site_icd10_o2 %in% cancer_codelist_outcome$icd10_4dig) 

cancer_filter_single <-
  cancer_filter |>
  left_join(final |> select(epatid, case, start, af, af_date, death_date, earliest_end), by = 'epatid') |>
  mutate(elig_cancer = cancer & diagnosisdatebest <= earliest_end) |>
  filter(elig_cancer) |> #keep only eligible cancers
  group_by(epatid, case, start, af, af_date, death_date, earliest_end, elig_cancer) |>
  summarise(diagnosisdatebest = as.Date(min(diagnosisdatebest, na.rm = T))) |> #get earliest eligible cancer date
  select(epatid, case, start, af, elig_cancer, cancer_date=diagnosisdatebest, af_date, death_date, earliest_end) |>
  distinct() |>
  ungroup()

table(cancer_filter_single$case, cancer_filter_single$elig_cancer, dnn = c("case", "elig cancer"), useNA = "ifany")

#join cancer death info with all patients info
final_cancer <-
  final |>
  left_join(cancer_filter_single |> select(epatid, case, cancer_date, elig_cancer), by = c('epatid', 'case')) |>
  mutate(elig_cancer = ifelse(is.na(elig_cancer), F, elig_cancer),
         elig_af = !case & af & af_date <= earliest_end & (is.na(cancer_date) | af_date <= cancer_date),
         outcome_date = if_else(elig_cancer,
                                pmin(cancer_date, earliest_end, na.rm = T),
                                earliest_end),
         outcome = case_when(
           cancer_date == outcome_date ~ "cancer",
           af_date == outcome_date & !case ~ "af",
           death_date == outcome_date ~ "death",
           end_date == outcome_date ~ "end"
         ),
         follow =  time_length(start %--% outcome_date, unit = "month")) |>
  select(-af) |>
  rename(cancer_outcome = elig_cancer,
         af_outcome = elig_af)

table(final_cancer$case, final_cancer$cancer_outcome, dnn = c("case", "elig cancer"), useNA = "ifany")
100*prop.table(table(final_cancer$case, final_cancer$cancer_outcome, dnn=c("case", "elig cancer"), useNA = "ifany"), margin = 1)
100*prop.table(table(final_cancer$case, final_cancer$outcome, dnn=c("case", "outcome"), useNA = "ifany"), margin = 1)

final_cancer |> group_by(epatid, case) |> summarise(n = n()) |> filter(n>1) #none

# organise columns in final dataset
final_cancer <- final_cancer |> 
  select(-alc_date)|>
  select(epatid, pracid, gender, age_start, case,
         diab, ht, smok, imd, alc,
         cancer_outcome, af_outcome, outcome,
         start, af_date, cancer_date, death_date, outcome_date, follow)

saveRDS(final_cancer, paste0(path_rw, "cancer_diag.RDS"))

#clean up
rm(cancer_filter, cancer_filter_single, final_cancer)
gc()

# ------ By sub categories -----
# only has patients with cancer
cancer_categ <-
  cancer |>
  inner_join(cancer_codelist_outcome |> select(icd10_4dig, body_region, body_region_desc,
                                               organ_system, organ_system_desc),
             by = c('site_icd10_o2'='icd10_4dig'))

## ------ By body region -----
write(paste(Sys.time(), "Getting cancer incidence - body region"), file = progress_file, append = T)

cancer_region <- cancer_categ |> 
  left_join(final |> select(epatid, case, earliest_end), by = 'epatid', multiple = "all") |>
  group_by(epatid, case, body_region) |> 
  filter(diagnosisdatebest == min(diagnosisdatebest, na.rm = T)) |> #keep the earliest diagnosis for each cancer category
  select(epatid, case, body_region, cancer_date=diagnosisdatebest, earliest_end) |>
  distinct() |>
  ungroup() |>
  mutate(elig_cancer = cancer_date <= earliest_end)  |>
  filter(elig_cancer) |> #keep only eligible cancers
  pivot_longer(cols = body_region, values_to = 'cancer_bodyregion') #long format

table(cancer_region$cancer_bodyregion, useNA = "ifany")
table(cancer_region$case, cancer_region$elig_cancer, dnn = c("case", "elig body region"), useNA = "ifany")
table(cancer_region$cancer_bodyregion, cancer_region$elig_cancer, dnn = c("body region", "elig body region"), useNA = "ifany")

# create "long" dataframe of patients repeated by the number of body regions there is
pat_long_region <- final |> 
  select(epatid, case) |> 
  slice(rep(1:n(), each = max(body_region$body_region_numb))) |>
  mutate(cancer_bodyregion = rep(1:max(body_region$body_region_numb),  times = nrow(final))) #add column for body region number

pat_long_region <- 
  pat_long_region |>
  left_join(cancer_region |> select(epatid, case, cancer_date, elig_cancer, cancer_bodyregion), by = c('epatid', 'case', 'cancer_bodyregion'))  |> #join with eligible cancer cases
  left_join(final, by = c('epatid', 'case'), multiple = "all") #join with patients info


#join cancer region info with all patients info
final_region <-
  pat_long_region |>
  mutate(elig_cancer = ifelse(is.na(elig_cancer), F, elig_cancer),
         af_outcome = !case & af & af_date <= earliest_end & (is.na(cancer_date) | af_date <= cancer_date),
         outcome_date = if_else(elig_cancer,
                                pmin(cancer_date, earliest_end, na.rm = T),
                                earliest_end),
         outcome = case_when(
           cancer_date == outcome_date ~ "cancer",
           af_date == outcome_date & !case ~ "af",
           death_date == outcome_date ~ "death",
           end_date == outcome_date ~ "end"
         ),
         outcome_detail=ifelse(outcome=="cancer",
                               cancer_bodyregion,
                               outcome),
         follow =  time_length(start %--% outcome_date, unit = "month")) |>
  select(-af) |>
  rename(cancer_outcome = elig_cancer)

final_region |>
  group_by(case, outcome_detail) |>
  summarise(n = n()) |>
  mutate(prop = 100*n/sum(n)) |>
  pivot_wider(id_cols = case, names_from = outcome_detail, values_from = prop)


setdiff(final$epatid, final_region$epatid) #all patients included
cancer_region |> filter(!is.na(cancer_bodyregion)) |> group_by(epatid, case) |> summarise(n = n()) |> filter(n>1)

#keep only required columns
final_region <-
  final_region |>
  select(epatid, pracid, gender, age_start, case,
         diab, ht, smok, imd, alc,
         af_outcome, cancer_outcome, cancer_bodyregion, outcome, outcome_detail,
         start, af_date, cancer_date, death_date, outcome_date, follow)

saveRDS(final_region, paste0(path_rw, "cancer_region.RDS"))

# ------ By organ system -----
write(paste(Sys.time(), "Getting cancer incidence - organ system"), file = progress_file, append = T)

cancer_sys <- cancer_categ |> 
  left_join(final |> select(epatid, case, earliest_end), by = 'epatid', multiple = "all") |>
  group_by(epatid, case, organ_system) |> 
  filter(diagnosisdatebest == min(diagnosisdatebest, na.rm = T)) |> #keep the earliest diagnosis for each cancer category
  select(epatid, case, organ_system, cancer_date=diagnosisdatebest, earliest_end) |>
  distinct() |>
  ungroup() |>
  mutate(elig_cancer = cancer_date <= earliest_end)  |>
  filter(elig_cancer) |> #keep only eligible cancers
  pivot_longer(cols = organ_system, values_to = 'cancer_organsys') #long format

table(cancer_sys$cancer_organsys, useNA = "ifany")
table(cancer_sys$case, cancer_sys$elig_cancer, dnn = c("case", "elig organ system"), useNA = "ifany")
table(cancer_sys$cancer_organsys, cancer_sys$elig_cancer, dnn = c("organ system", "elig organ system"), useNA = "ifany")

# create "long" dataframe of patients repeated by the number of body regions there is
pat_long_sys <- final |> 
  select(epatid, case) |> 
  slice(rep(1:n(), each = max(organ_system$organ_sys_numb))) |>
  mutate(cancer_organsys = rep(1:max(organ_system$organ_sys_numb),  times = nrow(final))) #add column for body region number

pat_long_sys <- 
  pat_long_sys |>
  left_join(cancer_sys |> select(epatid, case, cancer_date, elig_cancer, cancer_organsys), by = c('epatid', 'case', 'cancer_organsys'))  |> #join with eligible cancer cases
  left_join(final, by = c('epatid', 'case'), multiple = "all") #join with patients info

#join cancer region info with all patients info
final_sys <-
  pat_long_sys |>
  mutate(elig_cancer = ifelse(is.na(elig_cancer), F, elig_cancer),
         af_outcome = !case & af & af_date <= earliest_end & (is.na(cancer_date) | af_date <= cancer_date),
         outcome_date = if_else(elig_cancer,
                                pmin(cancer_date, earliest_end, na.rm = T),
                                earliest_end),
         outcome = case_when(
           cancer_date == outcome_date ~ "cancer",
           af_date == outcome_date & !case ~ "af",
           death_date == outcome_date ~ "death",
           end_date == outcome_date ~ "end"
         ),
         outcome_detail=ifelse(outcome=="cancer",
                               cancer_organsys,
                               outcome),
         follow =  time_length(start %--% outcome_date, unit = "month")) |>
  select(-af) |>
  rename(cancer_outcome = elig_cancer)

final_sys |>
  group_by(case, outcome_detail) |>
  summarise(n = n()) |>
  mutate(prop = 100*n/sum(n)) |>
  pivot_wider(id_cols = case, names_from = outcome_detail, values_from = prop)

setdiff(final$epatid, final_sys$epatid) #all patients included
cancer_sys |> filter(!is.na(cancer_organsys)) |> group_by(epatid, case) |> summarise(n = n()) |> filter(n>1)

#keep only required columns
final_sys <-
  final_sys |>
  select(epatid, pracid, gender, age_start, case,
         diab, ht, smok, imd, alc,
         af_outcome, cancer_outcome, cancer_organsys, outcome, outcome_detail,
         start, af_date, cancer_date, death_date, outcome_date, follow)

saveRDS(final_sys, paste0(path_rw, "cancer_system.RDS"))

write(paste(Sys.time(), "End 8_outcomes.R script"), file = progress_file, append = T)

rm(list=ls()); gc()

