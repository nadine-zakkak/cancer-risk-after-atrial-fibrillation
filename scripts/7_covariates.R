# Nadine Zakkak
# Identify covariates to be included in analysis for cases and controls

rm(list=ls()); gc()

library(dplyr)
library(stringr)
library(readr)

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 7_confounders.R script"), file = progress_file, append = T)

# read files
cases <- readRDS("./data/cases.rds")
controls <- readRDS("./data/controls.rds")
allpatients <- #join into 1 dataframe
  cases |>
  mutate(case = T) |>
  select(epatid, pracid, case, af=elig_af, gender, age_start, 
         start=start_case_fu, death_date=dor, af_date) |>
  bind_rows(
    controls |>
      mutate(case = F) |>
      select(epatid, pracid, case, af=elig_af, gender, age_start, 
             start=start_fu, 
             death_date=dor, af_date, end_control_fu)
  )

clinical <- readRDS("./data/final_data/clinical_final_cohort.rds") |> rename(epatid=e_patid)
referral <- readRDS("./data/final_data/referral_final_cohort.rds") |> rename(epatid=e_patid)
additional <- readRDS("./data/final_data/additional_final_cohort.rds") |> rename(epatid = e_patid)
therapy <- readRDS("./data/final_data/therapy_final_cohort.rds") |> rename(epatid = e_patid)
diag <- readRDS("./data/final_data/hes_diagnosis_epi_final_cohort.rds") |> rename(epatid = e_patid)
diag$icd_mod <- str_replace_all(diag$icd, "\\.", "") #replace all "." by nothing. (Ex: C44.1 --> C441)

# -------------------- DIABETES ------------------------- 
# -3: insulin dependent, 4:non-insulin dependent, 6:NOS
# -3: type I, 4: type II
write(paste(Sys.time(), "Getting diabetes status"), file = progress_file, append = T)

lookup_dm_gprd <- read.csv("./look_ups/dm_cprd_codelist_2.csv", header=TRUE, sep=",")
lookup_dm_gprd <- lookup_dm_gprd |> filter(category %in% c(3,4, 6)) |> select(category, medcode, readterm , readcode) 
lookup_dm_hes <- read.csv("./look_ups/dm_hes_codelist_1.csv", header=TRUE, sep=",")
lookup_dm_hes <- lookup_dm_hes |> filter(category %in% c(3,4,6)) |> select(category, icd_code, icd_term) 
lookup_diab_cprd <- read.csv("./look_ups/diabdiag_cprd_codelist_1.csv", header=TRUE, sep=",")
lookup_diab_cprd <- lookup_diab_cprd |> filter(category %in% c(3,4)) |> select(category, medcode, readterm, readcode)

# Get all events from CPRD
clinical_dm <- clinical |> filter(medcode %in% c(lookup_dm_gprd$medcode, lookup_diab_cprd$medcode))
referral_dm <- referral |> filter(medcode %in% c(lookup_dm_gprd$medcode, lookup_diab_cprd$medcode))

cprd_diab_df <- clinical_dm |> select(epatid, eventdate, medcode) |>
  bind_rows(referral_dm |> select(epatid, eventdate, medcode)) |>
  distinct()

# Only keep events that occurred before the patients start date
cprd_diab_df <- 
  cprd_diab_df |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate == min(eventdate, na.rm = T)) |>
  filter(eventdate <= start) |>
  distinct(case, epatid, eventdate) |>
  select(epatid,case, everything()) |>
  mutate(cprd_diab = T) |>
  ungroup()

# Add cprd_diab flag to main dataframe
all_patients <- allpatients |> 
  left_join(cprd_diab_df |> select(epatid, case, cprd_diab, cprd_eventdate=eventdate), by = c('epatid', 'case')) |>
  mutate(cprd_diab = ifelse(is.na(cprd_diab), F, cprd_diab))

# Get all events from HES DIAG
regex_pattern <- paste(paste0(lookup_dm_hes$icd_code, ".*"), collapse="|")
hes_diab_df <- diag |> filter(grepl(regex_pattern, icd_mod))

# Only keep events that occurred before the patient's start date
hes_diab_df <- 
  hes_diab_df |> 
  inner_join(allpatients |> select(epatid, case, start), by = 'epatid') |>
  group_by(case, epatid) |>
  filter(epistart == min(epistart, na.rm = T)) |>
  filter(epistart <= start) |>
  distinct(case, epatid, epistart) |>
  select(epatid,case, everything()) |>
  mutate(hes_diab = T) |>
  ungroup()

# Add hes_diab flag to main dataframe
all_patients <- all_patients |> 
  left_join(hes_diab_df |> select(epatid, case, hes_diab, hes_eventdate=epistart), by = c('epatid', 'case')) |>
  mutate(hes_diab = ifelse(is.na(hes_diab), F, hes_diab))

#final - flag as having diabetes of either had cprd_diab or hes_diab
all_patients <- all_patients |> mutate(diab = cprd_diab | hes_diab) #final diabetes

all_patients <- all_patients |> select(-cprd_diab, -hes_diab)
all_patients <- all_patients |> select(-cprd_eventdate, -hes_eventdate)

rm(clinical_dm, cprd_diab_df, hes_diab_df, referral_dm, 
   lookup_diab_cprd, lookup_dm_gprd, lookup_dm_hes, regex_pattern)

# -------------------- HYPERTENSION -------------------------
write(paste(Sys.time(), "Getting hypertension status"), file = progress_file, append = T)

# -3: ht, 4:secondary ht
lookup_ht_cprd <-  read.csv("./look_ups/ht_cprd_codelist_1.csv", header=TRUE, sep=",")
lookup_ht_cprd <- lookup_ht_cprd |> filter(category %in% c(3,4)) |> select(category, medcode, readterm , readcode)
lookup_ht_hes <- read.csv("./look_ups/ht_hes_codelist_1.csv", header=TRUE, sep=",")
lookup_ht_hes <- lookup_ht_hes |> filter(category %in% c(3,4)) |> select(category, icd_code, icd_term)

# Get all events from CPRD
clinical_ht <- clinical |> filter(medcode %in% lookup_ht_cprd$medcode)
referral_ht <- referral |> filter(medcode %in% lookup_ht_cprd$medcode)

cprd_ht_df <- clinical_ht |> select(epatid, eventdate, medcode) |>
  bind_rows(referral_ht |> select(epatid, eventdate, medcode)) |>
  distinct()

# Only keep events that occurred before the patient's start date
cprd_ht_df <- 
  cprd_ht_df |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate == min(eventdate, na.rm = T)) |>
  filter(eventdate <= start) |>
  distinct(case, epatid, eventdate) |>
  select(epatid,case, everything()) |>
  mutate(cprd_ht = T) |>
  ungroup()

# add flag cprd_ht to main dataframe
all_patients <- all_patients |> 
  left_join(cprd_ht_df |> select(epatid, case, cprd_ht, cprd_eventdate=eventdate), by = c('epatid', 'case')) |>
  mutate(cprd_ht = ifelse(is.na(cprd_ht), F, cprd_ht))

table(all_patients$cprd_ht, useNA = "ifany")

# Get all events from HES diagnosis
regex_pattern <- paste(paste0(lookup_ht_hes$icd_code, ".*"), collapse="|")
hes_ht_df <- diag |> filter(grepl(regex_pattern, icd_mod))
# Only keep events that occurred before the patient's start date
hes_ht_df <- 
  hes_ht_df |> 
  inner_join(allpatients |> select(epatid, case, start), by = 'epatid') |>
  group_by(case, epatid) |>
  filter(epistart == min(epistart, na.rm = T)) |>
  filter(epistart <= start) |>
  distinct(case, epatid, epistart) |>
  select(epatid,case, everything()) |>
  mutate(hes_ht = T) |>
  ungroup()

# add flag hes_ht to main dataframe
all_patients <- all_patients |> 
  left_join(hes_ht_df |> select(epatid, case, hes_ht, hes_eventdate=epistart), by = c('epatid', 'case')) |>
  mutate(hes_ht = ifelse(is.na(hes_ht), F, hes_ht))

##final - flag as having hypertension of either had cprd_ht or hes_ht
all_patients <- all_patients |> mutate(ht = cprd_ht | hes_ht) #final hypertension

all_patients <- all_patients |> select(-cprd_ht, -hes_ht)
all_patients <- all_patients |> select(-cprd_eventdate, -hes_eventdate)

rm(clinical_ht, referral_ht, cprd_ht_df, hes_ht_df,
   lookup_ht_cprd, lookup_ht_hes, regex_pattern)

# -------------------- SMOKING -------------------------
write(paste(Sys.time(), "Getting smoking status"), file = progress_file, append = T)

# assuming non-smokers if no record of smoking
#categories
#-- cprd: 1-non 2-ex 3-ever 4-current
#-- hes: 4-current
#-- additional: 1-yes, 2-no, 3-ex
lookup_smk_cprd <- read.csv("./look_ups/smoking_status_cprd_codelist_1.csv")
lookup_smk_cprd[13,"category"] <- 1 #changed "current non-smoker" to category 1 manually - check if keep that way or no
lookup_smk_hes <- read.csv("./look_ups/smoking_status_hes_codelist_2.csv")
lookup_smk_hes <- lookup_smk_hes |> filter(!is.na(category)) |> select(icd_code,icd_term, category) |>
  mutate(icd_regex = paste0(icd_code, ".*"))
regex_pattern <- paste(paste0(lookup_smk_hes$icd_code, ".*"), collapse="|")
lookup_bnf <- read_delim("./look_ups/bnfcodes.txt", col_types=c("ic"), delim="\t")
bnf_regex = "^041002.*" #nicotine dependence 
bnf_smk <- lookup_bnf |> filter(grepl(bnf_regex, bnf)) |> select(bnfcode) |> pull()

#record of smoking
# -2: ex, 3:ever, 4: current
lookup_smk_cprd_ever <- lookup_smk_cprd |> filter(category %in% c(2,3,4))
lookup_smk_hes_ever <- lookup_smk_hes |> filter(category %in% c(2,3,4))

# Same code structure as above
#smok in cprd
clinical_smk <- clinical |> filter(medcode %in% lookup_smk_cprd_ever$medcode)
referral_smk <- referral |> filter(medcode %in% lookup_smk_cprd_ever$medcode)

cprd_smk_df <- clinical_smk |> select(epatid, eventdate, medcode) |>
  bind_rows(referral_smk |> select(epatid, eventdate, medcode)) |>
  distinct()

cprd_smk_df <- 
  cprd_smk_df |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate == min(eventdate, na.rm = T)) |>
  filter(eventdate <= start) |> #get all events on or before start date
  distinct(case, epatid, eventdate) |>
  select(epatid,case, everything()) |>
  mutate(cprd_smk = T) |>
  ungroup()

all_patients <- all_patients |>
  left_join(cprd_smk_df |> select(epatid, case, cprd_smk, cprd_eventdate=eventdate), by = c('epatid', 'case')) |>
  mutate(cprd_smk = ifelse(is.na(cprd_smk), F, cprd_smk))

#diag code in hes
hes_smk_df <- diag |> filter(grepl(regex_pattern, icd_mod))
hes_smk_df <- 
  hes_smk_df |> 
  inner_join(allpatients |> select(epatid, case, start), by = 'epatid') |>
  group_by(case, epatid) |>
  filter(epistart == min(epistart, na.rm = T)) |>
  filter(epistart <= start) |>
  distinct(case, epatid, epistart) |>
  select(epatid,case, everything()) |>
  mutate(hes_smk = T) |>
  ungroup()

all_patients <- all_patients |> 
  left_join(hes_smk_df |> select(epatid, case, hes_smk, hes_eventdate=epistart), by = c('epatid', 'case')) |>
  mutate(hes_smk = ifelse(is.na(hes_smk), F, hes_smk))

# remove diag to free up space
rm(diag); gc()

#smok records by manual input (enttype=4)
clinical_ent_smk <- clinical |> filter(enttype == 4 & adid != 0)
add_smk_df <- additional |>
  inner_join(clinical_ent_smk |> select(-enttype), by = c('epatid', 'adid')) |>
  filter(data1 %in% c(1,3)) |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate == min(eventdate, na.rm = T)) |>
  filter(eventdate <= start) |>
  distinct(case, epatid, eventdate) |>
  select(epatid,case, everything()) |>
  mutate(add_smk = T) |>
  ungroup()

all_patients <- all_patients |> 
  left_join(add_smk_df |> select(epatid, case, add_smk, add_eventdate=eventdate), by = c('epatid', 'case')) |>
  mutate(add_smk = ifelse(is.na(add_smk), F, add_smk))

length(unique(add_smk_df$epatid)) #89595
table(all_patients$add_smk, useNA = "ifany")

#smok cessation prescriptions
therapy_smk_df <- therapy |> filter(bnfcode %in% bnf_smk)
therapy_smk_df <- therapy_smk_df |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate == min(eventdate, na.rm = T)) |>
  filter(eventdate <= start) |>
  distinct(case, epatid, eventdate) |>
  select(epatid,case, everything()) |>
  mutate(ther_smk = T) |>
  ungroup()

# therapy df no longer needed - remove to clear up space
rm(therapy); gc()

all_patients <- all_patients |> 
  left_join(therapy_smk_df |> select(epatid, case, ther_smk, ther_eventdate=eventdate), by = c('epatid', 'case')) |>
  mutate(ther_smk = ifelse(is.na(ther_smk), F, ther_smk))


#final
all_patients <- all_patients |> mutate(smok = cprd_smk | hes_smk | add_smk | ther_smk) #final smok

all_patients <- all_patients |> select(-c(cprd_smk, hes_smk, add_smk, ther_smk)) 
all_patients <- all_patients |> select(-c(cprd_eventdate, hes_eventdate, add_eventdate, ther_eventdate)) 

rm(clinical_smk, add_smk_df, cprd_smk_df, hes_smk_df, referral_smk, therapy_smk_df,
   clinical_ent_smk,
   regex_pattern, bnf_regex, bnf_smk, 
   lookup_smk_cprd, lookup_smk_cprd_ever, lookup_smk_hes, lookup_smk_hes_ever,
)

# -------------------- IMD -------------------------
write(paste(Sys.time(), "Getting imd"), file = progress_file, append = T)

imd <- readRDS("./data/final_data/practice_imd_final_cohort.rds")
imd <- imd |> mutate(e_practid = as.character(e_practid))

all_patients <-
  all_patients |>
  left_join(imd |> select(-country, pracid=e_practid, imd=e2015_imd_5), by='pracid')

table(all_patients$imd, useNA = "ifany")

rm(imd)

# -------------------- ALCOHOL -------------------------
write(paste(Sys.time(), "Getting alcohol drinking status"), file = progress_file, append = T)

# get most recent status before start date 
# split into :
# current, ex, excess, never
# maybe only for digestive? as a "second" analysis just to see have 
# 1:non, 2: ex, 3:occ, 4:current, 5:excess, 6:binge
# 1,2,3: drinker/current; 4,5: excess
# additional: 1: yes, 2:no, 3:ex

# New categories:
# 0: Missing
# 1: non-drinker
# 2: Ex drinker
# 3: Current drinker
lookup_alc_drink <- read.csv("./look_ups/alcohol_drinker_cprd_codelist_1.csv")
lookup_alc_unit <- read.csv("./look_ups/alcohol_units_cat_cprd_codelist_1.csv")
clinical_ent_alc <- clinical |> filter(enttype==5 & adid !=0 & medcode != 97126)

lookup_alc_dict <- data.frame(
  category_name = c("Missing", "Non", "Ex", "Current"),
  category_num = c(0,1,2,3))

write.csv(lookup_alc_dict, "./look_ups/alcohol_dictionary.csv")

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

lookup_alc_cprd <- lookup_alc_drink |> select(medcode, category_new) |>
  bind_rows(lookup_alc_unit |> select(medcode, category_new)) |>
  filter(!is.na(category_new))

lookup_add_ent <- data.frame(ent = as.character(c(1,2,3)), 
                             desc = c("yes", "no", "ex"),
                             category_new = c(3,1,2))

# get all at once 
clinical_alc_df <-
  clinical |>
  inner_join(lookup_alc_cprd, by = 'medcode') |>
  inner_join(allpatients |> select(case, epatid, start), by = 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate <= start) |> #before or on start date
  filter(eventdate == max(eventdate, na.rm = T)) |> #latest record before start date
  distinct(case, epatid, eventdate, category_new) |>
  select(epatid, case, everything()) |>
  ungroup()

alc_df <- allpatients |>
  select(epatid, case, start) |>
  left_join(clinical_alc_df |> select(epatid, case, clinical_alc=category_new, clinical_alc_eventdate=eventdate),
            by = c('epatid', 'case')) |> 
  distinct()

referral_alc_df <-
  referral |>
  inner_join(lookup_alc_cprd, by = 'medcode') |>
  inner_join(allpatients |> select(case, epatid, start), by = 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate <= start) |> #before or on start date
  filter(eventdate == max(eventdate, na.rm = T)) |> #latest record before start date
  distinct(case, epatid, eventdate, category_new) |>
  select(epatid, case, everything()) |>
  ungroup()

alc_df <- alc_df |>
  left_join(referral_alc_df |> select(epatid, case, referral_alc=category_new, referral_alc_eventdate=eventdate),
            by = c('epatid', 'case')) |> 
  distinct()

add_alc_df <- additional |>
  inner_join(clinical_ent_alc, by=c('epatid', 'adid')) |> #link with clinical alc data
  inner_join(lookup_add_ent |> select(ent, category_new), by = c('data1'='ent')) |>
  inner_join(allpatients |> select(case, epatid, start), by= 'epatid') |>
  group_by(case, epatid) |>
  filter(eventdate <= start) |>
  filter(eventdate == max(eventdate, na.rm = T)) |>
  distinct(case, epatid, eventdate, category_new) |>
  select(epatid,case, everything()) |>
  ungroup()

alc_df <- alc_df |>
  left_join(add_alc_df |> select(epatid, case, add_alc=category_new, add_alc_eventdate=eventdate),
            by = c('epatid', 'case')) |> 
  distinct()

#clinical file diagnoses takes priority
alc_df <- alc_df |> 
  mutate(alc_date = pmax(clinical_alc_eventdate, referral_alc_eventdate, add_alc_eventdate, na.rm = T)) |>
  mutate(alc = case_when(
    alc_date == clinical_alc_eventdate ~ clinical_alc,
    alc_date == referral_alc_eventdate ~ referral_alc,
    alc_date == add_alc_eventdate ~ add_alc
  )) |> 
  mutate(alc = ifelse(is.na(alc), 0, alc))
  distinct()

#clean df & merge same patients together
alc_df_clean <- 
  alc_df |> select(-c(clinical_alc, clinical_alc_eventdate, 
                      referral_alc, referral_alc_eventdate, 
                      add_alc, add_alc_eventdate)) |>
  distinct()

#if has >1 category on same date, always take the max one: 3-Current drinker, 2-Ex drinker,  1-non-drinker, 0-Missing  
alc_df_single <-
  alc_df_clean |>
  group_by(epatid, case, alc_date) |>
  summarise(alc = max(alc, na.rm = T)) |>
  ungroup()

write.csv(100*prop.table(table(alc_df_single$case, alc_df_single$alc, useNA = "ifany"), margin = 1), "alc.csv")

all_patients <- all_patients |>
  left_join(alc_df_single, by=c('epatid', 'case'))

rm(add_alc_df, alc_df, alc_df_clean, alc_df_single, clinical_alc_df, clinical_ent_alc, referral_alc_df,
   lookup_add_ent, lookup_alc_cprd, lookup_alc_dict, lookup_alc_drink, lookup_alc_unit)


# -------------------- GENDER -------------------------
# recoding gender: 0- male, 1-female
all_patients <- all_patients |> mutate(gender = gender - 1) 

# --------- final cohort w/o outcomes w/ missing data --------
# -- patid
# -- case
# -- start date
# -- all confounders
# -- end of af follow-up date (need to be adjusted accordingly for outcomes)

saveRDS(all_patients, "./data/final_data/confounders_final.rds")

write(paste(Sys.time(), "Confounders file saved"), file = progress_file, append = T)

write(paste(Sys.time(), "End 7_confounders.R script"), file = progress_file, append = T)


rm(list=ls()); gc()
source("./scripts/8_outcomes.R")



