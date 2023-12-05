# Fri 14/Apr/2023
# Initial preprocessing in 2_preprocess.R script and 3_pheno_AF.R
library(stringr)
library(dplyr)
library(lubridate)
library(ggplot2)

# Load data and add start and end dates
patients <- readRDS("./data/patients_potelig_af.rds")
controls_numb <- read.csv("./results/patients_numb.csv")

patients <-
  patients |>
  mutate(start_control = pmax(as.Date("1998-01-01", format = '%Y-%m-%d'), crd %m+% years(1), uts, age_18, na.rm=TRUE), #start of follow-up
         end_control = pmin(as.Date("2015-12-31"), dor, na.rm=TRUE)) #end of follow-up will be earliest of 31/12/2015, death, AF or outcome


# 1) Identify  eligible 'AF'  ----
# Eligible if AF diagnosis after start of follow-up

# keep only eligible AF cases (after start date) and patients without AF
af_diag <- 
  patients |>
  mutate(elig_af = potelig_af & af_date > start_control) |>
  filter(elig_af | !potelig_af) 

nrow(af_diag); length(which(af_diag$elig_af))
#ended with: 6356852 patients (223,221 with AF)

controls_numb <- rbind(controls_numb, list("diag before or on start of follow-up", "AF diagnosis before or on start of follow-up", 
                                           nrow(af_diag), length(which(af_diag$elig_af))))



write.csv(controls_numb, file="./results/controls_numb.csv", row.names = FALSE)

saveRDS(af_diag, "./data/controls_patients.rds")

years <- year(af_diag$af_date)
years <- years[!is.na(years)]
boxplot(years)
png("./results/af_years_controls.png")
ggplot() + geom_bar(aes(x = year(af_date)), data = af_diag) + labs(x = "Year of AF diagnosis")
dev.off()

# 2) exclude based on history of cancer -----
rm(list=ls()); gc()
library(dplyr)
library(lubridate)

controls_numb <- read.csv("./results/controls_numb.csv", stringsAsFactors = FALSE)
patients <- readRDS("./data/controls_patients.rds")
nrow(patients); length(which(patients$elig_af)) #6,356,852 patients (223,221 with AF)

cancer_full <- readRDS("./data/e_17_205_cancer_registration_1_final.rds")
cancer_full <- cancer_full |> mutate(site_icd10_o2 = trimws(site_icd10_o2))
cancer_full <- as.data.frame(cancer_full) |> rename(epatid = e_patid) |> filter(epatid %in% patients$epatid)

#icd-10 codes
cancer_pheno <- read.delim("./phenotypes/cancersite_v1.2_codelist_4digicd10.txt")
setdiff(cancer_full$site_icd10_o2, cancer_pheno$icd10_4dig) #C941, C141

#icd-9 codes
neoplasm_icd9 <- c(paste("^",setdiff(140:208, 173), ".*", sep=""), "^2252$", "^2254$")
neoplasm_icd9_regex <- paste(neoplasm_icd9, collapse="|")

cancer_icd9 <-
  cancer_full |>
  filter(grepl(neoplasm_icd9_regex, site_coded))

cancer_icd10 <- # only cancerous ICD-10 codes
  cancer_full |>
  inner_join(cancer_pheno |> 
               filter(cancer_site_number != 0) |>
               select(icd10_4dig, cancer_site_desc, cancer_group_desc),
             by = c("site_icd10_o2"="icd10_4dig"))

cancer_filter <- cancer_icd9 |> bind_rows(cancer_icd10)

# earliest cancer record for each patient 
cancer_filter |> filter(is.na(diagnosisdatebest)) #no missing diagnosis dates

cancer_filter_earliest <- cancer_filter |> 
  group_by(epatid) |> 
  filter(diagnosisdatebest == min(diagnosisdatebest)) |>
  distinct(epatid, diagnosisdatebest)

# did patient have a cancer history i.e. cancer on or before start of follow-up?
cancer_patients_hx <-
  cancer_filter_earliest |>
  left_join(patients |> select(epatid, start_control), by = "epatid") |>
  mutate(cancer_hx = diagnosisdatebest <= start_control) |> 
  filter(cancer_hx)

#get patient id of patients with a cancer history
cancer_hx_patid <- unique(cancer_patients_hx$epatid) #134,846 

patients <-
  patients |>
  mutate(cancer_hx = epatid %in% cancer_hx_patid)
table(patients$cancer_hx) #134,846 (2.1%) of all patients had a cancer history 
100*prop.table(table(patients$cancer_hx))

# keep only patients without a cancer history 
patients <-
  patients |>
  filter(!cancer_hx)
table(patients$cancer_hx)
table(patients$elig_af) #215,143 with AF from 223,221 originally

library(ggplot2)
png("./results/af_years_controls_excl_cancer.png")
ggplot() + geom_bar(aes(x=year(af_date)), data=patients) + labs(x = "Year of AF diagnosis")
dev.off()

controls_numb <- rbind(controls_numb, list("cancer history", "had a cancer history", 
                                           nrow(patients), length(which(patients$elig_af))))

write.csv(controls_numb, file="./results/controls_numb.csv", row.names = FALSE)
saveRDS(patients, "./data/controls_patients.rds")

#3)  Save earliest cancer date   ------
rm(list=ls()); gc()
controls_numb <- read.csv("./results/controls_numb.csv", stringsAsFactors = FALSE)
patients <- readRDS("./data/controls_patients.rds")

cancer <- readRDS("./data/e_17_205_cancer_registration_1_final.rds")
cancer <- cancer |> mutate(site_icd10_o2 = trimws(site_icd10_o2))

# codelists
source("./scripts/00_read_codelists.R")
cancer_codelist_outcome <- cancer_codelist |> filter(outcome == 1)

cancer_filter <-
  cancer |>
  rename(epatid = e_patid) |>
  mutate(cancer = site_icd10_o2 %in% cancer_codelist_outcome$icd10_4dig) |>
  filter(cancer)

# earliest cancer record for each patient 
cancer_filter |> filter(is.na(diagnosisdatebest)) #none missing
cancer_filter_earliest <- cancer_filter |> 
  group_by(epatid) |> 
  filter(diagnosisdatebest == min(diagnosisdatebest)) |>
  distinct(epatid, diagnosisdatebest) 

cancer_patients <-
  cancer_filter_earliest |>
  left_join(patients |> select(epatid, start_control), by = "epatid")

cancer_patients |> filter(diagnosisdatebest <= start_control) #should be none because previously excluded based on 'history' of cancer

cancer_patients <-
  cancer_patients |>
  mutate(cancer_diag = diagnosisdatebest > start_control)

table(cancer_patients$cancer_diag) #All true as it should be (264,402)

patients <-
  patients |>
  left_join(cancer_patients |> select(epatid, cancer_date=diagnosisdatebest, cancer=cancer_diag),
            by = 'epatid')

table(patients$cancer) # 264,402 

write.csv(controls_numb, file="./results/controls_numb.csv", row.names = FALSE)
saveRDS(patients, "./data/controls_patients.rds")

#4) Set new end date -----
rm(list=ls()); gc()
controls_numb <- read.csv("./results/controls_numb.csv", stringsAsFactors = FALSE)
patients <- readRDS("./data/controls_patients.rds")

patients <-
  patients |> 
  mutate(end_control_fu = pmin(as.Date("2015-12-31"), dor, af_date, cancer_date, na.rm=TRUE))

patients <-
  patients |>
  mutate(invalid_fu = end_control_fu <= start_control)

patients |> filter(invalid_fu) |> select(epatid, yob, crd, tod, deathdate, lcd, uts, dor, age_18, age_100, start_control, cancer_date, end_control_fu)

final_patients <-
  patients |>
  filter(!invalid_fu)

controls_numb <- rbind(controls_numb, list("end<=start", "max(1998,crd+1y,uts,18yo) to min(2015,dod,af,cancer)",
                                           nrow(final_patients), length(which(final_patients$elig_af))))
controls_numb <- controls_numb |> mutate(numb_diff = lag(numb)-numb, AF_diff = lag(AF)-AF)

write.csv(controls_numb, file="./results/controls_numb.csv", row.names = FALSE)
saveRDS(final_patients, "./data/controls_patients_final.rds")

rm(list=ls()); gc()

