# Fri 14/Apr/2023
# Initial preprocessing in 2_preprocess.R script and in 3_pheno_AF.R

library(stringr)
library(dplyr)
library(lubridate)
library(ggplot2)

# Load data and add start and end dates
patients <- readRDS("./data/patients_potelig_af.rds")
cases_numb <- read.csv("./results/patients_numb.csv")

patients <-
  patients |>
  mutate(start_case = pmax(as.Date("1998-01-01"), crd %m+% years(1), uts, age_18, na.rm=TRUE), #start of follow-up
         end_case = pmin(as.Date("2014-12-31"), tod, lcd, dor, age_100, na.rm=TRUE))

#5/May adaptations ----
# Eligible if AF diagnosis within CPRD follow-up period & not on death date
af_diag <-
  patients |>
  mutate(elig_af = potelig_af
         # within CPRD follow-up period (including lower and upper bounds)
         & af_date >= start_case & af_date <= end_case
         # didn't die or died but AF before death date
         & (is.na(dor) | af_date < dor)) |>
  # keep those with eligible AF and those without AF
  filter(elig_af | !potelig_af) |>
  distinct()
# End adaptations

cases_numb <- rbind(cases_numb, list("diag not in follow-up", "AF diagnosis was not within CPRD follow-up period", 
                                     nrow(af_diag), length(which(af_diag$potelig_af))))

nrow(af_diag); length(which(af_diag$elig_af))
#ended with: 6,266,905 patients (133,274 with AF)

#keep only AF cases
af_diag <-
  af_diag |>
  filter(elig_af) |>
  #set start of follow-up for case to be AF date
  mutate(start_case_fu = as.Date(af_date))

write.csv(cases_numb, file="./results/cases_numb.csv", row.names = FALSE)
saveRDS(af_diag, "./data/cases_patients_af.rds") #no exclusions based on hx of cancer 

years <- year(af_diag$af_date)
years <- years[!is.na(years)]
boxplot(years)
png("./results/af_years.png")
ggplot() + geom_bar(aes(x = year(af_date)), data = af_diag) + labs(x = "Year of AF diagnosis")
dev.off()


rm(list = ls())
gc()
# 2) exclude based on cancer history -----
library(dplyr)
library(lubridate)
library(ggplot2)

cases_numb <- read.csv("./results/cases_numb.csv", stringsAsFactors = FALSE)
af_diag <- readRDS("./data/cases_patients_af.rds") #133,354

cancer_full <- readRDS("./data/e_17_205_cancer_registration_1_final.rds")
cancer_full <- cancer_full |> mutate(site_icd10_o2 = trimws(site_icd10_o2))
cancer_full <- as.data.frame(cancer_full) |> rename(epatid = e_patid) |> filter(epatid %in% af_diag$epatid)

#icd-10 codes
cancer_pheno <- read.delim("./phenotypes/cancersite_v1.2_codelist_4digicd10.txt")
setdiff(cancer_full$site_icd10_o2, cancer_pheno$icd10_4dig) #NA
cancer_full |> filter(is.na(site_icd10_o2)) |> summarise(min = min(year(diagnosisdatebest)), max = max(year(diagnosisdatebest)))

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
cancer_filter |> filter(is.na(diagnosisdatebest)) #none missing
cancer_filter_earliest <- cancer_filter |> 
  group_by(epatid) |> 
  filter(diagnosisdatebest == min(diagnosisdatebest)) |>
  distinct(epatid, diagnosisdatebest)

cancer_af_patients <-
  cancer_filter_earliest |>
  # select(epatid, e_cr_patid, e_cr_id, diagnosisdatebest) |>
  left_join(af_diag |> select(epatid, start_case_fu), by = "epatid") |>
  mutate(cancer_hx = diagnosisdatebest <= start_case_fu) |>
  filter(cancer_hx)

af_diag <-
  af_diag |>
  mutate(cancer_hx = epatid %in% cancer_af_patients$epatid)
table(af_diag$cancer_hx) #16,101 (12.1%) of AF patients had a cancer history
100*prop.table(table(af_diag$cancer_hx))

af_diag <- 
  af_diag |> 
  filter(!cancer_hx)
table(af_diag$cancer_hx) 
nrow(af_diag) #117,173 AF patients from 133,274

library(ggplot2)
png("./results/af_years_final.png")
ggplot() + geom_bar(aes(x=year(start_case_fu)), data=af_diag) + labs(x = "Year of AF diagnosis")
dev.off()

cases_numb <- rbind(cases_numb, list("cancer history", "had a cancer history (AF patients)", 
                                     nrow(af_diag), length(which(af_diag$elig_af))))
cases_numb <- cases_numb |> mutate(numb_diff = lag(numb)-numb, AF_diff = lag(AF)-AF)

write.csv(cases_numb, file="./results/cases_numb.csv", row.names = FALSE)
saveRDS(af_diag, "./data/cases_patients_final.rds")



