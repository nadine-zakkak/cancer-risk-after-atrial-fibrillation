# Fri 14/Apr/2023
# Initial preprocessing in 2_preprocess.R script
# Get earliest AF diagnosis for all patients

library(dplyr)
library(stringr)
library(lubridate)

# Load data and add start and end dates
patients <- readRDS("./data/patients_preprocess.rds")
patients_numb <- read.csv("./results/patients_numb.csv")

# 2) Identify date of potentially eligible 'AF' diagnosis ----
###### potentially eligible AF: no history of AF
# clinical cprd
clinical <- readRDS(file="./data/e_17_205r_extract_clinical_final_1_50_r.rds") |>
  filter(e_patid %in% patients$epatid) |>
  rename(epatid = e_patid)

clinical_cont <- readRDS(file="./data/e_17_205r_extract_clinical_final_51_81_r.rds") |>
  filter(e_patid %in% patients$epatid) |>
  rename(epatid = e_patid)

#referral cprd
referral <- readRDS(file="./data/e_17_205r_extract_referral_final_r.rds") |>
  filter(e_patid %in% patients$epatid) |>
  rename(epatid = e_patid)

# diag hes
diag <- readRDS(file="./data/e_hes_diagnosis_epi_17_205r_final.rds") |>
  filter(e_patid %in% patients$epatid) |>
  rename(epatid = e_patid) #1938(?)-2017 mainly 1997-2017

# CALIBER AF codes
lookup_af_cprd <- read.csv(file="./look_ups//CPRD_AF.csv", header=TRUE, sep=",")
colnames(lookup_af_cprd) <- tolower(colnames(lookup_af_cprd))
lookup_af_cprd <-
  lookup_af_cprd |>
  mutate(cat = str_extract(category, '\\d{1,2}')) |> #extract categ numb
  filter(cat != 7) #remove atrial flutter

lookup_af_hes <- read.csv(file="./look_ups//ICD_AF.csv", header=TRUE, sep=",")
colnames(lookup_af_hes) <- tolower(colnames(lookup_af_hes))


#Step 1: Identify AF diag codes in CPRD or HES
## CPRD
### clinical
clinical_af_cprd <-
  clinical |> 
  filter(medcode %in% lookup_af_cprd$medcode) |> #any AF medcode
  select(epatid, eventdate, medcode) |> 
  inner_join(lookup_af_cprd |> select(medcode, cat), by = "medcode") #category of AF code

rm(clinical)
gc()

clinical_af_cprd_cont <-
  clinical_cont |> 
  filter(medcode %in% lookup_af_cprd$medcode) |>
  select(epatid, eventdate, medcode) |>
  inner_join(lookup_af_cprd |> select(medcode, cat), by="medcode") #AF codes

rm(clinical_cont)
gc()

#merge the data frames together
clinical_af_cprd_merge <- rbind(clinical_af_cprd, clinical_af_cprd_cont)

rm(clinical_af_cprd, clinical_af_cprd_cont)
gc()

#may contain multiple rows for same patient (same patient, same date, diff diag)
clinical_af_cprd_merge |> select(-medcode) |> distinct() |> filter(is.na(eventdate)) #1147 with missing eventdates
clinical_af_cprd_merge |> select(-medcode) |> distinct() |> filter(cat %in% c(1,2,3), is.na(eventdate)) #114 with historical AF and missing eventdates

# do not count if ONLY have events with missing data 
clinical_af_hx <-
  clinical_af_cprd_merge |>
  select(-medcode) |> 
  distinct() |> 
  filter(cat %in% c(1,2,3)) |> #historical AF 
  group_by(epatid) |>
  filter(
    all(is.na(eventdate)) # all event dates are missing 
    | eventdate == min(eventdate, na.rm = T)) #or the earliest eventdate in case 1+ recorded event date
#earliest date if >1 ICD code/patient
nrow(clinical_af_hx); length(unique(clinical_af_hx$epatid)) #55218 with 54947 unique patients

# why do some patients appear more than once?
clinical_af_hx |> count(epatid) |> filter(n>1)
#explore some
clinical_af_hx |> filter(epatid %in% unique(clinical_af_hx$epatid[duplicated(clinical_af_hx$epatid)]))
# same date but different category

clinical_af_cprd_merge |> select(-medcode) |> distinct() |> filter(cat %in% c(4,5,6), is.na(eventdate)) #1033 AF diag with missing eventdates

clinical_af_diag <-
  clinical_af_cprd_merge |>
  select(-medcode) |> 
  distinct() |> 
  filter(cat %in% c(4,5,6)) |> #AF diag
  group_by(epatid) |>
  filter(all(is.na(eventdate)) | eventdate == min(eventdate, na.rm = T)) #earliest date if >1 diag

table(year(clinical_af_diag$eventdate), useNA = "always") #1916-2017 (>1000: 1990-2017) #618 missing

rm(clinical_af_cprd_merge)
gc()

### referral
referral_af_cprd <-
  referral |> 
  filter(medcode %in% lookup_af_cprd$medcode) |> #AF medcode
  select(epatid, eventdate, medcode) |>
  inner_join(lookup_af_cprd |> select(medcode, cat), by = "medcode") #AF categories

rm(referral)
gc()

#may contain multiple rows for same patient (same patient, same date, diff diag)
referral_af_cprd |> select(-medcode) |> distinct() |> filter(cat %in% c(1,2,3), is.na(eventdate)) #none

referral_af_hx <-
  referral_af_cprd |>
  select(-medcode) |>
  distinct() |>
  filter(cat %in% c(1,2,3)) |> #historical AF
  group_by(epatid) |>
  filter(all(is.na(eventdate)) 
         | eventdate == min(eventdate, na.rm = T)) #earliest record

referral_af_cprd |> select(-medcode) |> distinct() |> filter(cat %in% c(4,5,6), is.na(eventdate)) #none

referral_af_diag <-
  referral_af_cprd |>
  select(-medcode) |>
  distinct() |>
  filter(cat %in% c(4,5,6)) |> #diag AF
  group_by(epatid) |>
  filter(all(is.na(eventdate)) 
         | eventdate == min(eventdate)) #earliest record

table(year(referral_af_diag$eventdate), useNA = "always") #1971-2017 (>100: 1990-2017)

rm(referral_af_cprd)
gc()

# cprd all
#may contain multiple rows for same patient (same patient, same date, diff diag)
cprd_af_hx <- rbind(clinical_af_hx, referral_af_hx) #all hx from cprd

cprd_af_hx |> group_by(epatid) |> filter(all(is.na(eventdate))) #83 patients with only empty eventdates - keep in data for now

cprd_af_hx <-
  cprd_af_hx |> 
  group_by(epatid) |> 
  filter(all(is.na(eventdate))
         | eventdate == min(eventdate, na.rm = T)) |> #earliest eventdate
  distinct() #remove probable duplicates from clinical+referral (same patid, eventdate)

rm(clinical_af_hx, referral_af_hx)

cprd_af_diag <- rbind(clinical_af_diag, referral_af_diag) #all AF diag

cprd_af_diag |> group_by(epatid) |> filter(all(is.na(eventdate))) #613 patients with only empty eventdates - keep in data for now

cprd_af_diag <-
  cprd_af_diag |> 
  group_by(epatid) |>
  filter(all(is.na(eventdate)) | eventdate == min(eventdate, na.rm = T)) |> 
  distinct() #remove probable duplicates from clinical+referral

table(year(cprd_af_diag$eventdate), useNA = "always") #1916-2017 (>1000: 1990-2017) #613 missing

rm(clinical_af_diag, referral_af_diag)
gc()

## HES
diag_af_hes <-
  diag |>
  filter(grepl("I48.*", icd)) |>
  # select(epatid, epistart, icd) #AF ICD code
  select(epatid, spno, epistart, icd, d_order) #AF ICD code

table(year(diag_af_hes$epistart), useNA = "always") #1995-2017, mainly (>6000: 1997-2017) and none missing

rm(diag)
gc()

#may contain multiple rows for same patient (same patient, same date, diff diag)
diag_af_hes |> select(-icd) |> distinct() |> filter(is.na(epistart)) #none
diag_af_hes <-
  diag_af_hes |>
  select(-icd) |>
  distinct() |>
  group_by(epatid) |>
  filter(all(is.na(epistart))
         | epistart == min(epistart, na.rm = T)) #earliest diag date
#306,524 patients

# AF cases
af_diag <- patients 

#adding cprd diag
af_diag <-
  af_diag |>
  mutate(cprd = epatid %in% cprd_af_diag$epatid) |>
  left_join(cprd_af_diag |> 
              select(epatid, eventdate), by=c("epatid" = "epatid"), multiple = 'all') |>
  rename(cprd_date=eventdate) |>
  distinct() #to remove any possible duplicates from cprd_af_diag (same patient, same date)

#cprd hx
af_diag <-
  af_diag |>
  mutate(hx = epatid %in% cprd_af_hx$epatid) |>
  left_join(cprd_af_hx |>
              select(epatid, eventdate), by=c("epatid" = "epatid"), multiple = 'all') |>
  rename(hx_date=eventdate) |>
  distinct() #to remove potential dup from cprd_af_hx

#hes diag
af_diag <-
  af_diag |>
  mutate(hes = epatid %in% diag_af_hes$epatid) |>
  left_join(diag_af_hes |>
              select(epatid, epistart), by=c("epatid" = "epatid"), multiple = 'all') |>
  rename(hes_date=epistart) |>
  distinct() #to remove potential dup from diag_af_hes

length(unique(af_diag$epatid)) == nrow(af_diag)

af_diag |> filter(cprd, is.na(cprd_date), !hes) |> select(epatid, start, cprd, cprd_date, hx, hx_date, hes, hes_date) 
#175 only cprd diagnosis with missing cprd diagnosis date

# Set flag for AF diagnosis and get earliest AF diagnosis date
# af_diag <-
af_diag <-
  af_diag |>
  mutate(af_date = pmin(cprd_date, hes_date, na.rm=TRUE)) |> #earliest date between cprd, hes diag
  mutate(potelig_af = cprd | hes)


table(af_diag$potelig_af) #334,912
table(af_diag$potelig_af, af_diag$hx)
patients_numb <- rbind(patients_numb, list("has AF related record", "has AF related record", 
                                           nrow(af_diag), length(which(af_diag$potelig_af | af_diag$hx))))

af_diag |> filter(potelig_af & is.na(af_date))
#175 potelig_af but with MISSING AF_DATE - EXCLUDE FROM DATA BECAUSE SHOULD NOT BE USED IN CONTROL AS 'NO-AF'
af_diag <- af_diag |> filter(!(potelig_af & is.na(af_date)))
table(af_diag$potelig_af) #334,737 

patients_numb <- rbind(patients_numb, list("has AF diagnosis record with known date or AF history record", 
                                           "has AF diagnosis record with known date or AF history record", 
                                           nrow(af_diag), length(which(af_diag$potelig_af | af_diag$hx))))

# # 72 patients with this case - exclude for now on the basis of having a history of AF
af_diag |> filter(hx & is.na(hx_date) & potelig_af) |> select(epatid, cprd, cprd_date, hx, hx_date, hes, hes_date, af_date, potelig_af)

# Incl/Excl based on AF records (hx/diag) and AF diagnosis date ----
table(af_diag$hes, af_diag$cprd, dnn  = c("hes", "cprd")) #180,261 in hes not in cprd, 28,213 in cprd not hes and 126,263 in both

table(af_diag$hx, af_diag$potelig_af, dnn  = c("history", "pheno")) #4,709 only with historical terms

#remove patients only with hx or monitoring terms
af_diag <-
  af_diag |>
  filter(!(hx & !potelig_af)) 

patients_numb <- rbind(patients_numb, list("only hx", "only medcodes for history of AF", 
                                           nrow(af_diag), length(which(af_diag$potelig_af))))

#remove patients with hx or monitoring before diag
af_diag <-
  af_diag |>
  filter(!(potelig_af & hx & hx_date <= af_date)) 
table(af_diag$hes, af_diag$cprd, dnn  = c("hes", "cprd")) #172,598 in hes not in cprd, 27,138 in cprd not hes and 120,332 in both

patients_numb <- rbind(patients_numb, list("hx<diag or hx unknown, diag known", "history of AF record before diagnosis record or history of AF date not recorded and has AF diagnosis record", 
                                           nrow(af_diag), length(which(af_diag$potelig_af))))

# drop unused columns
af_diag <-
  af_diag |> 
  select(-c(cprd, cprd_date, hx, hx_date, hes, hes_date))

write.csv(patients_numb, file="./results/patients_numb.csv", row.names = FALSE)
saveRDS(af_diag, "./data/patients_potelig_af.rds") #Potentially eligible AF flag
