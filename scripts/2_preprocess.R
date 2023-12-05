# Nadine Zakkak
# Fri 14/April/2023
# 
# Exclude if patid identified to be the same person via CR
# Join with pracid, death data from ONS
# 
# Steps:
# 1) Identify start of follow-up based on study dates ----
## Case
###### start: max(18yo, uts, crd+1y, 1998)
###### end of CPRD follow-up: min(100yo, tod, dod, lcd, 2014) 
###### censor at 2015, death or cancer/outcome
## Control
###### start: max(18yo, uts, crd+1y, 1998)
###### end of CPRD follow-up: min(dod, 2015) 
###### censor at 2015, AF, death or cancer/outcome
library(stringr)
library(dplyr)
library(lubridate)

#QC already applied by CPRD:
##"acceptable"
##eligible for linkage to HES and NCIN
##Gender: M or F
##Age: 18+ (at start or during study)
##Registered during 1/1/1998-31/5/2016

cprd_patients <- readRDS("./data/e_17_205r_extract_patient_final_r.rds") #original

#Extract pracid and patid from combined patid (epatid)
patients <-
  
  cprd_patients |>
  mutate(pracid = str_extract(epatid, '\\d{1,3}$')) |>
  mutate(pracid = as.character(str_remove(pracid, "^0+"))) |>
  select(epatid, pracid, everything())

all_prac <- readRDS(file = "./data/e_17_205r_extract_practice_final_r.rds")
all_prac$epracid <- as.character(all_prac$epracid)

patients <-
  patients |>
  left_join(all_prac, by = c("pracid"="epracid"))

death <- readRDS(file="./data/e_death_patient_17_205r_final.rds")

patients <-
  patients|>
  left_join(death |> select(e_patid, dor), by=c('epatid'='e_patid')) 

patients <-
  patients |>
  mutate(age_18 = as.Date(paste0(yob, "-06-15")) %m+% years(18), #age calculated at ~mid-year
         age_100 = as.Date(paste0(yob, "-06-15")) %m+% years(100),
         start = pmax(as.Date("1998-01-01"), crd %m+% years(1), uts, age_18, na.rm=TRUE)) # start date is the same for cases and controls

patients_numb <- data.frame(summary="original", description="original",
                         numb=nrow(patients), AF=NA, stringsAsFactors = FALSE) #to save pop numbers


# 2) Exclude patients that share the same tumour ID identified from cancer registry ------
cancer_full <- readRDS("./data/e_17_205_cancer_registration_1_final.rds")

cancer_tumour_mult <- cancer_full |>
  group_by(e_cr_id) |>
  filter(n()>1) #all rows with repeated tumours

cancer_mult_id <- unique(cancer_tumour_mult$e_patid) #patients with at least 1 repeated tumour with another patid (56,130 patid)

patients <- patients |> filter(!(epatid %in% cancer_mult_id)) #remove patients with 1+ repeated tumour with another patid

patients_numb <- rbind(patients_numb, list("diff pat id traced back to same person",
                                     "diff patid share same tumour id so all patid excluded",
                                     nrow(patients), NA))

write.csv(patients_numb, file="./results/patients_numb.csv", row.names = FALSE)
saveRDS(patients, "./data/patients_preprocess.rds") #exclusions based on initial criteria
