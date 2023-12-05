# Sunday 21/April/23
# (6) 'Clean' files by only including patients identified for study 
# potentially don't need all files loaded here so could be removed

rm(list=ls()); gc()

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 6_clean_files.R script"), file = progress_file, append = T)

library(dplyr)
library(lubridate)
library(stats)

path_read <- "./data/"
path_save <- "./data/final_data/"
if(!dir.exists(path_save)) {dir.create(path_save)}

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
patid <- unique(allpatients$epatid)
pracid <- unique(allpatients$pracid)

# explore
cases_nd_controls <- controls |> filter(epatid %in% cases$epatid)
allpatients_case_control <- allpatients |> filter(epatid %in% cases_nd_controls$epatid) |> arrange(epatid)


#cohort
cohort <- readRDS(paste0(path_read, "e_17_205r_cohort_final_r.rds"))
cohort <- cohort[cohort$e_patid %in% patid, ]
saveRDS(cohort, paste0(path_save,"e_17_205r_cohort_final_matched.rds"))
rm(cohort)

# CPRD
print("CPRD patients..")
cprd_patients <- readRDS(paste0(path_read, "e_17_205r_extract_patient_final_r.rds"))
cprd_patients <- cprd_patients[cprd_patients$epatid %in% patid, ]
saveRDS(cprd_patients, paste0(path_save,"patient_final_cohort.rds"))
rm(cprd_patients)

##clinical
print("CPRD clinical..")
clinical1 <- readRDS(paste0(path_read, "e_17_205r_extract_clinical_final_1_50_r.rds"))
clinical1 <- clinical1[clinical1$e_patid %in% patid, ]

clinical2 <- readRDS(paste0(path_read, "e_17_205r_extract_clinical_final_51_81_r.rds"))
clinical2 <- clinical2[clinical2$e_patid %in% patid,]

clinical <- rbind(clinical1, clinical2)
saveRDS(clinical, paste0(path_save,"clinical_final_cohort.rds"))
rm(clinical1, clinical2, clinical)

##additional
print("CPRD additional..")
additional <- readRDS(paste0(path_read, "e_17_205r_extract_additional_final_r.rds"))
additional <- additional[additional$e_patid %in% patid,] 
saveRDS(additional, paste0(path_save,"additional_final_cohort.rds"))
rm(additional)
gc()

##referral
print("CPRD referral..")
referral <- readRDS(paste0(path_read, "e_17_205r_extract_referral_final_r.rds"))
referral <- referral[referral$e_patid %in% patid,]
saveRDS(referral, paste0(path_save,"referral_final_cohort.rds"))
rm(referral)
gc()

##test
print("CPRD test..")
test1 <- readRDS(paste0(path_read, "e_17_205r_extract_test_final_1_35_r.rds"))
test1 <- test1[test1$e_patid %in% patid,]

test2 <- readRDS(paste0(path_read, "e_17_205r_extract_test_final_36_60_r.rds"))
test2 <- test2[test2$e_patid %in% patid,]

test3 <- readRDS(paste0(path_read, "e_17_205r_extract_test_final_61_77_r.rds"))
test3 <- test3[test3$e_patid %in% patid,]

test <- rbind(test1, test2, test3)
saveRDS(test, paste0(path_save,"test_final_cohort.rds"))
rm(test1, test2, test3, test)
gc()

##therapy
print("CPRD therapy..")
therapy1 <- readRDS(paste0(path_read, "e_17_205r_extract_therapy_final_1_40_r.rds"))
therapy1 <- therapy1[therapy1$e_patid %in% patid,]

therapy2 <- readRDS(paste0(path_read, "e_17_205r_extract_therapy_final_41_80_r.rds"))
therapy2 <- therapy2[therapy2$e_patid %in% patid,]

therapy3 <- readRDS(paste0(path_read, "e_17_205r_extract_therapy_final_81_120_r.rds"))
therapy3 <- therapy3[therapy3$e_patid %in% patid,]

therapy4 <- readRDS(paste0(path_read, "e_17_205r_extract_therapy_final_121_160_r.rds"))
therapy4 <- therapy4[therapy4$e_patid %in% patid,]

therapy5 <- readRDS(paste0(path_read, "e_17_205r_extract_therapy_final_161_189_r.rds"))
therapy5 <- therapy5[therapy5$e_patid %in% patid,]

therapy <- rbind(therapy1, therapy2, therapy3, therapy4, therapy5)
saveRDS(therapy, paste0(path_save,"therapy_final_cohort.rds"))
rm(therapy1, therapy2, therapy3, therapy4, therapy5, therapy)
gc()

write(paste(Sys.time(), "Done CPRD 'cleaning"), file = progress_file, append = T)

# HES
print("HES..")
hes_patients <- readRDS(paste0(path_read, "e_hes_patient_17_205r_final.rds"))
hes_patients <- hes_patients[hes_patients$e_patid %in% patid,]
saveRDS(hes_patients, paste0(path_save,"hes_patient_final_cohort.rds"))
rm(hes_patients)
gc()

hosp <- readRDS(paste0(path_read, "e_hes_hospital_17_205r_final.rds"))
hosp <- hosp[hosp$e_patid %in% patid,]
saveRDS(hosp, paste0(path_save,"hes_hospital_final_cohort.rds"))
rm(hosp)
gc()

epi <- readRDS(paste0(path_read, "e_hes_episodes_17_205r_final.rds"))
epi <- epi[epi$e_patid %in% patid,]
saveRDS(epi, paste0(path_save,"hes_episodes_final_cohort.rds"))
rm(epi)
gc()

diag <- readRDS(paste0(path_read, "e_hes_diagnosis_epi_17_205r_final.rds"))
diag <- diag[diag$e_patid %in% patid,]
saveRDS(diag, paste0(path_save,"hes_diagnosis_epi_final_cohort.rds"))
rm(diag)
gc()

proc <- readRDS(paste0(path_read, "e_hes_procedures_epi_17_205r_final.rds"))
proc <- proc[proc$e_patid %in% patid,]
saveRDS(proc, paste0(path_save,"hes_procedures_epi_final_cohort.rds"))
rm(proc)
gc()

acp <- readRDS(paste0(path_read, "e_hes_acp_17_205r_final.rds"))
acp <- acp[acp$e_patid %in% patid,]
saveRDS(acp, paste0(path_save,"hes_acp_final_cohort.rds"))
rm(acp)
gc()

cc <- readRDS(paste0(path_read, "e_hes_ccare_17_205r_final.rds"))
cc <- cc[cc$e_patid %in% patid,]
saveRDS(cc, paste0(path_save,"hes_ccare_final_cohort.rds"))
rm(cc)
gc()

write(paste(Sys.time(), "Done HES 'cleaning"), file = progress_file, append = T)

# Cancer
print("cancer..")
cancer <- readRDS(paste0(path_read, "e_17_205_cancer_registration_1_final.rds"))
cancer <- cancer[cancer$e_patid %in% patid,]
saveRDS(cancer, paste0(path_save,"cancer_registration_final_cohort.rds"))
rm(cancer)
gc()

write(paste(Sys.time(), "Done Cancer 'cleaning"), file = progress_file, append = T)

# deaths
print("deaths..")
death <- readRDS(paste0(path_read, "e_death_patient_17_205r_final.rds"))
death <- death[death$e_patid %in% patid,]
saveRDS(death, paste0(path_save,"death_patient_final_cohort.rds"))
rm(death)
gc()

write(paste(Sys.time(), "Done ONS death 'cleaning"), file = progress_file, append = T)

# imd
print("imd..")
imd <- readRDS(paste0(path_read, "e_practice_imd_17_205r_final.rds"))
imd <- imd[imd$e_practid %in%pracid,]
saveRDS(imd, paste0(path_save,"practice_imd_final_cohort.rds"))
rm(imd)
gc()

write(paste(Sys.time(), "Done IMD 'cleaning"), file = progress_file, append = T)

write(paste(Sys.time(), "End 6_clean_files.R script"), file = progress_file, append = T)


rm(list=ls()); gc()
source("./scripts/7_covariates.R")


