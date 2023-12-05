library("readr")

read_files <- function(file_name, begin, end, types){
  #Function to read files in a sequence
  #Takes as input the file name, beginning of the seq and end of the sequence
  #Returns a dataframe with all the files contents
  
  sequence <- seq(begin, end, 1)
  file_name <- file_name
  file_names <- sprintf(file_name, sequence)
  df <- data.frame()
  for (file in file_names){
    print(file)
    df <- rbind(df, read_delim(file, delim="\t",
                               col_types = types))
  }
  
  return(df)
}

create_name <- function(filename, from, to){
  return(paste0(filename,"_", from, "_", to, "_r.rds"))
}

#CPRD
#Additional CPRD ----
additional_types <- cols_only(e_patid = col_integer(),
                              enttype = col_integer(),
                              adid = col_integer(), 
                              #all data# could be date or integer 
                              #depending on enttype --> listed as character for now
                              data1 = col_character(), 
                              data2 = col_character(),
                              data3 = col_character(),
                              data4 = col_character(),
                              data5 = col_character(), 
                              data6 = col_character(),
                              data7 = col_character())
additional <- read_files("e_17_205r_extract_additional_%03d_r.txt", 1, 12, 
                         additional_types)
saveRDS(additional, "../../R_Data/e_17_205r_extract_additional_final_r.rds")

rm(additional_types)
gc()

#Clinical CPRD ----
filename="../../R_Data/e_17_205r_extract_clinical_final"
clinical_types <- cols_only(e_patid = col_integer(),
                            eventdate = col_date(format="%d/%m/%Y"),
                            medcode = col_integer(),
                            episode = col_integer(),
                            enttype = col_integer(),
                            adid = col_integer())

#saved: 1-50 (as rds)
from=1
to=50
clinical <- read_files("e_17_205r_extract_clinical_%03d_r.txt",
                        from, to, clinical_types)
saveRDS(clinical, create_name(filename, from, to))
rm(clinical)
gc()
 
#saved: 51-81 (as rds)
from=51
to=81
clinical1 <- read_files("e_17_205r_extract_clinical_%03d_r.txt",
                        from, to, clinical_types)
saveRDS(clinical1, create_name(filename, from, to))
rm(clinical1)
gc()

rm(clinical_types)

#Referral CPRD ---- 
referral_types <- cols_only(e_patid = col_integer(),
                            eventdate = col_date(format="%d/%m/%Y"),
                            medcode = col_integer(),
                            source = col_integer(),
                            inpatient = col_integer(),
                            attendance = col_integer())
referral <- read_files("e_17_205r_extract_referral_%03d_r.txt",
                       1, 4, referral_types)
saveRDS(referral, "../../R_Data/e_17_205r_extract_referral_final_r.rds")

rm(referral)
rm(referral_types)
gc()

#Test CPRD -----
filename = "../../R_Data/e_17_205r_extract_test_final"
test_types <- cols_only(e_patid = col_integer(),
                        eventdate = col_date(format="%d/%m/%Y"),
                        medcode = col_integer(),
                        data1 = col_integer(), #int for all enttypes
                        data2 = col_double(), #double for all
                        data3 = col_double(), #double or int
                        data4 = col_integer(), #double or int
                        data5 = col_double(), #double or non-existent
                        data6 = col_double(), #double or non-existent
                        data7 = col_integer(), #int or non-existent
                        data8 = col_character()) #date or int or non-existent


#saved: 1-35 (as rds)
from=1
to=35
test1 <- read_files("e_17_205r_extract_test_%03d_r.txt",
                    from, to, test_types) 
saveRDS(test1, create_name(filename, from, to))
rm(test1)
gc()

#saved:36-60 (as rds)
from=36
to=60
test2 <- read_files("e_17_205r_extract_test_%03d_r.txt",
                    from, to, test_types)
saveRDS(test2, create_name(filename, from, to))
rm(test2)
gc()

#saved: 61-77 (as rds)
from=61
to=77
test3 <- read_files("e_17_205r_extract_test_%03d_r.txt",
                    from, to, test_types)
saveRDS(test3, create_name(filename, from, to))
rm(test3)
rm(test_types)
gc()

#Therapy CPRD -----
filename="../../R_Data/e_17_205r_extract_therapy_final"
therapy_types <- cols_only(e_patid = col_integer(),
                           eventdate = col_date(format="%d/%m/%Y"),
                           consid = col_integer(),
                           prodcode = col_integer(),
                           dosageid = col_character(),
                           bnfcode = col_integer(),
                           qty = col_double(), #int in doc
                           numdays = col_integer())

#saved: 1-40 (as rds)
from=1
to=40
therapy1 <- read_files("e_17_205r_extract_therapy_%03d_r.txt", 
                       from, to, therapy_types)
saveRDS(therapy1, create_name(filename, from, to))
rm(therapy1)
gc()

#saved: 41-80 (as rds)
from=41
to=80
therapy2 <- read_files("e_17_205r_extract_therapy_%03d_r.txt", 
                       from, to, therapy_types)
saveRDS(therapy2, create_name(filename, from, to))
rm(therapy2)
gc()

#saved: 81-120 (as rds)
from=81
to=120
therapy3 <- read_files("e_17_205r_extract_therapy_%03d_r.txt", 
                       from, to, therapy_types)
saveRDS(therapy3, create_name(filename, from, to))
rm(therapy3)
gc()

#saved: 121-160 (as rds)
from=121
to=160
therapy4 <- read_files("e_17_205r_extract_therapy_%03d_r.txt", 
                       from, to, therapy_types)
saveRDS(therapy4, create_name(filename, from, to))
rm(therapy4)
gc()

#saved: 161-189 (as rds)
from=161
to=189
therapy5 <- read_files("e_17_205r_extract_therapy_%03d_r.txt", 
                       from, to, therapy_types)
saveRDS(therapy5, create_name(filename, from, to))
rm(therapy5)
gc()

rm(therapy_types)

#Patient CPRD ----
patient_types <- cols_only(epatid = col_integer(),
                           gender = col_integer(),
                           yob = col_integer(),
                           frd = col_date(format="%d/%m/%Y"),
                           crd = col_date(format="%d/%m/%Y"),
                           tod = col_date(format='%d/%m/%Y'),
                           deathdate = col_date(format="%d/%m/%Y"),
                           accept = col_integer())
patient <- read_files("e_17_205r_Extract_Patient_%03d_r.txt", 
                      1, 1, patient_types)
saveRDS(patient, "../../R_Data/e_17_205r_extract_patient_final_r.rds")

#Practice CPRD ----
prac_types <- cols_only(epracid = col_integer(),
                        region = col_integer(),
                        lcd = col_date(format="%d/%m/%Y") ,
                        uts = col_date(format="%d/%m/%Y"))
practice <- read_files("e_17_205r_Extract_Practice_%03d_r.txt",
                       1, 1, prac_types)
saveRDS(practice, "../../R_Data/e_17_205r_extract_practice_final_r.rds")


#Follow-up CPRD -----
follow_types <- cols_only(e_patid = col_integer(),
                          start = col_date(format="%d/%m/%Y"),
                          end = col_date(format="%d/%m/%Y"))
print("e_17_205r_cohort_file.txt")
follow <- read_delim("e_17_205r_cohort_file.txt", delim="\t",
                     col_types = follow_types)
saveRDS(follow, "../../R_Data/e_17_205r_cohort_final_r.rds")


#Linked Data ----
#HES ----
##Patients ----
patient_types <- cols_only(e_patid = col_integer(),
                           e_practid = col_integer(),
                           e_gen_hesid = col_integer(),
                           n_patid_hes = col_integer(),
                           gen_ethnicity = col_character(),
                           match_rank = col_integer())

patient <- read_delim("e_hes_patient_17_205r.txt", delim="\t",
                     col_types = patient_types)
saveRDS(patient, "../../R_Data/e_hes_patient_17_205r_final.rds")

rm(patient, patient_types)
gc()

##Hospitalisations -----
hosp_types <- cols_only(e_patid = col_integer(),
                        spno = col_integer(),
                        admidate = col_date(format="%d/%m/%Y"),
                        discharged = col_date(format="%d/%m/%Y"),
                        admimeth = col_character(),
                        admisorc = col_integer(),
                        disdest = col_integer(),
                        dismeth = col_integer(),
                        duration = col_integer())

hosp <- read_delim("e_hes_hospital_17_205r.txt", delim="\t",
                   col_types = hosp_types)

saveRDS(hosp, "../../R_Data/e_hes_hospital_17_205r_final.rds")

rm(hosp, hosp_types)
gc()

##Episodes -----
epi_types <- cols_only( e_patid = col_integer(),
                        spno = col_integer(),
                        epikey = col_integer(),
                        admidate = col_date(format="%d/%m/%Y"),
                        discharged = col_date(format="%d/%m/%Y"),
                        eorder = col_integer(),
                        disdest = col_integer(),
                        dismeth = col_integer(),
                        mainspef = col_character(),
                        tretspef = col_character(),
                        classpat = col_integer(),
                        ethnos = col_character())

epi <- read_delim("e_hes_episodes_17_205r.txt", delim="\t",
                   col_types = epi_types)

saveRDS(epi, "../../R_Data/e_hes_episodes_17_205r_final.rds")

rm(epi, epi_types)
gc()

##Diag -----
diag_types <- cols_only(e_patid = col_integer(),
                        spno = col_integer(),
                        epikey = col_integer(),
                        epistart = col_date(format="%d/%m/%Y"),
                        epiend = col_date(format="%d/%m/%Y"),
                        icd = col_character(),
                        icdx = col_character(),
                        d_order = col_integer())

diag <- read_delim("e_hes_diagnosis_epi_17_205r.txt", delim="\t",
                   col_types = diag_types)

saveRDS(diag, "../../R_Data/e_hes_diagnosis_epi_17_205r_final.rds")

rm(diag, diag_types)
gc()

##procedures -----
proc_types = cols_only(e_patid = col_integer(),
                       spno = col_integer(),
                       epikey = col_integer(),
                       admidate = col_date(format="%d/%m/%Y"),
                       epistart = col_date(format="%d/%m/%Y"),
                       epiend = col_date(format="%d/%m/%Y"),
                       discharged = col_date(format="%d/%m/%Y"),
                       opcs = col_character(),
                       evdate = col_date(format="%d/%m/%Y"),
                       p_order = col_integer())

proc <- read_delim("e_hes_procedures_epi_17_205r.txt", delim="\t",
                   col_types = proc_types)

saveRDS(proc, "../../R_Data/e_hes_procedures_epi_17_205r_final.rds")
rm(proc, proc_types)
gc()

##Augmented Care Period (acp) -----
acp_types = cols_only(e_patid = col_integer(),
                      spno = col_integer(),
                      epikey = col_integer(),
                      epistart = col_date(format="%d/%m/%Y"),
                      epiend = col_date(format="%d/%m/%Y"),
                      eorder = col_integer(),
                      epidur = col_integer(),
                      numacp = col_integer(),
                      acpn = col_integer(),
                      acpstar = col_date(format="%d/%m/%Y"),
                      acpend = col_date(format="%d/%m/%Y"),
                      acploc = col_integer(),
                      acpout = col_integer(),
                      acpplan = col_character(),
                      acpspef = col_character(),
                      orgsup = col_integer())

acp <- read_delim("e_hes_acp_17_205r.txt", delim="\t",
                   col_types = acp_types)

saveRDS(acp, "../../R_Data/e_hes_acp_17_205r_final.rds")
rm(acp, acp_types)
gc()

##Critical Care (cc) ------
cc_types = cols_only( e_patid = col_integer(),
                       spno = col_integer(),
                       epikey = col_integer(),
                       admidate = col_date(format="%d/%m/%Y"),
                       discharged = col_date(format="%d/%m/%Y"),
                       epistart = col_date(format="%d/%m/%Y"),
                       epiend = col_date(format="%d/%m/%Y"),
                       eorder = col_integer(),
                       ccstartdate = col_date(format="%d/%m/%Y"),
                       ccstarttime = col_character(),
                       ccdisdate = col_date(format="%d/%m/%Y"),
                       ccdistime = col_character(),
                      ccadmitype = col_integer(),
                      ccadmisorc = col_integer(),
                      ccdisstat = col_integer(),
                      ccdisdest = col_integer(),
                      cclev2days = col_integer(),
                      cclev3days = col_integer(),
                      bcardsupdays = col_integer(),
                      acardsupdays = col_integer(),
                      bressupdays = col_integer(),
                      aressupdays = col_integer(),
                      gisupdays = col_integer(),
                      liversupdays = col_integer(),
                      neurosupdays = col_integer(),
                      rensupdays = col_integer(),
                      dermsupdays = col_integer(),
                      orgsupmax = col_integer(),
                      ccunitfun = col_integer(),
                      unitbedconfig = col_integer())

cc <- read_delim("e_hes_ccare_17_205r.txt", delim="\t",
                  col_types = cc_types)

saveRDS(cc, "../../R_Data/e_hes_ccare_17_205r_final.rds")
rm(cc, cc_types)
gc()

#Cancer ------
cancer_types <- cols(
  .default = col_character(),
  e_patid = col_integer(), #number
  e_cr_patid = col_integer(), #number
  e_cr_id = col_integer(), #number
  age = col_double(), #number
  diagnosisdatebest = col_date(format="%d/%m/%Y"),
  diagnosisdateflag = col_integer(), #number
  basisdiagnosis = col_integer(), #some fields have "&-->NA". number
  behaviour_coded = col_integer(), #number
  coding_system = col_integer(), #number
  tumoursize = col_double(), #number
  stage_best_system = col_integer() #number
)

cancer <- read_delim("e_17_205_cancer_registration_1.txt", 
                     delim="\t",
                     col_types = cancer_types)

saveRDS(cancer, "../../R_Data/e_17_205_cancer_registration_1_final.rds")
rm(cancer, cancer_types)
gc()

#ONS - death -----
death_types <- cols(
  .default = col_character(),
  e_patid = col_integer(),
  e_practid = col_integer(),
  e_gen_death_id = col_integer(),
  n_patid_death = col_integer(),
  match_rank = col_integer(),
  dor = col_date(format='%d/%m/%Y'),
  dod = col_date(format='%d/%m/%Y')
)
death <- read_delim("e_death_patient_17_205r.txt", delim="\t",
                    col_types = death_types)
saveRDS(death, "../../R_Data/e_death_patient_17_205r_final.rds")
rm(death, death_types)
gc()

#imd -----
imd_types <- cols(.default=col_integer(),
                  country = col_character())
imd <- read_delim("e_practice_imd_17_205r.txt", delim="\t",
                  col_types = imd_types)
saveRDS(imd, "../../R_Data/e_practice_imd_17_205r_final.rds")
rm(imd, imd_types)
gc()
