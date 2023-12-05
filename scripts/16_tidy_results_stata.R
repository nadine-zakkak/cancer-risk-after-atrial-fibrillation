# tidy results tables from stata to use for table and figures outputs in R

library(dplyr)

source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
source("./scripts/00_read_codelists.R")

edit_stata_results <- function(data, gender){
  data |>
    filter(var == "1.case") |>
    mutate(gender = gender, term = "case") |>
    rename(estimate = b, lower = ll, upper = ul, p = pvalue)
}

path_read <- "./results/stata/"
path_save <- "./results/final_results_med_long/"

# Cancer Death -----
med_male <- read.csv(paste0(path_read, "cancer_death_med_term_male1.csv")) 
med_female <- read.csv(paste0(path_read, "cancer_death_med_term_male0.csv")) 
med <- bind_rows(edit_stata_results(med_male, "Men"), edit_stata_results(med_female, "Women")) |>
  mutate(outcome = "cancer death", fup = "med")

write.csv(med, paste0(path_save, "med_cancer_death.csv"), row.names = F)

long_male <- read.csv(paste0(path_read, "cancer_death_long_term_male1.csv")) 
long_female <- read.csv(paste0(path_read, "cancer_death_long_term_male0.csv")) 
long <- bind_rows(edit_stata_results(long_male, "Men"), edit_stata_results(long_female, "Women")) |>
  mutate(outcome = "cancer death", fup = "long")

write.csv(long, paste0(path_save, "long_cancer_death.csv"), row.names = F)

rm(list=ls(pattern = "*(med|long)"))

# Cancer Diagnosis -----
med_male <- read.csv(paste0(path_read, "cancer_diag_med_term_male1.csv")) 
med_female <- read.csv(paste0(path_read, "cancer_diag_med_term_male0.csv")) 
med <- bind_rows(edit_stata_results(med_male, "Men"), edit_stata_results(med_female, "Women"))|>
  mutate(outcome = "cancer incidence", fup = "med")

write.csv(med, paste0(path_save, "med_cancer.csv"), row.names = F)

long_male <- read.csv(paste0(path_read, "cancer_diag_long_term_male1.csv")) 
long_female <- read.csv(paste0(path_read, "cancer_diag_long_term_male0.csv")) 
long <- bind_rows(edit_stata_results(long_male, "Men"), edit_stata_results(long_female, "Women")) |>
  mutate(outcome = "cancer incidence", fup = "long")

write.csv(long, paste0(path_save, "long_cancer.csv"), row.names = F)

rm(list=ls(pattern = "*(med|long)"))

# Cancer by body region -----
med <- data.frame()
long <- data.frame()
for(male in c(1, 0)){
  df <- do.call(rbind, lapply(list.files(path=path_read, pattern=paste0("^cancer_region_med_term_male", male)), 
                              FUN = function(x) read.csv(paste0(path_read, x))))
  df1 <- do.call(rbind, lapply(list.files(path=path_read, pattern=paste0("^cancer_region_long_term_male", male)), 
                               FUN = function(x) read.csv(paste0(path_read, x))))
  sex <- ifelse(male == 1, "Men", "Women")
  med <- bind_rows(med, edit_stata_results(df, sex) |> 
                     rename(cancer_bodyregion = subtype) |>
                     add_region_desc() |>
                     mutate(outcome = "cancer region", fup = "med")) 
  long <- bind_rows(long, edit_stata_results(df1, sex) |> 
                      rename(cancer_bodyregion = subtype) |>
                      add_region_desc() |>
                      mutate(outcome = "cancer region", fup = "long")) 
  
}

write.csv(med, paste0(path_save, "med_region.csv"), row.names = F)
write.csv(long, paste0(path_save, "long_region.csv"), row.names = F)

rm(list=ls(pattern = "*(med|long)"))

# Cancer by organ system -----
med <- data.frame()
long <- data.frame()
for(male in c(1, 0)){
  df <- do.call(rbind, lapply(list.files(path=path_read, pattern=paste0("^cancer_organsys_med_term_male", male)), 
                              FUN = function(x) read.csv(paste0(path_read, x))))
  df1 <- do.call(rbind, lapply(list.files(path=path_read, pattern=paste0("^cancer_organsys_long_term_male", male)),
                               FUN = function(x) read.csv(paste0(path_read, x))))
  sex <- ifelse(male == 1, "Men", "Women")
  med <- bind_rows(med, edit_stata_results(df, sex) |> 
                     rename(cancer_organsys = subtype) |>
                     add_organsys_desc() |>
                     mutate(outcome = "cancer organ sys", fup = "med")) 
  long <- bind_rows(long, edit_stata_results(df1, sex) |>
                      rename(cancer_organsys = subtype) |>
                      add_organsys_desc()  |>
                      mutate(outcome = "cancer organ sys", fup = "long"))
  
}

write.csv(med, paste0(path_save, "med_organsys.csv"), row.names = F)
write.csv(long, paste0(path_save, "long_organsys.csv"), row.names = F)

rm(list=ls(pattern = "*(med|long)"))



