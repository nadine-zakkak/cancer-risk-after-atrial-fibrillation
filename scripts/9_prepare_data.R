rm(list=ls()); gc()
# Need to update comments 23/June

# Need to check the proportional hazards assumption
library(plyr)
library(dplyr)
library(lubridate)

#load functions, global variables and codelists ----
source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
source("./scripts/00_read_codelists.R")

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 9_prepare_data_v1.R script"), file = progress_file, append = T)

# global var  ------
allcov <- c("case", "age_start", "diab", "ht", "smok", "alc_desc")
common_cols <- c("epatid", "gender", "follow_3m", "follow_5y")
save_data_dir <- "./data/final_data/"
suffix_short <- "_short_term"
suffix_med <- "_med_term"
suffix_long <- "_long_term"

# load data ----
alc_dict <- read.csv("./look_ups/alcohol_dictionary.csv")

flag_followup_time <- function(data){
  # define cut-offs for follow-up times
  # follow_xx: followed up until xx time?
  # short, med, long: wich analysis should they be included in?
  data |>
    mutate(
      follow_3m =  follow <= 3,
      follow_5y = follow/12 <= 5,
      short = T, 
      med = !follow_3m,
      long = !follow_3m & !follow_5y
    )
}


calc_followup_survival <- function(data){
  # Calculate follow up time for survival analysis
  # Medium-term (3m-5y): shift start date by 3 months and calculate follow-up time (until outcome)
  # if the 'original' follow-up continues beyond 5y then limit medium-term follow-up to 4.75 (4y 9m) because start date has been shifted by 3m
  # Long-term (Beyond 5y): shift start date by 5 years and calculate follow-up time 
  # 4 new columns:
  # start_med: "new" start of follow-up for medium-term analysis 
  # follow_sa_med: follow-up time (new start until outcome) for medium-term analysis 
  # and their equivalent for long-term analysis: start_long, follow_sa_long
  data |> 
    mutate(follow_short = ifelse(!follow_3m, 3, follow), #if followed for more than 3m, then limit to 3 months (0.25 years)
           start_med = if_else(med, start %m+% months(3), NA), #shift start date by 3 months if it is included in medium 
           follow_sa_med =   time_length(start_med %--% outcome_date, unit = "month"),
           follow_sa_med = ifelse(follow>5*12, 4.75*12, follow_sa_med),
           start_long = start %m+% years(5), #shift 5y
           follow_sa_long = time_length(start_long %--% outcome_date, unit = "month"))
}

# 
identify_period_specific_outcomes <- function(data) {
  # Identify the outcomes for the medium and long-term analysis
  # Medium-term: 
  ## if it is to be included but the follow-up continues beyond 5y, then censor at 5y
  # 4 new columns:
  # outcome_med: includes the outcome for the medium-term analysis as described above
  # outcome_num_med: outcome_med in numeric format (calls external function in functions.R)
  # & their equivalent for long-term follow-up: outcome_long & outcome_num_long
  data |> 
    mutate(outcome_short = ifelse(short & !follow_3m, "end", outcome_mod),
           outcome_num_short = num_outcomes(outcome_short),
           outcome_med = ifelse(med & !follow_5y, "end", outcome_mod),
           outcome_num_med = num_outcomes(outcome_med),
           outcome_long = outcome_mod,
           outcome_num_long = num_outcomes(outcome_long)
    )
}


# For all the loaded datasets,
# 1) calculate follow-up time for medium-term (3m-5y) and long-term analysis (Beyond 5y)
# 2) Identify the outcomes specific to each analysis (i.e. censor at 5y for medium-term analysis)
# 3) Select columns relevant to future analysis
# 4) Save 'updated' dataset

# Cancer Death dataset ------
# load data
cancer_death <- readRDS("./data/final_data/cancer_death.RDS")
endpoint <- "cancer_death"

# Flag follow-up time limits and specify into short, med or long follow-up
cancer_death_analysis_data <- cancer_death |>
  mutate(alc_desc = factor(add_alc_desc(alc, alc_dict), levels = alc_levels),
         outcome_mod = ifelse(outcome == "af", "end", outcome), # mark af as censored
         outcome_mod_numeric = num_outcomes(outcome_mod)) |>
  flag_followup_time() |>
  calc_followup_survival() |>
  identify_period_specific_outcomes()

# Short-term analysis dataset
cancer_death_short <-
  cancer_death_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), short, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         follow = follow_short, outcome = outcome_short, outcome_num = outcome_num_short) |>
  mutate(cancer_death_outcome =  outcome=="cancer death") |>
  filter(short)

#check data
summary(cancer_death_short$follow_orig)
summary(cancer_death_short$follow)

saveRDS(object = cancer_death_short, file = paste0(save_data_dir, endpoint, suffix_short, ".RDS"))

# Medium-term analysis dataset
cancer_death_med <-
  cancer_death_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), med, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_med,
         follow = follow_sa_med, outcome = outcome_med, outcome_num = outcome_num_med) |>
  mutate(cancer_death_outcome =  outcome =="cancer death") |>
  filter(med)

#check data
summary(cancer_death_med$follow_orig)
summary(cancer_death_med$follow)

saveRDS(object = cancer_death_med, file= paste0(save_data_dir, endpoint, suffix_med, ".RDS"))

# Long-term analysis dataset
cancer_death_long <-
  cancer_death_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), long, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_long,
         follow = follow_sa_long, outcome = outcome_long, outcome_num = outcome_num_long) |>
  mutate(cancer_death_outcome =  outcome =="cancer death") |>
  filter(long)

summary(cancer_death_long$follow_orig)
summary(cancer_death_long$follow)

saveRDS(object = cancer_death_long, file =  paste0(save_data_dir, endpoint, suffix_long, ".RDS"))

# clean-up
rm(list=ls(pattern="*death*"), endpoint)

# Cancer Incidence dataset ------
# load data 
cancer_diag <- readRDS("./data/final_data/cancer_diag.RDS")
endpoint <- "cancer_diag"

cancer_diag_analysis_data <- cancer_diag |>
  mutate(alc_desc = factor(add_alc_desc(alc, alc_dict), levels = alc_levels),
         outcome_mod = ifelse(outcome == "af", "end", outcome), # mark af as censored
         outcome_mod_numeric = num_outcomes(outcome_mod)) |>
  flag_followup_time() |>
  calc_followup_survival() |>
  identify_period_specific_outcomes()

cancer_diag_short <-
  cancer_diag_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), short, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         follow = follow_short, outcome = outcome_short, outcome_num = outcome_num_short) |>
  filter(short) |>
  mutate(cancer_outcome =  outcome=="cancer")


#check data
summary(cancer_diag_short$follow_orig)
summary(cancer_diag_short$follow)


saveRDS(object = cancer_diag_short, file= paste0(save_data_dir, endpoint, suffix_short, ".RDS"))

cancer_diag_med <-
  cancer_diag_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), med, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_med,
         follow = follow_sa_med, outcome = outcome_med, outcome_num = outcome_num_med) |>
  filter(med) |>
  mutate(cancer_outcome =  outcome =="cancer")


#check data
summary(cancer_diag_med$follow_orig)
summary(cancer_diag_med$follow)

saveRDS(object = cancer_diag_med, file= paste0(save_data_dir, endpoint, suffix_med, ".RDS"))


cancer_diag_long <-
  cancer_diag_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), long, 
         follow_orig = follow, start_orig = start, outcome_date,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_long,
         follow = follow_sa_long, outcome = outcome_long, outcome_num = outcome_num_long) |>
  filter(long) |>
  mutate(cancer_outcome =  outcome =="cancer")

#check data
summary(cancer_diag_long$follow_orig)
summary(cancer_diag_long$follow)

saveRDS(object = cancer_diag_long, file= paste0(save_data_dir, endpoint, suffix_long, ".RDS"))

# clean up
rm(list=ls(pattern="*diag"), endpoint)

# Cancer incidence by region dataset -----
# load data
cancer_region <- readRDS("./data/final_data/cancer_region.RDS")
endpoint <- "cancer_region"

cancer_region_analysis_data <-
  cancer_region |>
  add_region_desc() |>
  mutate(alc_desc = factor(add_alc_desc(alc, alc_dict), levels = alc_levels),
         outcome_mod = ifelse(outcome == "af", "end", outcome), # mark af as censored
         outcome_mod_numeric = num_outcomes(outcome_mod),
         excl_analysis = cancer_bodyregion %in% c(5,6)) |> #Exclude if Brain/CNS or Upper and lower limbs)
  flag_followup_time() |>
  calc_followup_survival() |>
  identify_period_specific_outcomes()


cancer_region_short <-
  cancer_region_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), short, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_bodyregion, region_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         follow = follow_short, outcome = outcome_short, outcome_num = outcome_num_short) |>
  filter(short) |>
  mutate(cancer_outcome =  outcome=="cancer")

summary(cancer_region_short$follow_orig)
summary(cancer_region_short$follow)

saveRDS(object = cancer_region_short, file= paste0(save_data_dir, endpoint, suffix_short, ".RDS"))

cancer_region_med <-
  cancer_region_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), med, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_bodyregion, region_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_med,
         follow = follow_sa_med, outcome = outcome_med, outcome_num = outcome_num_med) |>
  filter(med) |>
  mutate(cancer_outcome =  outcome=="cancer")

#check data
summary(cancer_region_med$follow_orig)
summary(cancer_region_med$follow)

saveRDS(object = cancer_region_med, file= paste0(save_data_dir, endpoint, suffix_med, ".RDS"))

cancer_region_long <-
  cancer_region_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), long, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_bodyregion, region_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_long,
         follow = follow_sa_long, outcome = outcome_long, outcome_num = outcome_num_long) |>
  filter(long) |>
  mutate(cancer_outcome =  outcome=="cancer")

#check data
summary(cancer_region_long$follow_orig)
summary(cancer_region_long$follow)

saveRDS(object = cancer_region_long, file= paste0(save_data_dir, endpoint, suffix_long, ".RDS"))

# clean up
rm(list=ls(pattern="*region"), endpoint)

# Cancer incidence by organ system dataset ------
# load data
cancer_organsys <- readRDS("./data/final_data/cancer_system.RDS")
endpoint <- "cancer_organsys"

cancer_organsys_analysis_data <-
  cancer_organsys |>
  add_organsys_desc() |>
  mutate(alc_desc = factor(add_alc_desc(alc, alc_dict), levels = alc_levels),
         outcome_mod = ifelse(outcome == "af", "end", outcome), # mark af as censored
         outcome_mod_numeric = num_outcomes(outcome_mod),
         excl_analysis = (cancer_organsys == 10 & gender == 0) | #Breast cancer in men
           cancer_organsys %in% c(1, 2, 3, 4)) |> #Cardiovascular, Endocrine, Nervous, Skeletal
  flag_followup_time() |>
  calc_followup_survival() |>
  identify_period_specific_outcomes()


cancer_organsys_short <-
  cancer_organsys_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), short, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_organsys, organsys_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         follow = follow_short, outcome = outcome_short, outcome_num = outcome_num_short) |>
  filter(short) |>
  mutate(cancer_outcome =  outcome=="cancer")

summary(cancer_organsys_short$follow_orig)
summary(cancer_organsys_short$follow)

saveRDS(object = cancer_organsys_short, file= paste0(save_data_dir, endpoint, suffix_short, ".RDS"))

cancer_organsys_med <-
  cancer_organsys_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), med, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_organsys, organsys_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_med,
         follow = follow_sa_med, outcome = outcome_med, outcome_num = outcome_num_med) |>
  filter(med) |>
  mutate(cancer_outcome =  outcome=="cancer")

#check data
summary(cancer_organsys_med$follow_orig)
summary(cancer_organsys_med$follow)

saveRDS(object = cancer_organsys_med, file= paste0(save_data_dir, endpoint, suffix_med, ".RDS"))

cancer_organsys_long <-
  cancer_organsys_analysis_data |>
  select(!!!syms(common_cols), !!!syms(allcov), long, 
         follow_orig = follow, start_orig = start, outcome_date,
         cancer_organsys, organsys_desc, excl_analysis,
         outcome_orig = outcome_mod, outcome_num_orig = outcome_mod_numeric,
         start_long,
         follow = follow_sa_long, outcome = outcome_long, outcome_num = outcome_num_long) |>
  filter(long) |>
  mutate(cancer_outcome =  outcome=="cancer")

#check data
summary(cancer_organsys_long$follow_orig)
summary(cancer_organsys_long$follow)

saveRDS(object = cancer_organsys_long, file= paste0(save_data_dir, endpoint, suffix_long, ".RDS"))

# clean up
rm(list=ls(pattern="*organsys"), endpoint)

write(paste(Sys.time(), "End 9_prepare_data_V1.R script"), file = progress_file, append = T)




