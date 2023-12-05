# Draw cumulative incidence curves starting after 3 months
rm(list=ls());gc()
library(plyr)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(forcats)
library(survival)
library(tidycmprsk)
library(ggplot2)

# load functions, global variables and codelists
source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
source("./scripts/00_read_codelists.R")

save_image <- list(device = "jpeg", unit = "mm", width = 140, height = 80, dpi=1000)

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 11_cumulative_incidence_descriptive.R script"), file = progress_file, append = T)

# function
shift_follow_up_3m <- function(data){
  data |>
    mutate(follow_sa = time_length(start_orig %m+% months(3) %--% outcome_date, unit = "month"))
}

change_outcome_to_factor <- function(data) {
  data |> mutate(outcome_fact = case_when(
    # Cancer death outcome
    outcome == "cancer death" ~ "Cancer death",
    outcome == "other cause death" ~ "Other cause death",
    # Cancer outcome
    outcome == "cancer" ~ "Cancer",
    outcome == "death" ~ "Any cause death",
    # Censoring,
    outcome == "end" ~ "Censored",
    TRUE ~ outcome
  ),
  outcome_fact = factor(outcome_fact, levels = c(
    "Censored",
    "Cancer death", "Other cause death",
    "Cancer", "Any cause death"
  )))
}


calculate_cuminc <- function(data){
  cuminc(Surv(follow_sa, outcome_fact) ~ case, 
         data=data, conf.level =.95)
}

# Death -----
# From medium-term follow-up, only include those up to 5 years of follow-up originally,
# if longer than 5 years, they will be covered in the long-term follow-up
# *** Same concept applies to the rest of the datasets below ***
cancer_death_med <- readRDS("./data/final_data/new_follow/cancer_death_med_term.RDS") |> filter(follow_5y)
cancer_death_long <- readRDS("./data/final_data/new_follow/cancer_death_long_term.RDS")
cancer_death <- cancer_death_med |> bind_rows(cancer_death_long) |> shift_follow_up_3m() |> change_outcome_to_factor()


# CIF 
cif_death <- 
  cancer_death |>
  mutate(case = recode_case(case)) |>
  group_by(gender) |>
  nest() |>
  mutate(cumincs = map(data, ~calculate_cuminc(.x)
  ))  |>
  select(gender, cumincs) |>
  ungroup()

additional_font_size <-   theme(axis.text.x = element_text(size = 6),
                                axis.text.y = element_text(size = 6),
                                axis.title = element_text(size = 8),
                                legend.text = element_text(size = 7),
                                legend.title = element_text(size = 7.5),
                                panel.spacing.x = unit(2, "lines"),
                                strip.text.x = element_text(size = 7))

saveRDS(cif_death, file = "./data/final_data/cif_death_withCI.rds")

# Cancer ----
cancer_med <- readRDS("./data/final_data/cancer_diag_med_term.RDS") |> filter(follow_5y)
cancer_long <- readRDS("./data/final_data/cancer_diag_long_term.RDS")
cancer <- cancer_med |> bind_rows(cancer_long) |> shift_follow_up_3m() |> change_outcome_to_factor()

# CIF of all events
cif_cancer <- 
  cancer |>
  mutate(case = recode_case(case)) |>
  group_by(gender) |>
  nest() |>
  mutate(cumincs = map(data,  ~calculate_cuminc(.x)))  |>
  select(gender, cumincs) |>
  ungroup()

saveRDS(cif_cancer, file="./data/final_data/cif_cancer_withCI.rds")

# By Body Region ----
cancer_region_med <- readRDS("./data/final_data/cancer_region_med_term.RDS") |> filter(follow_5y)
cancer_region_long <- readRDS("./data/final_data/cancer_region_long_term.RDS")
cancer_region <- cancer_region_med |> bind_rows(cancer_region_long) |> shift_follow_up_3m() |> filter(!excl_analysis) |>
  change_outcome_to_factor()

# CIF 
cif_region  <- 
  cancer_region |>
  mutate(case = recode_case(case)) |>
  group_by(gender, region_desc) |>
  nest() |>
  mutate(cumincs = map(data, ~calculate_cuminc(.x)))  |>
  select(region_desc, cumincs) |>
  ungroup()

saveRDS(cif_region, file="./data/final_data/cif_region_all_withCI.rds")

# By Organ System ----
cancer_organsys_med <- readRDS("./data/final_data/cancer_organsys_med_term.RDS") |> filter(follow_5y)
cancer_organsys_long <- readRDS("./data/final_data/cancer_organsys_long_term.RDS")
cancer_organsys <- cancer_organsys_med |> bind_rows(cancer_organsys_long) |> shift_follow_up_3m()|> filter(!excl_analysis)  |>
  change_outcome_to_factor()

# CIF 
cif_organsys <- 
  cancer_organsys |>
  mutate(case = recode_case(case)) |>
  group_by(gender, organsys_desc) |>
  nest() |>
  mutate(cumincs = map(data, ~calculate_cuminc(.x)))  |>
  select(gender, organsys_desc, cumincs) |>
  ungroup()

saveRDS(cif_organsys, file="./data/final_data/cif_organsys_all_withCI.rds")

write(paste(Sys.time(), "End 11_cumulative_incidence_descriptive.R script"), file = progress_file, append = T)

