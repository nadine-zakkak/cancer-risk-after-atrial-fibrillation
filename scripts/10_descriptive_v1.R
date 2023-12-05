rm(list=ls()); gc()

library(plyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(flextable)
library(stringr)

# all patients
allpatients <- readRDS("./data/final_data/confounders_final.rds")
alc_dict <- read.csv("./look_ups/alcohol_dictionary.csv")

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 10_descriptive_v1.R script"), file = progress_file, append = T)

# load functions, global variables and codelists
source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
source("./scripts/00_read_codelists.R")

# covariates -----
allpatients_long <-
  allpatients |>
  mutate(alc_desc = factor(add_alc_desc(alc, alc_dict), levels = alc_levels),
         imd = ifelse(imd == 1, "1 - Least", ifelse(imd == 5, "5 - Most", imd))) |>
  mutate(imd = as.character(imd),
         diab = as.character(diab),
         ht = as.character(ht),
         alc_desc = as.character(alc_desc),
         smok = as.character(smok)) |>
  dplyr::rename(`IMD (2015)` = imd,
                `Diabetes` = diab,
                `Hypertension` = ht,
                `Alcohol status` = alc_desc,
                `Ex/current smoker` = smok) |>
  pivot_longer(cols = c(`IMD (2015)`, `Diabetes`, `Hypertension`, `Alcohol status`, `Ex/current smoker`), 
               names_to = "var")


summary_age <-
  allpatients |> 
  group_by(gender, case) |>
  summarise(median_age = median(age_start),
            q1_age = quantile(age_start, .25),
            q3_age = quantile(age_start, .75)) |>
  mutate(median_age = format_numb(median_age),
         q1_age = format_numb(q1_age),
         q3_age = format_numb(q3_age),
         n = sprintf("%s (%s, %s)", median_age, q1_age, q3_age)) |>
  mutate(var = "Age in years",
         value = "Median (IQR)") |>
  ungroup() |>
  mutate(gender = recode_gender(gender),
         case = recode_case(case))

summary <-
  allpatients_long |>
  group_by(gender, case, var, value) |>
  summarise(n = n()) |>
  mutate(N = sum(n),
         prop = 100*n/N) |>
  bind_rows(allpatients |>
              group_by(gender, case) |>
              summarise(n = n()) |>
              mutate(
                # prop = 100*n/sum(n),
                var = "Total",
                value = "Total") |>
              ungroup()) |>
  ungroup() |>
  select(gender, case, var, value, N, n, prop) |>
  mutate(gender = recode_gender(gender),
         case = recode_case(case),
         n = format_numb(n),
         prop = format_numb(prop))

summary_wide <-
  summary_age |>
  bind_rows(summary) |>
  mutate(summary = ifelse(is.na(prop),
                          n,
                          stringr::str_glue("{n} ({prop}%)"))) |>
  filter(value != "FALSE") |>
  pivot_wider(id_cols = c(var, value), names_from = c(gender, case), values_from = c(summary)) |>
  mutate(var = factor(var,
                      levels = c("Age in years", "IMD (2015)", 
                                 "Diabetes", "Hypertension",
                                 "Alcohol status", "Ex/current smoker", 
                                 "Total"))) |>
  arrange(var)

write.csv(summary_wide, "./results/covar_desc.csv", row.names = F)

## save to word doc ----
summary_wide |>
  flextable() |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  merge_v(j=1) |>
  hline(i=c(1,6,7,8,12,13)) |>
  add_table_theme() |>
  save_as_docx(path = "./results/covar_desc.docx")

# outcomes -----
outcome_summary <- function(data, outcome, follow_time, var_name){
  outcome <- enquo(outcome)
  return(data |>
           group_by(gender, case, !!outcome, .drop = F) |>
           summarise(n = n()) |>
           mutate(
             follow = follow_time,
             var = var_name,
             value = as.character(!!outcome)) |>
           ungroup())
  
}

clean_outcome <- function(data) {
  data |>
    mutate(gender = recode_gender(gender),
           case = recode_case(case),
           prop = format_numb(100*n/N),
           prop_suppress = ifelse(n<5, format_numb(100*5/N), format_numb(100*n/N)),
           summary = ifelse(n<5, str_glue('n<5 ({prop_suppress}%)'), str_glue('{n} ({prop_suppress}%)'))
    )
}

count_N_patients <- function(data, ...){
  data |>
    group_by(!!!ensyms(...)) |>
    summarise(N = n_distinct(epatid))
}

#' calculate_ci (95% CI) - from Matt
#'
#' @param n - number of samples (numerator)
#' @param N - total number of samples (denominator)
#'
#' @return: a list with 2 values: "lower" = lower CI, "upper" = upper CI
calculate_ci <- function(data){
  data |>
    mutate(lb = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2))),
           ub = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2))),
           lb = ifelse(n<5, "-", format_numb(100*lb)),
           ub = ifelse(n<5, "-", format_numb(100*ub))
    )
}

# cancer death ----
cancer_death_short <- readRDS("./data/final_data/cancer_death_short_term.RDS")
cancer_death_med <- readRDS("./data/final_data/cancer_death_med_term.RDS")
cancer_death_long <- readRDS("./data/final_data/cancer_death_long_term.RDS")

N_counts <- 
  cancer_death_short |> count_N_patients(gender, case) |> mutate(follow = "short") |>
  bind_rows(cancer_death_med |> count_N_patients(gender, case) |> mutate(follow = "med")) |>
  bind_rows(cancer_death_long |> count_N_patients(gender, case) |> mutate(follow = "long"))

cancer_death_summary <-
  outcome_summary(cancer_death_short, cancer_death_outcome, "short", "Cancer Death") |>
  left_join(N_counts, by = c('gender', 'case', 'follow')) |>
  clean_outcome() |>
  calculate_ci() |>
  bind_rows(
    outcome_summary(cancer_death_med, cancer_death_outcome, "med", "Cancer Death") |>
      left_join(N_counts, by = c('gender', 'case', 'follow')) |>
      clean_outcome() |>
      calculate_ci()
  ) |>
  bind_rows(
    outcome_summary(cancer_death_long, cancer_death_outcome, "long", "Cancer Death") |>
      left_join(N_counts, by = c('gender', 'case', 'follow')) |>
      clean_outcome()  |>
      calculate_ci()
  ) |>
  filter(cancer_death_outcome) |>
  pivot_wider(id_cols = c(var, follow, value), names_from = c(gender, case), values_from = c(summary, N, n, prop, lb, ub))

rm(cancer_death_short, cancer_death_med, cancer_death_long, N_counts); gc()

# cancer
cancer_short <- readRDS("./data/final_data/cancer_diag_short_term.RDS")
cancer_med <- readRDS("./data/final_data/cancer_diag_med_term.RDS")
cancer_long <- readRDS("./data/final_data/cancer_diag_long_term.RDS")

N_counts <- 
  cancer_short |> count_N_patients(gender, case) |> mutate(follow = "short") |>
  bind_rows(cancer_med |> count_N_patients(gender, case) |> mutate(follow = "med")) |>
  bind_rows(cancer_long |> count_N_patients(gender, case) |> mutate(follow = "long"))

cancer_summary <-
  outcome_summary(cancer_short, cancer_outcome, "short", "Cancer Incidence") |>
  left_join(N_counts, by = c('gender', 'case', 'follow')) |>
  clean_outcome() |>
  calculate_ci() |>
  bind_rows(
    outcome_summary(cancer_med, cancer_outcome, "med", "Cancer Incidence") |>
      left_join(N_counts, by = c('gender', 'case', 'follow')) |>
      clean_outcome()
  ) |>
  bind_rows(
    outcome_summary(cancer_long, cancer_outcome, "long", "Cancer Incidence") |>
      left_join(N_counts, by = c('gender', 'case', 'follow')) |>
      clean_outcome() |>
      calculate_ci() 
  ) |>
  filter(cancer_outcome) |>
  pivot_wider(id_cols = c(var, follow, value), names_from = c(gender, case), values_from =  c(summary, N, n, prop, lb, ub))

rm(cancer_short, cancer_med, cancer_long, N_counts); gc()

# cancer by body region
cancer_region_short <- readRDS("./data/final_data/cancer_region_short_term.RDS")
cancer_region_med <- readRDS("./data/final_data/cancer_region_med_term.RDS")
cancer_region_long <- readRDS("./data/final_data/cancer_region_long_term.RDS")

N_counts <- 
  cancer_region_short |> count_N_patients(region_desc, gender, case) |> mutate(follow = "short") |>
  bind_rows(cancer_region_med |> count_N_patients(region_desc, gender, case) |> mutate(follow = "med")) |>
  bind_rows(cancer_region_long |> count_N_patients(region_desc, gender, case) |> mutate(follow = "long"))

cancer_region_summary <-
  outcome_summary(cancer_region_short |> filter(cancer_outcome), region_desc, "short", "Body Region") |>
  left_join(N_counts, by = c('region_desc', 'gender', 'case', 'follow')) |>
  clean_outcome() |>
  calculate_ci() |>
  bind_rows(
    outcome_summary(cancer_region_med  |> filter(cancer_outcome), region_desc, "med", "Body Region") |>
      left_join(N_counts, by = c('region_desc', 'gender', 'case', 'follow')) |>
      clean_outcome() |>
      calculate_ci() 
  ) |>
  bind_rows(
    outcome_summary(cancer_region_long |> filter(cancer_outcome), region_desc, "long", "Body Region") |>
      left_join(N_counts, by = c('region_desc', 'gender', 'case', 'follow')) |>
      clean_outcome() |>
      calculate_ci() 
  )  |>
  pivot_wider(id_cols = c(var, follow, region_desc), names_from = c(gender, case), values_from = c(summary, N, n, prop, lb, ub)) |>
  rename(value = region_desc)

rm(cancer_region_short, cancer_region_med, cancer_region_long, N_counts); gc()

# cancer by organ system
cancer_organsys_short <- readRDS("./data/final_data/cancer_organsys_short_term.RDS")
cancer_organsys_med <- readRDS("./data/final_data/cancer_organsys_med_term.RDS")
cancer_organsys_long <- readRDS("./data/final_data/cancer_organsys_long_term.RDS")

N_counts <- 
  cancer_organsys_short |> count_N_patients(organsys_desc, gender, case) |> mutate(follow = "short") |>
  bind_rows(cancer_organsys_med |> count_N_patients(organsys_desc, gender, case) |> mutate(follow = "med")) |>
  bind_rows(cancer_organsys_long |> count_N_patients(organsys_desc, gender, case) |> mutate(follow = "long"))

cancer_organsys_summary <-
  outcome_summary(cancer_organsys_short |> filter(cancer_outcome), organsys_desc, "short", "Organ System") |>
  left_join(N_counts, by = c('organsys_desc', 'gender', 'case', 'follow')) |>
  clean_outcome() |>
  calculate_ci() |>
  bind_rows(
    outcome_summary(cancer_organsys_med  |> filter(cancer_outcome), organsys_desc, "med", "Organ System") |>
      left_join(N_counts, by = c('organsys_desc', 'gender', 'case', 'follow')) |>
      clean_outcome()  |>
      calculate_ci()
  ) |>
  bind_rows(
    outcome_summary(cancer_organsys_long |> filter(cancer_outcome), organsys_desc, "long", "Organ System") |>
      left_join(N_counts, by = c('organsys_desc', 'gender', 'case', 'follow')) |>
      clean_outcome()  |>
      calculate_ci()
  )  |>
  pivot_wider(id_cols = c(var, follow, organsys_desc), names_from = c(gender, case), values_from = c(summary, N, n, prop, lb, ub)) |>
  rename(value = organsys_desc)

outcomes_summary  <-
  cancer_death_summary |>
  bind_rows(cancer_summary) |>
  bind_rows(cancer_region_summary) |>
  bind_rows(cancer_organsys_summary) |>
  mutate(follow = factor(case_when(
    follow == "short" ~ "<= 3 months",
    follow == "med" ~ "3 months - 5 years",
    follow == "long" ~ "> 5 years"
  ),
  c("<= 3 months", "3 months - 5 years", "> 5 years")),
  var = factor(var, levels = c("Cancer Death", "Cancer Incidence", "Body Region", "Organ System"))) |>
  select(follow, var, value, everything()) |>
  arrange(follow, var, value) |>
  group_by(follow, var) |>
  select(follow, var, value, 
         starts_with("summary_Men"), starts_with("summary_Women"),
         starts_with("N_Men"), starts_with("N_Women"),
         starts_with("n_Men"), starts_with("n_Women"),
         starts_with("prop_Men"), starts_with("prop_Women"),
         starts_with("lb_Men"), starts_with("ub_Men"),
         starts_with("lb_Women"), starts_with("ub_Women"))

write.csv(outcomes_summary, "./results/outcome_desc_timesplit.csv", row.names = F)


outcomes_summary |>
  select(!c(starts_with("N"), starts_with("n"), starts_with("prop"), starts_with("lb"), starts_with("ub"))) |>
  flextable() |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(3,2,2)) |>
  merge_v(j=c(1,2)) |>
  hline(i=c(1,2,10,22,23,24,32,44,45,46,54)) |>
  add_table_theme() |>
  save_as_docx(path = "./results/outcome_desc_timesplit.docx")

# follow up ----
## Cancer death
cancer_death <- readRDS("./data/final_data/cancer_death.RDS")

cancer_death_followup <- 
  cancer_death |>
  mutate(cohort=recode_case(case),
         outcome=factor(outcome, c("af", "cancer death", "other cause death", "end"))) |>
  arrange(outcome) |>
  ggplot() +
  geom_boxplot(aes(y=follow/12, x=outcome, fill=cohort), width = .5) +
  scale_fill_manual(values = c("#AA92DD", "#67CBA6")) +
  scale_y_continuous(breaks=seq(0,20,2), labels = seq(0,20,2)) +
  ylab("Follow-up time in years") +
  theme_minimal()
ggsave(plot = cancer_death_followup,
       filename="./results/cancer_death_followup.png", width = 25, height = 13, unit = "cm")

cancer_death_follow_summary <-
  cancer_death |> 
  group_by(gender, case, cancer_died_outcome) |> 
  summarise(median_follow = median(follow/12)) |>
  ungroup() |>
  pivot_wider(id_cols = c(cancer_died_outcome), names_from = c(gender, case), values_from = c(median_follow)) |>
  mutate(var = "Cancer death follow-up",
         cancer_died_outcome = as.factor(cancer_died_outcome)) |>
  rename(value = cancer_died_outcome)

# Cancer diagnosis
cancer <- readRDS("./data/final_data/cancer_diag.RDS")

cancer_followup <-
  cancer |>
  mutate(cohort=recode_case(case),
         outcome=factor(outcome, c("af", "cancer", "death", "end"))) |>
  ggplot() +
  geom_boxplot(aes(y=follow/12, x=outcome, fill=cohort), width = .5) +
  scale_fill_manual(values = c("#AA92DD", "#67CBA6")) +
  scale_y_continuous(breaks=seq(0,20,2), labels = seq(0,20,2)) +
  ylab("Follow-up time in years") +
  theme_minimal()
ggsave(plot=cancer_followup,
       filename="./results/cancer_followup.png", width = 25, height = 13, unit = "cm")

cancer_follow_summary <-
  cancer |> 
  group_by(gender, case, cancer_outcome) |> 
  summarise(median_follow = median(follow/12)) |>
  ungroup() |>
  pivot_wider(id_cols = c(cancer_outcome), names_from = c(gender, case), values_from = c(median_follow)) |>
  mutate(var = "Cancer diagnosis follow up",
         cancer_outcome = as.factor(cancer_outcome)) |>
  rename(value = cancer_outcome)

# By body region
cancer_region <- readRDS("./data/final_data/cancer_region.RDS") |>   
  add_region_desc()

cancer_region_followup <-
  cancer_region |>
  mutate(cohort=recode_case(case),
         outcome_detail =  plyr::mapvalues(outcome_detail,
                                           from = body_region$body_region_numb,
                                           to = body_region$body_region_desc),
         outcome_detail = factor(outcome_detail, 
                                 c(region_order, "af", "death", "end")))|>
  ggplot() +
  geom_boxplot(aes(y=follow/12, x=outcome_detail, fill=cohort), width = .5) +
  scale_fill_manual(values = c("#AA92DD", "#67CBA6")) +
  scale_y_continuous(breaks=seq(0,20,2), labels = seq(0,20,2)) +
  ylab("Follow-up time in years") +
  theme_minimal()
ggsave(plot=cancer_region_followup,
       filename="./results/cancer_region_followup.png", width = 35, height = 20, unit = "cm")

cancer_region_follow_summary <-
  cancer_region |> 
  filter(cancer_outcome) |>
  group_by(gender, case, region_desc) |> 
  summarise(median_follow = median(follow/12)) |>
  ungroup() |>
  pivot_wider(id_cols = c(region_desc), names_from = c(gender, case), values_from = c(median_follow)) |>
  mutate(var = "Cancer Diagnosis by Body Region follow up",
         region_desc = factor(region_desc, levels = region_order)) |>
  arrange(region_desc) |>
  rename(value = region_desc)

## By organ system
cancer_organsys <- readRDS("./data/final_data/cancer_system.RDS")  |>
  add_organsys_desc()

cancer_organsys_followup <-
  cancer_organsys |>
  mutate(cohort=recode_case(case),
         outcome_detail =  plyr::mapvalues(outcome_detail,
                                           from = organ_system$organ_sys_numb,
                                           to = organ_system$organ_sys_desc),
         outcome_detail = factor(outcome_detail, 
                                 c(organsys_order, "af", "death", "end")))|>
  ggplot() +
  geom_boxplot(aes(y=follow/12, x=outcome_detail, fill=cohort), width = .5) +
  scale_fill_manual(values = c("#AA92DD", "#67CBA6")) +
  scale_y_continuous(breaks=seq(0,20,2), labels = seq(0,20,2)) +
  ylab("Follow-up time in years") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
ggsave(plot=cancer_organsys_followup,
       filename="./results/cancer_organsys_followup.png", width = 35, height = 20, unit = "cm")

cancer_organsys_follow_summary <-
  cancer_organsys |> 
  filter(cancer_outcome) |>
  group_by(gender, case, organsys_desc) |> 
  summarise(median_follow = median(follow/12)) |>
  ungroup() |>
  pivot_wider(id_cols = c(organsys_desc), names_from = c(gender, case), values_from = c(median_follow)) |>
  mutate(var = "Cancer diagnosis by organ system follow up",
         organsys_desc = factor(organsys_desc, levels = organsys_order)) |>
  arrange(organsys_desc) |>
  rename(value = organsys_desc)

followup_summary <-
  cancer_death_follow_summary |>
  bind_rows(cancer_follow_summary) |>
  bind_rows(cancer_region_follow_summary) |>
  bind_rows(cancer_organsys_follow_summary) |>
  select(var, value, everything()) |>
  mutate(across(where(is.numeric), ~format_numb(.))) |>
  group_by(var)

write.csv(followup_summary, "./results/followup_desc_detail.csv", row.names = F)

followup_summary |>
  flextable() |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  add_header_row(values = c("", "Median follow-up time"), colwidths = c(2, 4)) |>
  merge_v(j=1) |>
  hline(i=c(2,4,12)) |>
  add_table_theme() |>
  save_as_docx(path = "./results/followup_desc_detail.docx")


followup_summary_cdeath_incidence <- 
  cancer_death |> 
  group_by(gender, case) |> 
  summarise(median_follow = format_numb(median(follow/12)),
            IQR_follow = format_numb(IQR(follow/12)),
            Q1_follow = format_numb(quantile(follow/12, .25)),
            Q3_follow = format_numb(quantile(follow/12, .75)),
            median_IQR = str_glue('{median_follow} ({Q1_follow}-{Q3_follow})')) |>
  ungroup() |>
  mutate(var = "Cancer death analysis") |> 
  bind_rows(
    cancer |> 
      group_by(gender, case) |> 
      summarise(median_follow = format_numb(median(follow/12)),
                IQR_follow = format_numb(IQR(follow/12)),
                Q1_follow = format_numb(quantile(follow/12, .25)),
                Q3_follow = format_numb(quantile(follow/12, .75)),
                median_IQR = str_glue('{median_follow} ({Q1_follow}-{Q3_follow})')) |>
      ungroup() |>
      mutate(var = "Cancer incidence analysis")) |>
  mutate(gender = recode_gender(gender),
         case = recode_case(case)) |>
  pivot_wider(id_cols = c(var), names_from = c(gender, case), values_from = c(median_IQR)) |>
  select(var, everything())

write.csv(followup_summary_cdeath_incidence, "./results/followup_desc.csv", row.names = F)

followup_summary_cdeath_incidence |>
  flextable() |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(1,2,2)) |>
  add_header_row(values = c("", "Median follow-up time (IQR) to any of the possible outcomes"), colwidths = c(1,4)) |>
  add_table_theme() |>
  save_as_docx(path = "./results/followup_desc.docx")

write(paste(Sys.time(), "End 10_descriptive_v1.R script"), file = progress_file, append = T)


