rm(list=ls()); gc()
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(broom)
library(sandwich) 
library(lmtest)
library(ggplot2)

#load functions, global variables and codelists
source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
source("./scripts/00_read_codelists.R")

# global var
allcov <- c("case", "age_start", "diab", "ht", "smok")

common_cols <- c("gender", "case", "follow_3m", "outcome")
path_save <- "./results/models/"
if(!dir.exists(path_save)) {dir.create(path_save)}

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 12_analysis_short_v1.R script"), file = progress_file, append = T)

# Read files
cancer_death <- readRDS("./data/final_data/cancer_death_short_term.RDS") 
cancer <- readRDS("./data/final_data/cancer_diag_short_term.RDS") 
cancer_region <- readRDS("./data/final_data/cancer_region_short_term.RDS")
cancer_organsys <- readRDS("./data/final_data/cancer_organsys_short_term.RDS")

# keep only required columns for analysis
cancer_death <- cancer_death |> 
  select(any_of(common_cols), any_of(allcov), outcome_bin=cancer_death_outcome)

cancer <- cancer |>
  select(any_of(common_cols), any_of(allcov),  outcome_bin=cancer_outcome)

cancer_region <- cancer_region |>
  filter(!excl_analysis) |>
  select(any_of(common_cols), any_of(allcov), cancer_bodyregion, region_desc, outcome_bin=cancer_outcome)

cancer_organsys <- cancer_organsys |>
  filter(!excl_analysis) |>
  select(any_of(common_cols), any_of(allcov), alc_desc, cancer_organsys, organsys_desc, outcome_bin=cancer_outcome)
gc()

# Functions
make_formula <- function(outcome, cov){
  # Convert into a formula that can be used for glm regression
  # input, outcome (string) and covariates (as a string or vector)
  # output: outcome ~ cov1 + cov2 + ...
  as.formula(paste(outcome, "~", paste(cov, collapse = "+")))
}

get_cov <- function(df){
  # Extract the covariates included in the model from the tidied dataframe
  str_replace_all(string = 
                    paste(setdiff(unique(df$term), c("(Intercept)", "caseTRUE")), 
                          collapse = ", "),
                  pattern = "TRUE", 
                  "")
}

forest_plot <- function(df, y, x = log(RR), xmin = log(lower), xmax = log(upper), facet, xlab, ylab) {
  # Draw a scatter plot
  # Input: dataframe, x, y, xmin, xmax, upper and lower bar range, facet grid formula, x and y labels
  y <- enquo(y)
  x <- enquo(x)
  xmin <- enquo(xmin)
  xmax <- enquo(xmax)
  df |>
    ggplot() +
    geom_pointrange(aes(y = !!y, x = !!x, xmin = !!xmin, xmax = !!xmax)) +
    facet_grid(as.formula(facet)) +
    labs(x=xlab, y = ylab) +
    theme_bw()
}

#Exponentiate a number and format it to be displayed neatly
format_exp <- function(x) {format_numb(exp(x))}

# Poisson Regression
# Cancer Death ----
cancer_death_model <-
  cancer_death |>
  group_by(gender) |>
  nest() |>
  mutate(
    # Adjusted for age only
    glms = map(data, ~glm(make_formula("outcome_bin",  c("case", "age_start")),
                          family = poisson(link = 'log'),
                          data = .x)),
    # Adjusted for all covariates
    glms_adjust = map(data, ~glm(make_formula("outcome_bin", allcov),
                                 family = poisson(link = 'log'),
                                 data = .x))
  ) |>
  select(gender, glms, glms_adjust) |>
  rowwise() |>
  mutate(gender = recode_gender(gender)) |>
  ungroup()

save(cancer_death_model, file=paste0(path_save, "short_cancer_death.RData"))


cancer_death_model_tidy <- tidy_poisson_model(cancer_death_model, "gender")
write.csv(cancer_death_model_tidy, paste0(path_save, "short_cancer_death.csv"), row.names = F)

cancer_death_model_tidy_rbst <- tidy_poisson_model_robust(cancer_death_model, "gender")
write.csv(cancer_death_model_tidy_rbst, paste0(path_save, "short_cancer_death_rbst.csv"), row.names = F)

# Adjusted for all covariates
cancer_death_adj_model <- cancer_death_model |> select(-glms, glms=glms_adjust)

cancer_death_adj_model_tidy <- tidy_poisson_model(cancer_death_adj_model, "gender") 
write.csv(cancer_death_adj_model_tidy, paste0(path_save, "short_cancer_death_adj.csv"), row.names = F)

cancer_death_adj_model_tidy_rbst <- tidy_poisson_model_robust(cancer_death_adj_model, "gender")
write.csv(cancer_death_adj_model_tidy_rbst, paste0(path_save, "short_cancer_death_adj_rbst.csv"), row.names = F)

cancer_death_adj_model_tidy_rbst_plt <- 
  forest_plot(cancer_death_adj_model_tidy_rbst |>
                filter(grepl("case", term)) |>
                mutate(outcome = ""),
              y = outcome, x = log(RR), xmin = log(lower), xmax = log(upper),
              facet = "~gender",
              xlab = "Relative risk of cancer deaths in cases vs controls ", ylab = "Cancer Death") + labs(
                caption = paste("adjusted for", get_cov(cancer_death_adj_model_tidy_rbst))) +
  scale_x_continuous(labels = format_exp)
ggsave(cancer_death_adj_model_tidy_rbst_plt, file = paste0(path_save, "short_cancer_death_adj_rbst.png"),
       width = 20, height = 15, unit = "cm")

# Cancer Incidence ----
# Adjusted only for age 
cancer_model <-
  cancer |>
  group_by(gender) |>
  nest() |>
  mutate(
    # Adjusted only for age 
    glms = map(
      data, ~glm(make_formula("outcome_bin", c("case", "age_start")),
                 family = poisson(link = "log"),
                 data = .x)),
    # Adjusted for all covariates
    glms_adjust = map(
      data, ~glm(make_formula("outcome_bin", allcov),
                 family = poisson(link = "log"),
                 data = .x))
  ) |>
  select(gender, glms, glms_adjust) |>
  rowwise() |>
  mutate(gender = recode_gender(gender)) |>
  ungroup()

save(cancer_model, file=paste0(path_save, "short_cancer.RData"))

cancer_model_tidy <- tidy_poisson_model(cancer_model, "gender")
write.csv(cancer_model_tidy, paste0(path_save, "short_cancer.csv"), row.names = F)

cancer_model_tidy_rbst <- tidy_poisson_model_robust(cancer_model, "gender")
write.csv(cancer_model_tidy_rbst, paste0(path_save, "short_cancer_rbst.csv"), row.names = F)

# Adjusted for all covariates
cancer_adj_model <- cancer_model |> select(-glms, glms=glms_adjust)

cancer_adj_model_tidy <- tidy_poisson_model(cancer_adj_model, "gender") 
write.csv(cancer_adj_model_tidy, paste0(path_save, "short_cancer_adj.csv"), row.names = F)

cancer_adj_model_tidy_rbst <- tidy_poisson_model_robust(cancer_adj_model, "gender")
write.csv(cancer_adj_model_tidy_rbst, paste0(path_save, "short_cancer_adj_rbst.csv"), row.names = F)

cancer_adj_model_tidy_rbst_plt <- 
  forest_plot(cancer_adj_model_tidy_rbst |>
                filter(grepl("case", term)) |>
                mutate(outcome = ""),
              y = outcome,
              facet = "~gender",
              xlab = "Relative risk of cancer in cases vs controls ", ylab = "Cancer") + 
  labs(caption = paste("adjusted for", get_cov(cancer_adj_model_tidy_rbst))) +
  scale_x_continuous(labels = format_exp)
ggsave(cancer_adj_model_tidy_rbst_plt, file = paste0(path_save, "short_cancer_adj_rbst.png"),
       width = 20, height = 15, unit = "cm")

# Body Region ----
cancer_region_model <-
  cancer_region |>
  filter(!cancer_bodyregion %in% c(5,6)) |> #Exclude if Brain/CNS or Upper and lower limbs
  group_by(gender, cancer_bodyregion) |>
  nest() |>
  mutate(
    # Adjusted only for age 
    glms = map(
      data, ~glm(make_formula("outcome_bin", c("case", "age_start")),
                 family = poisson(link = "log"),
                 data = .x)),
    # Adjusted for all covariates
    glms_adjust = map(
      data, ~glm(make_formula("outcome_bin", 
                              # setdiff(allcov, "alc_desc")
                              allcov),
                 family = poisson(link = "log"),
                 data = .x)
    )
  ) |>
  select(gender, cancer_bodyregion, glms, glms_adjust) |>
  rowwise() |>
  mutate(gender = recode_gender(gender)) |>
  ungroup()

save(cancer_region_model, file=paste0(path_save, "short_region.RData"))

cancer_region_model_tidy <- tidy_poisson_model(cancer_region_model, c("gender", "cancer_bodyregion")) |>
  add_region_desc() |>
  select(gender, region_desc, everything()) |>
  select(-cancer_bodyregion)

write.csv(cancer_region_model_tidy, paste0(path_save, "short_region.csv"), row.names = F)

cancer_region_model_tidy_rbst <- tidy_poisson_model_robust(cancer_region_model, c("gender", "cancer_bodyregion")) |>
  add_region_desc() |>
  select(gender, region_desc, everything()) |>
  select(-cancer_bodyregion)

write.csv(cancer_region_model_tidy_rbst, paste0(path_save, "short_region_rbst.csv"), row.names = F)

# Adjusted for all covariates
cancer_region_adj_model <- cancer_region_model |> select(-glms, glms=glms_adjust)

cancer_region_adj_model_tidy <- tidy_poisson_model(cancer_region_adj_model, c("gender", "cancer_bodyregion")) |>
  add_region_desc() |>
  select(gender, region_desc, everything()) |>
  select(-cancer_bodyregion)

write.csv(cancer_region_adj_model_tidy, paste0(path_save, "short_region_adj.csv"), row.names = F)

cancer_region_adj_model_tidy_rbst <- tidy_poisson_model_robust(cancer_region_adj_model, c("gender", "cancer_bodyregion")) |>
  add_region_desc() |>
  select(gender, region_desc, everything()) |>
  select(-cancer_bodyregion)

write.csv(cancer_region_adj_model_tidy_rbst, paste0(path_save, "short_region_adj_rbst.csv"), row.names = F)

cancer_region_adj_model_tidy_rbst_plt <- 
  forest_plot(cancer_region_adj_model_tidy_rbst |>
                filter(grepl("case", term)) |>
                mutate(region_desc = factor(region_desc, levels = region_order)),
              y = region_desc,
              facet = "~gender",
              xlab = "Relative risk of cancer in cases vs controls ", ylab = "Body Region") + 
  labs(caption = paste("adjusted for", get_cov(cancer_region_adj_model_tidy_rbst))) +
  scale_x_continuous(labels = format_exp)

ggsave(cancer_region_adj_model_tidy_rbst_plt, file = paste0(path_save, "short_region_adj_rbst.png"),
       width = 20, height = 15, unit = "cm")

# Organ System ----
cancer_orgaynsys_model <-
  cancer_organsys |>
  filter(!(cancer_organsys %in% c(1, 3, 4)   #Exclude skeletal, endocrine and cardiovascular from analyses
           | (gender == 0 & cancer_organsys == 10) #Exclude men with breast cancer from analyses
  ) 
  ) |> 
  group_by(gender, cancer_organsys) |>
  nest() |>
  mutate(
    # Adjusted only for age
    glms = map(data, ~glm(make_formula("outcome_bin", c("case", "age_start")),
                          family = poisson(link = "log"),
                          data = .x)),
    # Adjusted for all covariates
    glms_adjust = map_if(.x = data,
                         .p = cancer_organsys != 7,
                         .f = ~glm(make_formula("outcome_bin", allcov),
                                   family = poisson(link = "log"),
                                   data = .x),
                         .else = ~glm(make_formula("outcome_bin", c(allcov, "alc_desc")),
                                      family = poisson(link = "log"),
                                      data = .x)
    )
  ) |>
  select(gender, cancer_organsys, glms, glms_adjust) |>
  rowwise() |>
  mutate(gender = recode_gender(gender)) |>
  ungroup()

save(cancer_orgaynsys_model, file=paste0(path_save, "short_organsys.RData"))

cancer_orgaynsys_model_tidy <- tidy_poisson_model(cancer_orgaynsys_model, c("gender", "cancer_organsys")) |>
  add_organsys_desc() |>
  select(gender, organsys_desc, everything()) |>
  select(-cancer_organsys)

write.csv(cancer_orgaynsys_model_tidy, paste0(path_save, "short_organsys.csv"), row.names = F)

cancer_orgaynsys_model_tidy_rbst <- tidy_poisson_model_robust(cancer_orgaynsys_model, c("gender", "cancer_organsys")) |>
  add_organsys_desc() |>
  select(gender, organsys_desc, everything()) |>
  select(-cancer_organsys)

write.csv(cancer_orgaynsys_model_tidy_rbst, paste0(path_save, "short_organsys_rbst.csv"), row.names = F)

# Adjusted for all covariates
cancer_orgaynsys_adj_model <- cancer_orgaynsys_model |> select(-glms, glms=glms_adjust)

cancer_organsys_adj_model_tidy <- tidy_poisson_model(cancer_orgaynsys_adj_model, c("gender", "cancer_organsys")) |>
  add_organsys_desc() |>
  select(gender, organsys_desc, everything()) |>
  select(-cancer_organsys)

write.csv(cancer_organsys_adj_model_tidy, paste0(path_save, "short_organsys_adj.csv"), row.names = F)

cancer_organsys_adj_model_tidy_rbst <- tidy_poisson_model_robust(cancer_orgaynsys_adj_model, c("gender", "cancer_organsys")) |>
  add_organsys_desc() |>
  select(gender, organsys_desc, everything()) |>
  select(-cancer_organsys)

write.csv(cancer_organsys_adj_model_tidy_rbst, paste0(path_save, "short_organsys_adj_rbst.csv"), row.names = F)

cancer_organsys_adj_model_tidy_rbst_plt <- 
  forest_plot(cancer_organsys_adj_model_tidy_rbst |>
                filter(grepl("case", term)) |>
                mutate(organsys_desc = factor(organsys_desc, levels = organsys_order)),
              y = organsys_desc, x = log(RR), xmin = log(lower), xmax = log(upper),
              facet = "~gender",
              xlab = "Relative risk of cancer in cases vs controls ", ylab = "Organ System") + labs(
                caption = paste("adjusted for", get_cov(cancer_organsys_adj_model_tidy_rbst))) +
  scale_x_continuous(labels = format_exp)

ggsave(cancer_organsys_adj_model_tidy_rbst_plt, file = paste0(path_save, "short_organsys_adj_rbst.png"),
       width = 20, height = 15, unit = "cm")

write(paste(Sys.time(), "End 12_analysis_short_v1.R script"), file = progress_file, append = T)







