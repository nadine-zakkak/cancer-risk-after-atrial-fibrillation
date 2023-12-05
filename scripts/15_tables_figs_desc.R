library(dplyr)
library(flextable)

source("./scripts/00_global_var.R")
source("./scripts/00_functions.R")

path_save <- "./tables_figs/"

# General descriptives ----
allpatients <- readRDS("./data/final_data/confounders_final.rds")

median_age <- allpatients |>
  summarise(median_age = median(age_start),
            q1_age = quantile(age_start, .25),
            q3_age = quantile(age_start, .75)) |>
  mutate(cohort = "All") |>
  bind_rows(allpatients |>
              mutate(cohort = recode_gender(gender)) |>
              group_by(cohort) |>
              summarise(median_age = median(age_start),
                        q1_age = quantile(age_start, .25),
                        q3_age = quantile(age_start, .75))) |> 
  bind_rows(allpatients |>
              mutate(cohort = recode_case(case)) |>
              group_by(cohort) |>
              summarise(median_age = median(age_start),
                        q1_age = quantile(age_start, .25),
                        q3_age = quantile(age_start, .75))) |>
  bind_rows(allpatients |>
              mutate(gender_lb = recode_gender(gender),
                     case_lb = recode_case(case),
                     cohort = paste0(gender_lb, "-", case_lb)) |>
              group_by(cohort) |>
              summarise(median_age = median(age_start),
                        q1_age = quantile(age_start, .25),
                        q3_age = quantile(age_start, .75))) |>
  mutate(median_age = format_numb(median_age),
         q1_age = format_numb(q1_age),
         q3_age = format_numb(q3_age),
         n = sprintf("%s (%s, %s)", median_age, q1_age, q3_age)) |>
  select(cohort, n) |>
  rename("Age (years) - Median (IQR)"= n)

write.csv(median_age, paste0(path_save, "median_age.csv"))

proportion_sex <- 
  allpatients |>
  mutate(gender_lb = recode_gender(gender)) |>
  group_by(gender_lb) |>
  summarise(n=n()) |> mutate(N=sum(n), prop = 100*n/N)

write.csv(proportion_sex, paste0(path_save, "proportion_sex.csv"))



# Table 1 -----
# Descriptive by covariates
# csv file produce in 10_descriptive_v1.R
cov_desc <- read.csv("./results/covar_desc.csv")

cov_desc |>
  flextable() |>
  set_header_labels(var= "", value = "",
                    Men_Control = "Control",
                    Women_Control = "Control",
                    Men_Case = "Case",
                    Women_Case = "Case") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  merge_v(j=1) |>
  hline(i=c(1,6,7,8,12,13)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "Table1.docx"))

# Table 2  -----
followup_summary <- read.csv("./results/followup_desc.csv") |>
  mutate(val = "Median follow-up time (IQR) to any of the possible outcomes") 
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")
outcome_desc <- outcome_desc |>  select(!c(starts_with("N"), starts_with("n"), starts_with("prop"))) |>
  filter(var %in% c("Cancer Death", "Cancer Incidence")) |>
  mutate(var = case_when(
    var == "Cancer Death" ~ "Cancer death analysis",
    var == "Cancer Incidence" ~ "Cancer incidence analysis"
  )) |>
  select(-value) |>
  rename(val = follow) |>
  rename_with(.fn = ~gsub("summary_", "", .x), 
              cols = starts_with("summary_"))

table2 <-
  followup_summary |>
  bind_rows(outcome_desc) |>
  arrange(var) |>
  select(var, val, everything())

table2 |>
  select(var, val, starts_with("Men"), starts_with("Women")) |>
  flextable() |>
  merge_v(j=c(1,2)) |>
  set_header_labels(var= "", value = "",
                    Men_Control = "Control",
                    Women_Control = "Control",
                    Men_Case = "Case",
                    Women_Case = "Case") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "Table2.docx"))

# Cumulative incidence curves -----
# Draw cumulative incidence curves starting after 3 months
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

# cancer death and incidence
save_image <- list(device = "jpeg", unit = "mm", width = 140, height = 80, dpi=1000)
additional_font_size <-   theme(axis.text.x = element_text(size = 6),
                                axis.text.y = element_text(size = 6),
                                axis.title = element_text(size = 8),
                                legend.text = element_text(size = 7),
                                legend.title = element_text(size = 7.5),
                                panel.spacing.x = unit(2, "lines"),
                                strip.text.x = element_text(size = 7))

# cancer region and organ system
save_image_large <- list(device = "jpeg", unit = "mm", width = 180, height = 70, dpi=1000)
additional_theme <-
  theme(axis.text.x = element_text(size = 6, margin=margin(t=0)),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 5),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        panel.spacing = unit(1, "lines"))

single_cif_to_df <- function(cif) {
  # change cuminc object into a dataframe
  cif$tidy
}

plot_cif <- function(cif, group_vars, rows=NULL, cols=NULL, eventtype=NULL) {
  cif_to_df <- function(cif, group_vars) {
    # change nested dataframe with cuminc objects into an unnested dataframe
    t <- cif |> 
      group_by(!!!syms(group_vars)) |>
      mutate(cif = map(cumincs, .f = ~single_cif_to_df(.x))) |>
      select(!!!syms(group_vars), cif) |>
      unnest(cif) |>
      mutate(gender = recode_gender(gender),
             case = factor(strata, levels = c("Case", "Control")),
             outcome_mod = case_when(
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
             outcome = factor(outcome, levels = c(
               "Cancer death", "Other cause death",
               "Cancer", "Any cause death",
               "Censored"
             )))
    
    region_exists <- cif |> select(any_of("region_desc")) |> ncol()
    if(region_exists == 1) {
      t <- t |> mutate(region_desc = factor(region_desc, levels = region_order))
    }
    organsys_exists <- cif |> select(any_of("organsys_desc")) |> ncol()
    if(organsys_exists == 1){
      t <- t |> mutate(organsys_desc = factor(organsys_desc, levels = organsys_order))
    }
    return(t)
  }
  
  cif_df <- cif_to_df(cif, group_vars)
  if(!missing(eventtype)){
    cif_df <- cif_df |> mutate(lower_outcome = tolower(outcome)) |> filter(lower_outcome == tolower(eventtype)) |> select(-lower_outcome)
  }
  max_time <- 12*20
  
  cif_df |>
    mutate(time_shifted = time + time_shift) |>
    ggplot() +
    # shift time 3 months back because shifted the start date by 3 months
    geom_line(aes(x=time_shifted, y = 100*estimate, col=outcome, linetype=case), linewidth=.2) + 
    geom_ribbon(aes(x = time_shifted, ymin=100*conf.low, ymax=100*conf.high, fill=outcome, linetype=case), alpha = .2) +
    scale_x_continuous(
      limits = c(0, max_time),
      minor_breaks = 12*c(seq(5, max_time, 5)),
      breaks = 12*c(0.25, 10, 20),
      labels = c("3m", 10, 20),
      expand = c(0,0)
    ) +
    scale_linetype_manual(values = c("solid", "dotdash"), name = "Case") +
    scale_colour_manual(values = c(outcome_colours$outcome, outcome_colours$compete), name = "Event type") +
    scale_fill_manual(values = c(outcome_colours$outcome, outcome_colours$compete), name = "Event type") +
    scale_y_continuous(expand = expansion(mult=c(0.01, 0.1))) +
    guides(linetype=guide_legend(override.aes=list(fill=NA))) +
    facet_grid(cols = cols, rows = rows, labeller = label_wrap_gen(width = 15)) +
    labs(x = "Time since start of follow-up (years)",
         y = "Cumulative Incidence (%)",
         col = "Event type", 
         fill = "Event type",
         linetype = "Case") +
    theme(
      panel.background = element_rect(fill="white"),
      panel.grid.major = element_line(linewidth=.2, colour = "#F5F5F5"),
      panel.grid.minor = element_line(linewidth=.1, colour = "#F5F5F5"),
      axis.line = element_line(colour = "black", linewidth=.1),
      legend.key = element_rect(fill = "white"),
      legend.spacing.y = unit(1, "mm"),
      strip.background  = element_rect(fill = "#F5F5F5", color = "#F5F5F5"),
      strip.text = element_text(colour = "black"),
      axis.title.x = element_text(vjust = -1),
      axis.ticks.x = element_line(linewidth=.1),
      axis.ticks.length.x = unit(.5, "mm"),
      axis.ticks.y = element_blank()) 
}                                                                                                                                                                 

## Death -----
cif_death <- readRDS(file = "./data/final_data/cif_death_withCI.rds")

cif_death_plot_censor <-
  plot_cif(cif_death, group_vars="gender", cols=vars(gender)) +
  additional_font_size

do.call(ggsave, c(list(plot = cif_death_plot_censor, filename = paste0(path_save, "cif_death_withCI.jpeg")), save_image))

## Cancer ----
cif_cancer <- readRDS(file="./data/final_data/cif_cancer_withCI.rds")

cif_cancer_plot_censor <- 
  plot_cif(cif_cancer, group_vars="gender", cols=vars(gender)) +
  additional_font_size

do.call(ggsave, c(list(plot = cif_cancer_plot_censor, filename = paste0(path_save, "cif_cancer_withCI.jpeg")), save_image))

## By Body Region ----
cif_region <- readRDS(file="./data/final_data/cif_region_all_withCI.rds")

cif_region_all_plot <- plot_cif(cif_region, group_vars = c("gender", "region_desc"),  cols=vars(region_desc), rows=vars(gender),)  + 
  additional_theme
do.call(ggsave, c(list(plot = cif_region_all_plot, filename = paste0(path_save, "cif_region_all_withCI.jpeg")), save_image_large))

cif_region_cancer_plot <- plot_cif(cif_region, 
                                   group_vars = c("gender", "region_desc"), 
                                   cols=vars(region_desc),
                                   rows=vars(gender),
                                   eventtype = "cancer")  + additional_theme
do.call(ggsave, c(list(plot = cif_region_cancer_plot, filename = paste0(path_save, "cif_region_cancer_withCI.jpeg")), save_image_large))

## By Organ System ----
cif_organsys <- readRDS(file="./data/final_data/cif_organsys_all_withCI.rds")

cif_organsys_all_plot <- plot_cif(cif_organsys, 
                                  group_vars = c("gender", "organsys_desc"), cols=vars(organsys_desc), rows=vars(gender)) +
  additional_theme
do.call(ggsave, c(list(plot = cif_organsys_all_plot, filename = paste0(path_save, "cif_organsys_all_withCI.jpeg")), save_image_large))

cif_organsys_cancer_plot <- plot_cif(cif_organsys, 
                                     group_vars = c("gender", "organsys_desc"), 
                                     cols=vars(organsys_desc),
                                     rows=vars(gender),
                                     eventtype = "cancer") +
  additional_theme
do.call(ggsave, c(list(plot = cif_organsys_cancer_plot, filename = paste0(path_save, "cif_organsys_cancer_withCI.jpeg")), save_image_large))

# Table 3  -----
# Outcomes proportion by time-split
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")

outcome_desc <- outcome_desc |>  select(!c(starts_with("N"), starts_with("n"), starts_with("prop"), starts_with("lb"), starts_with("ub")))

outcome_desc |>
  flextable() |>
  set_header_labels(follow = "Follow-up time",
                    var= "", value = "",
                    summary_Men_Control = "Control",
                    summary_Women_Control = "Control",
                    summary_Men_Case = "Case",
                    summary_Women_Case = "Case") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(3,2,2)) |>
  merge_v(j=c(1,2)) |>
  hline(i=c(1,2,10,22,23,24,32,44,45,46,54)) |>
  add_table_theme() |>
  save_as_docx(path =  paste0(path_save, "Table3.docx"))





