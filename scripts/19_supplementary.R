# Supplementary figures 
#"Descriptive"
#
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

path_save <- "./tables_figs/supplementary/"

# Supplementary table 1  -----
# Outcomes proportion by time-split
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")

outcome_desc <- outcome_desc |>  select(!c(starts_with("N"), starts_with("n"), starts_with("prop"), starts_with("lb"), starts_with("ub"))) |>
  filter(!var %in% c("Cancer Death", "Cancer Incidence")) 

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
  hline(i=c(8,20,28,40,48)) |>
  add_table_theme() |>
  save_as_docx(path =  paste0(path_save, "supplementary_table1.docx"))

# % of cancer at different time points not stratified by sex -----
#  3m - 5 years:
#  > 5 years:
# function
shift_follow_up_3m <- function(data){
  data |>
    mutate(follow_sa = time_length(start_orig %m+% months(3) %--% outcome_date, unit = "month"))
}

# Cancer ----
cancer_med <- readRDS("./data/final_data/cancer_diag_med_term.RDS") |> filter(follow_5y)
cancer_long <- readRDS("./data/final_data/cancer_diag_long_term.RDS")
cancer <- cancer_med |> bind_rows(cancer_long) |> 
  shift_follow_up_3m() |>
  mutate(case = recode_case(case),
         outcome_fact = case_when(outcome == "cancer" ~ "Cancer", 
                                  outcome == "death" ~ "Any cause death",
                                  outcome == "end" ~ "Censored",
                                  TRUE ~ outcome),
         outcome_fact = as.factor(outcome_fact),
         outcome_fact = factor(outcome_fact, levels = c("Censored", "Cancer", "Any cause death"))) |>
  select(follow_sa, outcome_fact, case)

cif_cancer <-  cuminc(Surv(follow_sa, outcome_fact) ~ case, 
                      data=cancer, conf.level =.95)
cif_cancer_df <- cif_cancer$tidy |> 
  mutate(case = factor(strata, levels = c("Case", "Control")),
         outcome_mod = case_when( outcome == "cancer" ~ "Cancer",
                                  outcome == "death" ~ "Any cause death",
                                  outcome == "end" ~ "Censored",
                                  TRUE ~ outcome),
         outcome = factor(outcome, levels = c("Cancer", "Any cause death", "Censored")))

# Get cumulative incidence estimates at 5 and 10 years
cif_estimates <- cif_cancer_df |> 
  filter(outcome_mod == "Cancer", case   == "Case") |>
  mutate(time_shifted = time + time_shift,
         estimate  = 100*estimate ) |>
  filter(abs(time_shifted - 60) == min(abs(time_shifted - 60)) |
           abs(time_shifted - 120) == min(abs(time_shifted - 120))) |>
  mutate(time_shifted = round(time_shifted)) |>
  group_by(time_shifted, outcome_mod, case) |>
  summarise(est  = max(estimate )) |>
  arrange(time_shifted)

# Proportion of cancer diagnosed within AF patients across different follow-up periods
# < 3m: proportion of all AF patients diagnosed with cancer
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")
summary_df <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Cancer Incidence", follow == "<= 3 months") |>
  mutate(var = case_when(
    var == "Cancer Incidence" ~ "Cancer incidence analysis"
  )) |>
  mutate(est = (n_Women_Case+n_Men_Case)/(N_Women_Case+N_Men_Case)*100) |>
  select(var, est)

summary_df <- 
  summary_df |>
  mutate(follow = "3 months") |>
  bind_rows(cif_estimates |>
              mutate(follow = case_when(
                time_shifted == 60 ~ "5 years",
                time_shifted == 120 ~ "10 years"
              )) |>
              select(-time_shifted, -case))


summary_df <- summary_df |> 
  mutate(follow = factor(follow, levels = rev(c("3 months", "5 years", "10 years"))))

# Supplementary figure 1 ----

p1 <- 
  ggplot() +
  geom_bar(aes(y = follow, x = est, fill = follow), data = summary_df, 
           stat = "identity", position=position_dodge(), width = .8) +
  geom_text(aes(label = sprintf("%g%%", round(est, 1)), x = est, y = follow), 
            data = summary_df, hjust=-0.1, size = 3) +
  theme_minimal() +
  scale_x_continuous(expand=expansion(mult=c(0,0.1)))+
  xlim(c(0,50)) +
  scale_fill_brewer(palette = "Blues", type = "seq") +
  ylab("Follow-up time") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        aspect.ratio = .5/1)

do.call(ggsave, c(list(plot = p1, filename = paste0(path_save, "supp1.jpeg")), list(device="jpeg",
                                                                                    unit="mm",
                                                                                    width = 160,
                                                                                    height = 120,
                                                                                    dpi = 1000)))

# Supplementary figure 2 -----
colours <- c("#2C4C63", "#809BCE", "#53A276", "#B8E0D2",
             "#C45A8A", "#ECC9D9", "#5E37A6", "#C3B5DE")
# Cancer < 3 months after new-onset AF by body region
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")
N_cancer_3m <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Cancer Incidence", follow == "<= 3 months") |>
  mutate(est = n_Women_Case+n_Men_Case) |>
  select(var, est) |> pull(est)

# Proportion of cancer diagnosed within AF patients across different follow-up periods
# < 3m
region_summary_df <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Body Region", follow == "<= 3 months") |> 
  mutate(est = (n_Women_Case+n_Men_Case)/(N_cancer_3m)*100) |>
  select(var, value, est) |>
  arrange(desc(est)) |>
  mutate(value = factor(value)) |>
  arrange(value)


p2 <- 
  ggplot() +
  geom_bar(aes(y = reorder(value, est), x = est, fill = reorder(value, est)), data = region_summary_df, 
           stat = "identity", position=position_dodge(), width = .8) +
  geom_text(aes(label = sprintf("%g%%", round(est, 1)), x = est, y = reorder(value, est)),
            data = region_summary_df, hjust=-0.1, size = 3) +
  theme_minimal() +
  scale_x_continuous(expand=expansion(mult=c(0,0.1)))+
  xlim(c(0,50)) +
  scale_fill_manual(values=colours) +
  ylab("Body Region") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        aspect.ratio = .5/1)

do.call(ggsave, c(list(plot = p2, filename = paste0(path_save, "supp2.jpeg")), list(device="jpeg",
                                                                                    unit="mm",
                                                                                    width = 160,
                                                                                    height = 120,
                                                                                    dpi = 1000)))

# Supplementary figure 3 -----
# Cancer < 3 months after new-onset AF by body region
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")
N_cancer_3m <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Cancer Incidence", follow == "<= 3 months") |>
  mutate(est = n_Women_Case+n_Men_Case) |>
  select(var, est) |> pull(est)

# Proportion of cancer diagnosed within AF patients across different follow-up periods
# < 3m
organsys_summary_df <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Organ System", follow == "<= 3 months") |> 
  mutate(est = (n_Women_Case+n_Men_Case)/(N_cancer_3m)*100) |>
  select(var, value, est) |>
  arrange(desc(est)) |>
  mutate(value = factor(value)) |>
  arrange(value)

colours <- c("#2C4C63", "#809BCE", "#FDBF6F", "#FF7F00", 
             "#B92929", "#F19393", "#53A276", "#B8E0D2",
             "#C45A8A", "#ECC9D9", "#5E37A6", "#C3B5DE")

p3 <- 
  ggplot() +
  geom_bar(aes(y = reorder(value, est), x = est, fill = reorder(value, est)), data = organsys_summary_df, 
           stat = "identity", position=position_dodge(), width = .8) +
  geom_text(aes(label = sprintf("%g%%", round(est, 1)),
                x = est, y = reorder(value, est)), data = organsys_summary_df, 
            hjust=-0.1,
            size = 3) +
  theme_minimal() +
  scale_x_continuous(expand=expansion(mult=c(0,0.1)))+
  xlim(c(0,50)) +
  scale_fill_manual(values=colours) +
  ylab("Organ System") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        aspect.ratio = .5/1)

do.call(ggsave, c(list(plot = p3, filename = paste0(path_save, "supp3.jpeg")), list(device="jpeg",
                                                                                    unit="mm",
                                                                                    width = 160,
                                                                                    height = 120,
                                                                                    dpi = 1000)))
