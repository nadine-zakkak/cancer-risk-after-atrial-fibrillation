# October/2023
# Figures for graphical abstract

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(flextable)

source("./scripts/00_global_var.R")
source("./scripts/00_functions.R")

path_short <- "./results/models/"
path_med_long <- "./results/final_results_med_long/"
path_save <- "./tables_figs/stata/"

# Analysis figure -----
# Cancer death 
short <- read.csv(paste0(path_short, "short_cancer_death_adj_rbst.csv"))
med <- read.csv(paste0(path_med_long, "med_cancer_death.csv")) 
long <- read.csv(paste0(path_med_long, "long_cancer_death.csv")) 

# Merge results into 1 dataframe (only "case" term)
all_death <- short |> filter(grepl("case.*", term)) |> mutate(followup = "<= 3 months", estimand="RR") |> select(followup, gender, term, estimand, estimate=RR, lower, upper, p) |>
  bind_rows(med |> mutate(followup = "3 months to 5 years", estimand="sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  bind_rows(long |> mutate(followup = "> 5 years", estimand = "sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  mutate(outcome="Cancer death",
         followup = factor(followup, levels = rev(c("<= 3 months", "3 months to 5 years", "> 5 years"))))

# Cancer incidence 
short <- read.csv(paste0(path_short, "short_cancer_adj_rbst.csv"))
med <- read.csv(paste0(path_med_long, "med_cancer.csv")) 
long <- read.csv(paste0(path_med_long, "long_cancer.csv")) 

# Merge results into 1 dataframe (only "case" term)
all_cancer <- short |> filter(grepl("case.*", term)) |> mutate(followup = "<= 3 months", estimand="RR") |> select(followup, gender, term, estimand, estimate=RR, lower, upper, p) |>
  bind_rows(med |> mutate(followup = "3 months to 5 years", estimand="sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  bind_rows(long |>  mutate(followup = "> 5 years", estimand="sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  mutate(outcome="Cancer incidence",
         followup = factor(followup, levels = rev(c("<= 3 months", "3 months to 5 years", "> 5 years")))) |>
  arrange(desc(followup))

# combined death and incidence plot
all_death_cancer <- all_death |> bind_rows(all_cancer) |>
  mutate(outcome = factor(outcome, levels = rev(c("Cancer death", "Cancer incidence"))),
         followup = factor(followup, c("<= 3 months", "3 months to 5 years", "> 5 years"))) |>
  arrange(outcome)


# Plot all by time period as facet and all outcomes on y-axis
library(grid); library(gtable)

all_plt <- ggplot() +
  geom_errorbarh(aes(xmin=log(lower), xmax=log(upper), y = outcome, colour=gender), 
                 data = all_death_cancer, height=.3, size=.3,  position=position_dodge(width=.3)) +
  geom_point(aes(x=log(estimate), y = outcome, colour=gender), 
             data = all_death_cancer, size = .6, position=position_dodge(width=.3))  +
  scale_colour_manual(values = c("#7A28CB", "#6F9B1C"), name = "Sex", na.translate = F) +
  geom_vline(xintercept = log(1), linetype="dashed", colour="black", linewidth=.2) +
  facet_grid(~followup, switch = "y",  drop = T, scales="free_x") +
  scale_x_continuous(breaks = log(2^(c(-2,-1,0,1,2,3,4))), 
                     labels = 2^(c(-2,-1,0,1,2,3,4)),
                     minor_breaks = NULL) +
  xlab(paste(c(rep(" ", 15), "Risk ratio", rep(" ", 35), "Subdistribution hazard ratio", rep(" ", 30)), collapse="")) +
  theme_bw()+
  theme(strip.text = element_text(size = 6, colour = "black", margin=margin(.5,.5,0.5,.5, "mm")),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.y = element_blank(),
        panel.border= element_rect(linewidth=.4, colour = "black"),
        panel.grid.major = element_line(linewidth=.3, colour = "#d3d3d3"),
        strip.background = element_rect(fill = "#d3d3d3", color = "#d3d3d3"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 7, colour = "black", vjust=.5),
        legend.position = "bottom")
gt <- ggplot_gtable(ggplot_build(all_plt))
gt$widths[6] <- 3*gt$widths[6]
grid.draw(gt)
do.call(ggsave, c(list(plot = gt, filename = paste0(path_save, "graphical_abstract1.jpeg")), list(
  device = "jpeg",
  unit = "mm",
  width = 120,
  height = 60,
  dpi = 1000
)))

# Relative risks ----
# Proportion of cancer diagnosed within AF and non-AF patients across different follow-up periods
# < 3m: proportion of all AF patients diagnosed with cancer
outcome_desc <- read.csv("./results/outcome_desc_timesplit.csv")
summary_df <- outcome_desc |>  select(!c(starts_with("summary"))) |>
  filter(var == "Cancer Incidence", follow == "<= 3 months") |>
  mutate(var = case_when(
    var == "Cancer Incidence" ~ "Cancer incidence analysis"
  )) |>
  mutate(Case = (n_Women_Case+n_Men_Case)/(N_Women_Case+N_Men_Case)*100,
         Control = (n_Women_Control+n_Men_Control)/(N_Women_Control+N_Men_Control)*100) |>
  select(var, Case, Control) |>
  pivot_longer(!var, names_to = "Case", values_to = "prop") |>
  mutate(Case = case_when(
    Case == "Case" ~ "New-onset atrial fibrillation",
    Case == "Control" ~ "Age-sex matched controls"
  )) |>
  mutate(Case = factor(Case, levels = rev(c("New-onset atrial fibrillation", "Age-sex matched controls")))) |>
  arrange(Case)

p2 <- ggplot() +
  geom_bar(aes(y = Case, x =prop), data = summary_df, 
           stat = "identity", position=position_dodge(), width = .8, fill = outcome_colours$outcome, alpha=.7) +
  geom_text(aes(label = sprintf("%g%%", round(prop, 1)), x = prop, y = Case), 
            data = summary_df, hjust=-0.1, size = 5/.pt) +
  theme_minimal() +
  ylab("")+
  xlab("Cancer incidence (%)") +
  scale_x_continuous(expand=expansion(mult=c(0,0.2)))+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6, colour = "black", hjust=0),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y.left  = element_line(colour = "#F5F5F5", linewidth = .5),
        plot.margin = grid::unit(c(0,0,0,0), "mm")
        ,
        aspect.ratio = .5/1
  )


do.call(ggsave, c(list(plot = p2, filename = paste0(path_save, "graphical_abstract2.jpeg")), list(device="jpeg",
                                                                                                  unit="cm",
                                                                                                  width = 5.6,
                                                                                                  height = 1.83,
                                                                                                  dpi = 1000)))

# Cancer cumulative incidence -----
rm(list=ls());gc()
library(lubridate)
library(purrr)
library(forcats)
library(survival)
library(tidycmprsk)

# load functions, global variables and codelists
source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")

path_save <- "./tables_figs/stata/"

convert_cif_to_df <- function(cif, group_vars) {
  single_cif_to_df <- function(cif) {
    # change cuminc object into a dataframe
    cif$tidy
  }
  # change nested dataframe with cuminc objects into an unnested dataframe
  t <- cif |> 
    group_by(!!!syms(group_vars)) |>
    mutate(cif = map(cumincs, .f = ~single_cif_to_df(.x))) |>
    select(!!!syms(group_vars), cif) |>
    unnest(cif) |>
    mutate(gender = recode_gender(gender),
           case = strata,
           outcome_mod = case_when(
             # Cancer outcome
             outcome == "cancer" ~ "Cancer",
             outcome == "death" ~ "Any cause death",
             # Censoring,
             outcome == "end" ~ "Censored",
             TRUE ~ outcome
           ),
           outcome = factor(outcome, levels = c(
             "Cancer", "Any cause death",
             "Censored"
           )))
  return(t)
}

cif_cancer <- readRDS("./data/final_data/cif_cancer_withCI.rds")

cif_cancer_df <- convert_cif_to_df(cif_cancer, group_vars="gender")

write.csv(cif_cancer_df, "./results/cif_cancer_withCI.csv", row.names = F)

cif_cancer_df <- read.csv("./results/cif_cancer_withCI.csv") |>
  mutate(case = factor(case, levels = c("Case", "Control")),
         outcome = factor(outcome, levels = c("Cancer", "Any cause death", "Censored")))

max_time <- 12*20

p3 <- cif_cancer_df |>
  mutate(time_shifted = time + time_shift) |>
  ggplot() +
  # shift time 3 months back because shifted the start date by 3 months
  geom_ribbon(aes(x = time_shifted, ymin=100*conf.low, ymax=100*conf.high, fill=outcome, linetype=case), alpha = .2) +
  geom_line(aes(x=time_shifted, y = 100*estimate, col=outcome, linetype=case), linewidth=.2) + 
  scale_x_continuous(
    limits = c(0, max_time),
    minor_breaks = 12*c(seq(5, max_time, 5)),
    breaks = 12*c(0.25, 10, 20),
    labels = c("3m", 10, 20),
    expand = c(0,0.1)
  ) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "Case") +
  scale_colour_manual(values = c(outcome_colours$outcome, outcome_colours$compete), name = "Event type") +
  scale_fill_manual(values = c(outcome_colours$outcome, outcome_colours$compete), name = "Event type") +
  scale_y_continuous(expand = expansion(mult=c(0, 0.1))) +
  guides(linetype=guide_legend(override.aes=list(fill=NA))) +
  facet_grid(~gender,labeller = label_wrap_gen(width = 15)) +
  labs(x = "Time since start of follow-up (years)",
       y = "Cumulative Incidence (%)",
       col = "Event type", 
       fill = "Event type",
       linetype = "Case") +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(),
        legend.key = element_rect(fill = "white"),
        legend.spacing.x= unit(1, "mm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10.5),
        panel.background = element_rect(fill="white"),
        panel.grid.major = element_line(linewidth=.3, colour = "#d3d3d3"),
        panel.grid.minor = element_line(linewidth=.1, colour = "#d3d3d3"),
        panel.spacing.x = unit(2, "lines"),
        axis.line = element_line(colour = "black", linewidth=.1),
        strip.background  = element_rect(fill = "#d3d3d3", color = "#d3d3d3"),
        strip.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(vjust = -1),
        axis.ticks.x = element_line(linewidth=.1),
        axis.ticks.length.x = unit(.5, "mm"),
        axis.ticks.y = element_blank(),axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 11))


do.call(ggsave, c(list(plot = p3, filename = paste0(path_save, "graphical_abstract3.jpeg")), list(device = "jpeg",
                                                                                                  unit = "cm",
                                                                                                  width = 12,
                                                                                                  height = 10,
                                                                                                  dpi=1000)))


