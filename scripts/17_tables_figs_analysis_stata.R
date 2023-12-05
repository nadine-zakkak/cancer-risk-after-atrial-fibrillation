# 7/September/2023
# Forest plots and summary of results

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

forest_plot <- function(data, facets) {
  ggplot() +
    geom_errorbarh(aes(xmin=log(lower), xmax=log(upper), y = followup), data = data, height=.1, size=.2, colour = "#1F1F1F") +
    geom_point(aes(x=log(estimate), y = followup, shape = estimand), data = data, size = 1)  +
    scale_shape_manual(values = rev(c(17, 16)), name = "Estimand") +
    geom_vline(xintercept = log(1), linetype="dashed", colour="#999999")+
    facet_grid(as.formula(facets), switch = "y") +
    scale_x_continuous(limits = c(log(min(data$lower, na.rm = T)), log(max(data$upper, na.rm = T))),
                       breaks = log(2^(c(-2,-1,0,1,2,3,4))), 
                       labels = 2^(c(-2,-1,0,1,2,3,4)),
                       minor_breaks = NULL) +
    theme_bw()+
    theme(strip.text = element_text(size = 7, colour = "black"),
          axis.text.x = element_text(size = 6, colour = "black"),
          axis.text.y = element_text(size = 6, colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.background = element_rect(colour="#d3d3d3", fill = "#d3d3d3"))
}

forest_plot_by_time <- function(data, facets) {
  ggplot() +
    geom_errorbarh(aes(xmin=log(lower), xmax=log(upper), y = outcome, colour=gender), 
                   data = data, height=.1, size=.2,  position=position_dodge(width=.3)) +
    geom_point(aes(x=log(estimate), y = outcome, colour=gender), data = data, size = 1, position=position_dodge(width=.3))  +
    scale_colour_manual(values = c("#7A28CB", "#6F9B1C"), name = "Sex", na.translate = F) +
    geom_vline(xintercept = log(1), linetype="dashed", colour="#999999")+
    facet_grid(as.formula(facets), switch = "y", scales = "free", drop = T) +
    scale_x_continuous(breaks = log(2^(c(-2,-1,0,1,2,3,4))), 
                       labels = 2^(c(-2,-1,0,1,2,3,4)),
                       minor_breaks = NULL) +
    theme_bw()+
    theme(strip.text = element_text(size = 8, colour = "black"),
          axis.text.x = element_text(size = 8, colour = "black"),
          axis.text.y = ggtext::element_markdown(size =  rev(c(rep(8, 2), 9, rep(8, 6), 9, rep(8, 8)))
                                                 , colour = "black"
                                                 , face = rev(c(rep("plain", 2), "bold", rep("plain", 6), "bold", rep("plain", 8)))
          ),
          axis.ticks.y = element_line(colour = c(rev(c(rep("black", 2), "transparent", rep("black", 6), "transparent", rep("black", 8))))),
          axis.title.y = element_blank(),
          strip.background = element_rect(colour="#d3d3d3", fill = "#d3d3d3"))
}

# Cancer death ----
short <- read.csv(paste0(path_short, "short_cancer_death_adj_rbst.csv"))
med <- read.csv(paste0(path_med_long, "med_cancer_death.csv")) 
long <- read.csv(paste0(path_med_long, "long_cancer_death.csv")) 

# Merge results into 1 dataframe (only "case" term)
all_death <- short |> filter(grepl("case.*", term)) |> mutate(followup = "<= 3 months", estimand="RR") |> select(followup, gender, term, estimand, estimate=RR, lower, upper, p) |>
  bind_rows(med |> mutate(followup = "3 months to 5 years", estimand="sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  bind_rows(long |> mutate(followup = "> 5 years", estimand = "sHR") |> select(followup, gender, term, estimand, estimate, lower, upper, p)) |>
  mutate(outcome="Cancer death",
         followup = factor(followup, levels = rev(c("<= 3 months", "3 months to 5 years", "> 5 years"))))

# Cancer incidence ----
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

# combined death and incidence plot -----
all_death_cancer <- all_death |> bind_rows(all_cancer)

death_cancer_plt <- forest_plot(all_death_cancer, "outcome~gender")
do.call(ggsave, c(list(plot = death_cancer_plt, filename = paste0(path_save, "death_cancer_plt.jpeg")), save_image_var_large))
# forest_plot_by_time(all_death_cancer, "gender~followup")

all_death_cancer |>
  mutate(estimate = format_numb(estimate),
         lower = format_numb(lower),
         upper = format_numb(upper),
         p = ifelse(p<0.01, "<0.01", format_numb(p))) |>
  mutate(estimate_ci = paste0(estimate, " (", lower, ", ", upper, ")")) |>
  pivot_wider(id_cols = c(outcome, followup), names_from = gender, values_from = c(estimate_ci, p), names_vary = "slowest") |>
  flextable() |>
  set_header_labels(outcome = "", followup = "Time Period",
                    estimate_ci_Men = "Estimate (95% CI)", p_Men = "p-value",
                    estimate_ci_Women = "Estimate (95% CI)", p_Women = "p-value") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  merge_v(j=1) |>
  hline(i=c(3)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "death_cancer_est.docx"))

# Cancer incidence by region -----
short <- read.csv(paste0(path_short, "short_region_adj_rbst.csv"))
med <- read.csv(paste0(path_med_long, "med_region.csv")) 
long <- read.csv(paste0(path_med_long, "long_region.csv")) 

# Merge results into 1 dataframe (only "case" term)
all_region <- short |> filter(grepl("case.*", term)) |> mutate(followup = "<= 3 months", estimand="RR") |> select(followup, region=region_desc, gender, term, estimand, estimate=RR, lower, upper, p) |>
  bind_rows(med |> filter(grepl("case.*", term)) |> mutate(followup = "3 months to 5 years", estimand="sHR") |> select(followup, region=region_desc, gender, term, estimand, estimate, lower, upper, p)) |>
  bind_rows(long |> filter(grepl("case.*", term)) |> mutate(followup = "> 5 years", estimand="sHR") |> select(followup, region=region_desc, gender, term, estimand, estimate, lower, upper, p)) |>
  mutate(outcome="Cancer incidence by region",
         followup = factor(followup, levels = rev(c("<= 3 months", "3 months to 5 years", "> 5 years"))),
         region = factor(region, levels = region_order)) |>
  arrange(desc(followup))

region_plt <- forest_plot(all_region, "region~gender")
do.call(ggsave, c(list(plot = region_plt, filename = paste0(path_save, "region_plt.jpeg")), save_image_var_large))

all_region |>
  mutate(estimate = format_numb(estimate),
         lower = format_numb(lower),
         upper = format_numb(upper),
         p = ifelse(p<0.01, "<0.01", format_numb(p))) |>
  mutate(estimate_ci = paste0(estimate, " (", lower, ", ", upper, ")")) |>
  pivot_wider(id_cols = c(region, followup), names_from = gender, values_from = c(estimate_ci, p), names_vary = "slowest") |>
  arrange(region) |>
  flextable() |>
  set_header_labels(region = "", followup = "Time Period",
                    estimate_ci_Men = "Estimate (95% CI)", p_Men = "p-value",
                    estimate_ci_Women = "Estimate (95% CI)", p_Women = "p-value") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  merge_v(j=1) |>
  hline(i=c(3,6,9,12,15)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "region_est.docx"))

# Cancer incidence by organ system -----
short <- read.csv(paste0(path_short, "short_organsys_adj_rbst.csv"))
med <- read.csv(paste0(path_med_long, "med_organsys.csv")) 
long <- read.csv(paste0(path_med_long, "long_organsys.csv")) 

# Merge results into 1 dataframe (only "case" term)
all_organsys <- short |> filter(grepl("case.*", term)) |> mutate(followup = "<= 3 months", estimand="RR") |> select(followup, organsys=organsys_desc, gender, term, estimand, estimate=RR, lower, upper, p) |>
  bind_rows(med |> filter(grepl("case.*", term)) |> mutate(followup = "3 months to 5 years", estimand="sHR") |> select(followup, organsys=organsys_desc, gender, term, estimand, estimate, lower, upper, p)) |>
  bind_rows(long |> filter(grepl("case.*", term)) |> mutate(followup = "> 5 years", estimand="sHR") |> select(followup, organsys=organsys_desc, gender, term, estimand, estimate, lower, upper, p)) |>
  mutate(outcome="Cancer incidence by organ system",
         followup = factor(followup, levels = rev(c("<= 3 months", "3 months to 5 years", "> 5 years"))),
         organsys = factor(organsys, levels = organsys_order)) |>
  arrange(desc(followup))

organsys_plt <- forest_plot(all_organsys, "organsys~gender")
do.call(ggsave, c(list(plot = organsys_plt, filename = paste0(path_save, "organsys_plt.jpeg")), save_image_var_large))

all_organsys |>
  mutate(estimate = ifelse(is.na(estimate), "-", format_numb(estimate)),
         lower = ifelse(is.na(lower), "-", format_numb(lower)),
         upper = ifelse(is.na(upper), "-", format_numb(upper)),
         p = ifelse(p<0.01, "<0.01", format_numb(p))) |>
  mutate(estimate_ci = paste0(estimate, " (", lower, ", ", upper, ")")) |>
  pivot_wider(id_cols = c(organsys, followup), names_from = gender, values_from = c(estimate_ci, p), names_vary = "slowest") |>
  arrange(organsys) |>
  flextable() |>
  set_header_labels(organsys = "", followup = "Time Period",
                    estimate_ci_Men = "Estimate (95% CI)", p_Men = "p-value",
                    estimate_ci_Women = "Estimate (95% CI)", p_Women = "p-value") |>
  add_header_row(values = c("", "Men", "Women"), colwidths = c(2,2,2)) |>
  merge_v(j=1) |>
  hline(i=c(3,6,9,12,15,18,21)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "organsys_est.docx"))

# Table of all estimates -----
all_estimates_df <- all_death_cancer |> mutate(value=outcome) |>
  bind_rows(all_region |> select(-outcome, value=region) |> mutate(outcome="Body Region")) |>
  bind_rows(all_organsys |> select(-outcome, value=organsys) |> mutate(outcome = "Organ System",
                                                                       value = paste0(" ", value))) |>
  mutate(value = factor(value, levels = c("Cancer death", "Cancer incidence", 
                                          region_order,
                                          paste0(" ", organsys_order)
  )),
  followup = factor(followup, c("<= 3 months", "3 months to 5 years", "> 5 years"))) |>
  arrange(value) |>
  mutate(estimate = ifelse(is.na(estimate), "-", format_numb(estimate)),
         lower = ifelse(is.na(lower), "-", format_numb(lower)),
         upper = ifelse(is.na(upper), "-", format_numb(upper)),
         p = ifelse(p<0.01, "<0.01", format_numb(p))) |>
  mutate(estimate_ci = paste0(estimate, " (", lower, ", ", upper, ")")) |>
  pivot_wider(id_cols = c(gender, outcome, value), names_from = c(followup), values_from = c(estimate_ci, p), names_vary = "slowest") |>
  arrange(gender, value) |>
  flextable() |>
  merge_v(j=c(1,2)) |>
  set_header_labels(gender = "", outcome="", value="", followup = "Time Period",
                    `estimate_ci_<= 3 months` = "Estimate (95% CI)", `p_<= 3 months` = "p-value",
                    `estimate_ci_3 months to 5 years` = "Estimate (95% CI)", `p_3 months to 5 years` = "p-value",
                    `estimate_ci_> 5 years` = "Estimate (95% CI)", `p_> 5 years` = "p-value") |>
  add_header_row(values = c("", "<= 3 months", "3 months to 5 years", "> 5 years"), colwidths = c(3,2,2, 2)) |>
  hline(i=c(1,2,8,15,16,17,23,30)) |>
  add_table_theme()  

all_estimates_df |>
  save_as_docx(path = paste0(path_save, "allmodel_estimates.docx"))

all_estimates_df |>
  save_as_docx(path = "./tables_figs/supplementary/supplementary_table2.docx")

# Plot all by time period as facet and all outcomes on y-axis -----
all <- all_death_cancer |>
  bind_rows(all_region |> select(-outcome, outcome=region)) |>
  bind_rows(all_organsys |> select(-outcome, outcome=organsys) |>
              mutate(outcome = paste0(" ", outcome))) |>
  bind_rows(data.frame(followup = rep(c("<= 3 months", "3 months to 5 years",  "> 5 years"), 2),
                       outcome = c(rep("Anatomical Region", 3), rep("System", 3)),
                       estimand = rep(c("RR", "sHR", "sHR"), 2))) |>
  mutate(outcome = factor(outcome, levels = rev(c("Cancer death", "Cancer incidence", 
                                                  "Anatomical Region",
                                                  region_order,
                                                  "System",
                                                  paste0(" ", organsys_order)
  ))),
  followup = factor(followup, c("<= 3 months", "3 months to 5 years", "> 5 years"))) |>
  arrange(outcome)

# Add space between (1) short and (2) med and long facets 
library(grid); library(gtable)
all_plt <- forest_plot_by_time(all, "~followup") + 
  xlab(paste(c(rep(" ", 8), "RR", rep(" ", 60), "sHR", rep(" ", 27)), collapse="")) +
  theme(axis.title.x = element_text(size = 10, colour = "black"))
gt <- ggplot_gtable(ggplot_build(all_plt))
gt$widths[6] <- 5*gt$widths[6]
grid.draw(gt)
do.call(ggsave, c(list(plot = gt, filename = paste0(path_save, "all_bytime_plt_extraspace.jpeg")), save_image_var_large))
