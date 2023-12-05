# Convert analysis data into csv files to be read into stata
# 
library(dplyr)

convert_rds_to_csv <- function(filename){
  df <- readRDS(sprintf("./data/final_data/%s.RDS", filename))
  
  stata_data <- df |>
    mutate(male = gender == 0) |>
    filter(across(any_of("excl_analysis"), ~.x==FALSE)) |>
    select(any_of(c("male", "case", "age_start", "diab", "ht", "smok","alc_desc", "follow", "outcome_num", "cancer_bodyregion", "cancer_organsys"))) |>
    mutate(male = as.integer(male),
           case = as.integer(case),
           diab = as.integer(diab),
           ht = as.integer(ht),
           smok = as.integer(smok),
           outcome_num = case_when(
             outcome_num == "end" ~ 0,
             outcome_num == "event" ~ 1,
             outcome_num == "compete" ~ 2
           )) |>
    rename(any_of(c(subtype = "cancer_bodyregion", subtype = "cancer_organsys")))
  
  write.csv(stata_data, sprintf("./data/final_data/stata/%s.csv", filename), row.names = F)
}

tryCatch(expr = {
  for (outcome in c("cancer_diag", "cancer_death", "cancer_region", "cancer_organsys")) {
    for (fup in c("med_term", "long_term")){
      convert_rds_to_csv(sprintf("%s_%s", outcome, fup))
    }
  }
},
error = function(e){
  e <- gsub("[\r\n]", "", e)
  message(sprintf("Error (%s) at %s_%s", e, outcome,  fup))
}
)




