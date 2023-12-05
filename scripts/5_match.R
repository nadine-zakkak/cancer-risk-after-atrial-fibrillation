# matching 1-1
# AF cases with controls 
# based on sex, yob, start of follow up

rm(list=ls());gc()

#load libraries
library(dplyr)
library(lubridate)

progress_file <- "./progress.txt"
write(paste(Sys.time(), "Begin 5b_match.R script"), file = progress_file, append = T)

cases <- readRDS("./data/cases_patients_final.rds")
controls <- readRDS("./data/controls_patients_final.rds")

cases <- cases |> 
  mutate(lower_yob = yob-1,
         upper_yob = yob+1)

cases_match <- cases |> select(epatid, gender, start_case_fu, lower_yob, upper_yob)
controls_match <- controls |> select(epatid, gender, start_control, end_control_fu,  yob)

write(paste(Sys.time(), "data prepared for matching"), file = progress_file, append = T)

# for each AF patients, get set of eligible patients to be matched with 
# randomly choose control for that patient and mark as unavailable for next match
# repeat until all AF patients have had controls

matched_df <- data.frame()
start <- Sys.time()
pb <- txtProgressBar(min = 0, max = nrow(cases_match), initial = 0)
set.seed(1705231653)
for(i in 1:nrow(cases_match)){
  setTxtProgressBar(pb, i)
  matched <- T
  
  lower_yob <- cases_match$lower_yob[i]
  upper_yob <- cases_match$upper_yob[i]
  case_date <- cases_match$start_case_fu[i]
  sex <- cases_match$gender[i]
  case_patid <- cases_match$epatid[i]
  
  eligible_patid <- controls_match |>
    mutate(valid_fu = 
             case_date >= start_control & case_date < end_control_fu) |>
    filter(gender == sex 
           & valid_fu
           & between(yob, lower_yob, upper_yob)) |> #bounds inclusive 
    select(epatid) |> 
    pull()
  
  #get only potential matches that have not already been selected
  possible_patid <- setdiff(eligible_patid, matched_df$control_patid)
  
  if(length(possible_patid)==1){
    select_patid <- possible_patid #only 1 possible patient
  }else{
    select_patid <- sample(possible_patid, size = 1)
  }
  
  if(is.na(select_patid)) matched <- F #didn't find eligible control
  matched_df <- bind_rows(matched_df, data.frame(case_patid = case_patid, #case patid
                                                 matched = matched, 
                                                 control_patid = select_patid,
                                                 start_fu = as.Date(case_date), #start of follow-up for both cases and controls
                                                 numb = length(possible_patid)))
}
end <- Sys.time()
close(pb)
end - start

write(paste(Sys.time(), "Matching complete in:", end-start), file = progress_file, append = T)

saveRDS(matched_df, file = "./data/matched_patid.rds")

write(paste(Sys.time(), "Matched patid saved"), file = progress_file, append = T)

temp_df <- matched_df |>
  inner_join(cases |> select(epatid, yob_af=yob, af_date, end_case),
             by = c("case_patid"="epatid")) |>
  inner_join(controls |> select(epatid, yob_control=yob, start_control, end_control_fu, elig_af, control_af_date=af_date),
             by = c("control_patid"="epatid"))

table(matched_df$matched) #all matched
table(temp_df$elig_af, dnn = "AF later") #12.5% (14,635) developed AF later

match_controls <-
  matched_df |>
  select(epatid = control_patid, start_fu) |>
  inner_join(controls, by = "epatid") |>
  select(epatid, pracid, gender, yob, start_control, 
         start_fu, end_control_fu,
         dor, af_date, cancer_date, elig_af, cancer) |>
  mutate(cancer = ifelse(is.na(cancer), F, cancer),
         age_start = time_length(difftime(start_fu, as.Date(paste0(yob, "-06-15"))), 
                                 unit = "years"))

# check time changes that occurred due to moving start date in controls
controls_timediff <-
  match_controls |>
  select(epatid, start_control, start_fu, end_control_fu) |>
  mutate(fu_orig = time_length(difftime(end_control_fu, start_control), unit = "years"),
         fu_new = time_length(difftime(end_control_fu, start_fu), unit = "years"),
         start_diff = time_length(difftime(start_fu, start_control), unit = "years"))

library(ggplot2)
library(tidyr)
p <- controls_timediff |>
  pivot_longer(cols = c(start_diff, fu_orig, fu_new), names_to = 'Period' , values_to = 'years') |>
  mutate(Period = factor(Period, levels = c('start_diff', 'fu_orig', 'fu_new'))) |>
  select(epatid, Period, years) |>
  ggplot() + geom_boxplot(aes(y = years, fill  = Period)) +
  scale_fill_discrete(labels = c("Difference in start dates", "Original follow up", "New follow up")) +
  theme_bw() + 
  theme(axis.text.x = element_blank()) +
  labs(y = "Time (years)")
ggsave(p, filename = "./results/controls_difftime.png", unit = 'cm', width = 20, height = 15)

match_cases <-
  cases |> select(-c(lower_yob, upper_yob)) |>
  mutate(age_start = time_length(difftime(start_case_fu, as.Date(paste0(yob, "-06-15"))), 
                                 unit = "years"))

# Sense check that matching is reasonable based on age
temp <- 
  matched_df |>
  inner_join(match_cases |> select(epatid, yob_af=yob, start_case_fu, age_case = age_start),
             by = c("case_patid"="epatid")) |>
  inner_join(match_controls |> select(epatid, yob_control=yob, start_control_fu=start_fu, end_control_fu, elig_af, age_control=age_start),
             by = c("control_patid"="epatid"))

p <- ggplot() + geom_point(aes(x=age_case, y = age_control), data = temp) + 
  # add age limits (5 years)
  geom_line(aes(x=age_case, y = age_case-1), data = temp, linetype = "dashed") +
  geom_line(aes(x=age_case, y = age_case+1), data = temp, linetype = "dashed")
ggsave(p, filename = "./results/match_age_scatter.png", unit = 'cm', width = 20, height = 20)

write(paste(Sys.time(), "Matching checks completed"), file = progress_file, append = T)

saveRDS(match_controls, file = "./data/controls.rds")
saveRDS(match_cases, file = "./data/cases.rds")

write(paste(Sys.time(), "Matching files saved"), file = progress_file, append = T)

write(paste(Sys.time(), "End 5b_match.R script"), file = progress_file, append = T)


rm(list=ls()); gc()
source("./scripts/6_clean_files.R")
