# Helper functions ----

#' format_numb
#'
#' @param numb - number to format 
#'
#' @return formatted number to include 3 digits. if <10: 2 decimal places, (10,100): 1 decimal place and >100 0 decimal places
#'
#' @examples: format_numb(0.128) --> 0.13, format_numb(90.4736) --> 90.5, format_numb(123.2737) --> 123
format_numb <- function(numb) {
  
  numb = ifelse(numb<10,
                sprintf(numb, fmt='%.2f'),
                ifelse(numb<100,
                       sprintf(numb, fmt = '%.1f'),
                       sprintf(numb, fmt = '%1.0f')))
}

recode_gender <- function(gender){
  case_when(gender == 0 ~ "Men",
            gender == 1 ~ "Women")
}

recode_case <- function(case){
  ifelse(case, 'Case', 'Control')
}

# requires plyr
add_alc_desc <- function(alc, lookup) {
  plyr::mapvalues(alc,
                  from = lookup$category_num,
                  to = lookup$category_name)
}

num_outcomes <- function(outcome_mod) {
  #change outcomes to factors
  outcome_mod_numeric = factor(case_when(
    # Cancer death outcomes
    outcome_mod =="cancer death" ~ 1,
    outcome_mod =="other cause death" ~ 2, # competing risk
    
    # Cancer outcomes
    outcome_mod=="cancer" ~ 1,
    outcome_mod=="death" ~ 2, # competing risk
    
    # End (censor)
    outcome_mod=="end" ~ 0  
  )
  , levels=c(0, 1, 2),
  labels = c("end", "event", "compete"))
}

# requires plyr and dplyr
add_region_desc <- function(data) {
  data |> mutate(region_desc = factor(plyr::mapvalues(cancer_bodyregion,
                                                from = body_region$body_region_numb,
                                                to = body_region$body_region_desc),
                                      levels = region_order))
}

# requires plyr and dplyr
add_organsys_desc <- function(data) {
  data |> mutate(organsys_desc = factor(plyr::mapvalues(cancer_organsys,
                                                  from = organ_system$organ_sys_numb,
                                                  to = organ_system$organ_sys_desc),
                                        levels = organsys_order))
}

# requires flextable
add_table_theme <- function(table) {
  table |>
    fontsize(part = "all", size = 9) |>
    bold(part = "header") |>
    font(part = "all", fontname =  "Calibri") 
}

# requires tidyr, purrr and broom
tidy_poisson_model <- function(model, group_cols) {
  # Tidy Poisson model to get its estimates in a dataframe along with its confidence intervals, p-value and AIC
  # Input is the model and the group_by variables used before nesting the data
  # Output is a dataframe
  model |>
    mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
           stat = map(glms, glance)) |>
    select(!!!syms(group_cols), estimates, stat) |>
    unnest(cols = c(estimates, stat)) |>
    select(!!!syms(group_cols), term = term, RR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)
}

# requires tidyr, purrr, broom, sandwich and lmtest
tidy_poisson_model_robust <- function(model, group_cols){
  # Tidy Poisson model using robust standard errors 
  # to get its estimates in a dataframe along with its confidence intervals, p-value and AIC
  # Input is the model and the group_by variables used before nesting the data
  # Output is a dataframe
  model |>
    mutate(estimates = map(glms, ~left_join(tidy(coeftest(.x, vcov = sandwich)) #by default uses HC3 .. I think?
                                            , as.data.frame(coefci(.x, vcov = sandwich)) |> 
                                              rownames_to_column("term"), 
                                            by = 'term')),
           stat = map(glms, glance)) |>
    select(!!!syms(group_cols), estimates, stat) |>
    unnest(cols = c(estimates, stat)) |>
    select(!!!syms(group_cols), term = term, RR = estimate, lower = `2.5 %`, upper = `97.5 %`, p = p.value, AIC) |>
    mutate(RR = exp(RR), lower = exp(lower), upper = exp(upper))
}

# requires tidyr, purrr and broom
tidy_fgray_model <- function(data, group_cols) {
  # tidy fine and gray model
  data |>
    mutate(estimates = map(f_gray, ~(tidy(.x, conf.int = T, exponentiate = TRUE))),
           stat = map(f_gray, ~glance(.x))) |>
    select(!!!syms(group_cols), any_of(allcov), estimates, stat) |>
    unnest(cols = estimates) |>
    select(!!!syms(group_cols), term = term, HR = estimate, lower = conf.low, upper = conf.high, p = p.value, stat) |>
    unnest(cols = stat) |>
    mutate(gender = recode_gender(gender)) |>
    select(!!!syms(group_cols), any_of(allcov), term, HR, lower, upper, p, logLik)
}

