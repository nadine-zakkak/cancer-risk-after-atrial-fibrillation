# Order of body regions
region_order <- c(
  "Brain/CNS"
  ,"Head and neck"
  ,"Chest"
  ,"Abdomen"
  ,"Pelvis"
  ,"Upper and lower limbs"
  ,"Other"
  ,"Unknown primary"
)

# Order of organ systems
organsys_order <- c(
  "Breast",
  "Cardiovascular",
  "Digestive",
  "Endocrine",
  "Immune & haematological",
  "Nervous",
  "Reproductive",
  "Respiratory",
  "Skeletal",
  "Urinary",
  "Other",
  "Unknown primary"
)

# Codes for sex
gender_codes <- c(Men=0, Women=1)

# Colours for cumulative incidence plots
outcome_colours <- list("outcome" = "#d95f02", "compete"="#1b9e77")

# order of alcohol status categories
alc_levels <- c("Non", "Ex", "Current", "Missing")

# Common image size
save_image_var_small <- list(device = "jpeg", unit = "mm", width = 140, height = 120, dpi=1000)
save_image_var_large <- list(device = "jpeg", unit = "mm", width = 190, height = 150, dpi=1000)

# Time shift for survival analysis
time_shift <- 12*.25

