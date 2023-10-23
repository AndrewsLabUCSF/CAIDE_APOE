setwd("/Users/elinorvelasquez/Desktop/HABSHD/v3")
library(tidyverse)
library(dplyr)
library(broom)
# library(fabricatr)
library(gtsummary)
library(gt)
# library(emmeans)
# theme_gtsummary_compact()

aa_raw <- read_csv("HD 1 African American 50+ Request 173.csv", 
            na = c("", "NA", "9999")) %>% janitor::clean_names() 
ma_raw <- read_csv("HD 1 Mexican American 50+ Request 173.csv", 
            na = c("", "NA", "9999")) %>% janitor::clean_names()
nhw_raw <- read_csv("HD 1 Non-Hispanic White 50+ Request 173.csv", 
                    na = c("", "NA", "9999")) %>% janitor::clean_names()

habshd.raw <- rbind(aa_raw, ma_raw, nhw_raw) %>%
  mutate(
    id_race_white = as.factor(id_race_white),
    id_race_black = as.factor(id_race_black),
    id_race_indian_alaska = as.factor(id_race_indian_alaska),
    id_race_asian = as.factor(id_race_asian),
    id_race_japanese = as.factor(id_race_japanese),
    id_race_korean = as.factor(id_race_korean),
    id_race_vietnamese = as.factor(id_race_vietnamese),
    id_race_native_hawaiian = as.factor(id_race_native_hawaiian),
    id_race_guam_chamorro = as.factor(id_race_guam_chamorro),
    id_race_samoan = as.factor(id_race_samoan),
    id_race_other_pacific = as.factor(id_race_other_pacific),
    id_race_other = as.factor(id_race_other),
    id_hispanic = as.factor(id_hispanic),
    id_hispanic_other = as.factor(id_hispanic_other)
  )

#### Construct sdoh_score ####

sdoh <- habshd.raw %>%  
  dplyr::select(id_residence, id_education, id_marital_status, id_income, 
         id_retire, id_retire_years, id_job, id_job_retire, 
         health_status, starts_with("hc"), 
         pswq_total, rapa_1_total, rapa_2_total,
         chronic_stress_total, 
         amnart_z_score, 
         adi_nat_rank, adi_state_rank, 
         starts_with('eds'), starts_with('chronic_stress'),
         starts_with('social_support')                          
         ) %>%
  mutate(
    occupation = case_when(
      id_retire == 1 ~ 0, # retired
      id_job %in% c(
          "Currently Disabled - was a Truck Driver",
          "Currently disabled - previous job: Custodian",
          "Currently on disability",
          "Currently on disability - was a painter",
          "Disability",
          "Disability (Construction)",
          "Disability (Formerly: Home healthcare)",
          "Disability (Formerly: Housekeeper)",
          "Disability (Formerly: McDonalds worker)",
          "Disability (Formerly: clinical staff)",
          "Disability (cashier)",
          "Disability (sales)",
          "Disability (waitress)",
          "Disability, primary occupation used to be in pipeline/drilling",
          "Disable 2016",
          "Disabled",
          "Disabled - Assembler",
          "Disabled - Factory Packer",
          "Disabled - Factory Worker",
          "Disabled - Housekeeper",
          "Disabled - Laborer",
          "Disabled - Leather cuttter",
          "Disabled - was a housekeeper",
          "Disabled 2012 - was a laborer",
          "Disabled(formerly: Supervisor)",
          "Disabled, previous occupation was housekeeping",
          "Disabled; worked as a custodian before that",
          "Is on disablitiy",
          "Driver, landscaping, currently on disability",
          "Electrician (disabled)",
          "Mechanic - currently on disability",
          "Medically disabled previous business owner",
          "On Disability",
          "On disability",
          "On disability since 2005",
          "Home Health, currently on disability",
          "Was on disability",
          "Was working in the collection department before filling disability",
          "disability",
          "disabled",
          "on disability",
          "Administrative Assistant, currently on disability",
          "\"Im work injured.\""
          ) ~ 0, # disability
      id_job %in% c("Cafeteria worker/housewife",
                    "Home Maker",
                    "Homekeeper",
                    "Homemaker",
                    "Homemaker/Housewife",
                    "House wife",
                    "Housewife",
                    "STAY AT HOME MOM",
                    "Stay at Home Mom",
                    "Stay at home Mom",
                    "ama de casa/housewife",
                    "house wife",
                    "housewife",
                    "housewife & volunteer",
                    "home",
                    "home maker",
                    "homemaker",
                    "house (husband pension)"
                    ) ~ 0, # homemaker
      id_job %in% c( 
                    "Cannot find job Caregiver", 
                    "Currently unemployed (general manager for a franchise restaurant)",
                    "Currently unemployed but was senior pc maintenance tech",
                    "Currently unemployed, but worked as a custodian",
                    "Curretly not working",
                    "Does not work",
                    "Hair Stylist (No longer working, but not retired)",
                    "In the process of retiring",
                    "Looking for work",
                    "Never worked",
                    "No",
                    "No longer working (takes care of son)",
                    "No longer works, was a hair stylist",
                    "Not Currently Working",
                    "Not working",
                    "She is currently unemployed.",
                    "Takes care of her sister but before worked in the assembly line",
                    "Unemployed",
                    "Unemployed but completes small side jobs",
                    "Waiting to file for retirement, but no longer working.",
                    "house wife (unemployed at the moment)",
                    "not currently working",
                    "not employed currently",
                    "not working right now",
                    "not working, but hasn't retired",
                    "taking care of son, not currently working",
                    "unemployed",
                    "unemployed at the moment (caregiver before)",
                    "unemployed, but normally Certified nurse assistant",
                    "unemployed, previously Patient care attendant",
                    "works part time but not at the moment"
                    ) ~ 1,
      id_retire == 0 & id_job != "NA" ~ 0, # employed
      TRUE ~ NA_real_ 
    ),
    income = case_when(
      id_income > 50000 ~ 0, 
      id_income <= 50000 ~ 1,
      TRUE ~ NA_real_
    ),
    residence = case_when(
      id_residence == 1 | id_residence == 3 ~ 0, # own, live rent-free
      id_residence == 2 | id_residence == 4 ~ 1, # rent, don't know
      id_residence == 8888 ~ NA_real_, # refused
      TRUE ~ NA_real_
    ),
    education = case_when(
      id_education < 16 ~ 1, 
      id_education >= 16 ~ 0,
      TRUE ~ NA_real_
    ),
    marital_status = case_when(
      id_marital_status == 1 ~ 0, # married
      id_marital_status == 2 | # divorced
      id_marital_status == 3 | # separated
      id_marital_status == 4 | # widowed
      id_marital_status == 5 ~ 1, # never married
      id_marital_status == 6 ~ NA_real_, # no information
      id_marital_status == 8888 ~ NA_real_, # refused
      TRUE ~ NA_real_
    ),
    # need to fix
  ) %>%
  rowwise() %>%
  mutate(
    test = sum(eds_1, eds_2, eds_3, eds_4, eds_5, eds_6, 
               eds_7, eds_8, eds_9, na.rm = TRUE),
    discrim_bucket = ntile(test, 3)
    # discrim_bucket     n
    # 1                  2715
    # 2                  0
    # 3                  0
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    discrim_buckets =  ntile(sum(eds_1, eds_2, eds_3, eds_4, eds_5, eds_6, 
                              eds_7, eds_8, eds_9, na.rm = TRUE), 3), # tertiles
    discrim = case_when(
          discrim_buckets == 1 ~ 0,
          discrim_buckets == 2 ~ 0,
          discrim_buckets == 3 ~ 1,
          TRUE ~ NA_real_
        ),
    cs_1b = case_when(
      chronic_stress_1b == 1 ~ 0,
      chronic_stress_1b == 2 | chronic_stress_1b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_2b = case_when(
      chronic_stress_2b == 1 ~ 0,
      chronic_stress_2b == 2 | chronic_stress_2b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_3b = case_when(
      chronic_stress_3b == 1 ~ 0,
      chronic_stress_3b == 2 | chronic_stress_3b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_4b = case_when(
      chronic_stress_4b == 1 ~ 0,
      chronic_stress_4b == 2 | chronic_stress_4b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_5b = case_when(
      chronic_stress_5b == 1 ~ 0,
      chronic_stress_5b == 2 | chronic_stress_5b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_6b = case_when(
      chronic_stress_6b == 1 ~ 0,
      chronic_stress_6b == 2 | chronic_stress_6b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_7b = case_when(
      chronic_stress_7b == 1 ~ 0,
      chronic_stress_7b == 2 | chronic_stress_7b == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_8c = case_when(
      chronic_stress_8c == 1 ~ 0,
      chronic_stress_8c == 2 | chronic_stress_8c == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    cs_total = sum(cs_1b, cs_2b, cs_3b, cs_4b, cs_5b, cs_6b, cs_7b, cs_8c, 
                   na.rm = TRUE),
    chronic_stress = case_when(
      cs_total < 2 ~ 0,
      cs_total >= 2 ~ 1,
      TRUE ~ NA_real_
    ),
    soc_supp_buckets =  ntile(sum(social_support_1, social_support_2, 
                                  social_support_3, social_support_4, 
                                  social_support_5, social_support_6, 
                                  social_support_7, social_support_8, 
                                  social_support_9, social_support_10,
                                  social_support_11, social_support_12, 
                                  na.rm = TRUE), 3), # tertiles
    soc_support = case_when(
      soc_supp_buckets == 1 ~ 0,
      soc_supp_buckets == 2 ~ 0,
      soc_supp_buckets == 3 ~ 1,
      TRUE ~ NA_real_
    ),
    area_development_index = case_when(
    0 <= adi_nat_rank & adi_nat_rank <= 75 ~ 0,
    75 < adi_nat_rank & adi_nat_rank <= 100 ~ 1,
    TRUE ~ NA_real_
    ),
    health_insurance = case_when(
      hc_insurance_no_insurance == 0 ~ 0, # No
      hc_insurance_no_insurance == 1 ~ 1, # Yes
      TRUE ~ NA_real_
    ),
    personal_care_provider = case_when(
      hc_pcp == 0 ~ 0,  # No
      hc_pcp > 0 ~ 1, # Yes
      TRUE ~ NA_real_
    ),
    health_care_cost = case_when(
      hc_12month == 0 ~ 0,  # No
      hc_12month == 1 ~ 1, # Yes
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    sdoh_score = sum(occupation, income, residence, education, marital_status, 
                     discrim, chronic_stress, soc_support, 
                     area_development_index, health_insurance, 
                     personal_care_provider, health_care_cost, na.rm = TRUE
    )
  ) %>%
  ungroup()

sdoh %>% count(discrim_buckets)
print(sdoh %>% count(sdoh_score), n=50)
## problematic variables:
sdoh %>% count(discrim)
sdoh %>% count(chronic_stress)
sdoh %>% count(soc_support)

#### Make table ####
sdoh %>%
  dplyr::select(
    occupation, income, residence, education, marital_status,
    discrim, chronic_stress, soc_support, area_development_index,
    health_insurance, personal_care_provider, health_care_cost, sdoh_score
  ) %>%
  tbl_summary(
    # by = health_insurance,
    type = list(occupation ~ "categorical", 
                income ~ "categorical",
                residence ~ "categorical", 
                education ~ "categorical",
                marital_status ~ "categorical", 
                discrim ~ "dichotomous", 
                chronic_stress ~ "dichotomous", 
                soc_support ~ "dichotomous", 
                # race ~ "categorical",
                area_development_index ~ "categorical", 
                health_insurance ~ "categorical",
                personal_care_provider ~ "categorical", 
                health_care_cost ~ "categorical",
                sdoh_score ~ "continuous"
    ),
    statistic = list(
      # all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    # digits = all_continuous() ~ 2,
    label = list(occupation ~ "Occupation", 
                 income ~ "Income", 
                 residence ~ "Residence", 
                 # race ~ "Race", 
                 education ~ "Education", 
                 marital_status ~ "Marital Status",
                 discrim ~ "Discrimination", 
                 chronic_stress ~ "Chronic Stress", 
                 soc_support ~ "Social Support", 
                 area_development_index ~ "Area Development Index",
                 health_insurance ~ "Health Insurance", 
                 personal_care_provider ~ "Personal Care Provider",
                 health_care_cost ~ "Health Care Cost",
                 sdoh_score ~ "SDOH Score"
    ),
    missing_text = "Missing"
  ) %>%
  modify_header(label = "**SDOH Score**") %>%
  bold_labels() %>%
  as_gt %>%
  gt::gtsave(
    filename = "~/Desktop/HABSHD/v3/table_habshd.png", 
    vwidth = 648, vheight = 288
  )


