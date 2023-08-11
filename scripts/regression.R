# Filename: tbl_regression.R
setwd("~/Desktop/Chi_data3/data")
library("tidyverse")
library("gtsummary")
library("broom")
library("MASS")

# upload data:
iidp_mcaide <- read_csv("mcaide_v2_last.csv")
iidp_caide <- read_csv("caide_v2_last.csv")

# prepare data:
iidp_mcaide$diagnosis <- as.factor(iidp_mcaide$diagnosis)
iidp_mcaide$mcaide_apoe <- as.factor(iidp_mcaide$mcaide_apoe)
iidp_caide$diagnosis <- as.factor(iidp_caide$diagnosis)
iidp_caide$caide_apoe <- as.factor(iidp_caide$caide_apoe)
iidp_caide$physical_act <- as.factor(iidp_caide$physical_act)
iidp_mcaide$physical_act <- as.factor(iidp_mcaide$physical_act)

# diagnosis is binary:
iidp_caide_binary <- iidp_caide %>%
  mutate(
    diag = case_when(
      diagnosis == 'N' ~ 0,
      diagnosis == 'CI' ~ 1,
      diagnosis == 'D' ~ 1,
    )
  )
iidp_mcaide_binary <- iidp_mcaide %>%
  mutate(
    diag = case_when(
      diagnosis == 'N' ~ 0,
      diagnosis == 'CI' ~ 1,
      diagnosis == 'D' ~ 1,
    )
  )

# Set reference:
iidp_caide_binary$apoe <- factor(iidp_caide_binary$apoe, ordered = FALSE )
iidp_caide_binary$apoe = relevel(iidp_caide_binary$apoe, ref = "e3/e3")
iidp_mcaide_binary$apoe <- factor(iidp_mcaide_binary$apoe, ordered = FALSE )
iidp_mcaide_binary$apoe = relevel(iidp_mcaide_binary$apoe, ref = "e3/e3")
iidp_caide_binary$gender <- factor(iidp_caide_binary$gender, ordered = FALSE )
iidp_mcaide_binary$gender <- factor(iidp_mcaide_binary$gender, ordered = FALSE )

# build logistic regression models with all the caide variables + apoe:
# model_log <- glm(outcome ~ predictor1 + predictor2, data, family = binomial)
# output <- tbl_regression(model_log, exponentiate = TRUE)

model_caide_logistic_all <- glm(diag ~ i_bmi + i_systbp + physical_act + i_age 
                                + gender + i_education + apoe, 
                                data=iidp_caide_binary, family = binomial)

add_glance_source_note(tbl_regression(model_caide_logistic_all, 
                                      exponentiate = TRUE)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 5. Logistic Regression for CAIDE 
  with All Variables") %>%
  gt::gtsave(
    filename = "model_caide_logistic_all.html"
  )

model_mcaide_logistic_all <- glm(diag ~ i_bmi + i_systbp + physical_act + 
                                   i_age + gender + i_education + apoe, 
                                 data=iidp_mcaide_binary, family = binomial)

add_glance_source_note(tbl_regression(model_mcaide_logistic_all, 
                                      exponentiate = TRUE)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 6. Logistic Regression for MCAIDE 
  with All Variables") %>%
  gt::gtsave(
    filename = "model_mcaide_logistic_all.html"
  )

# build linear regression model using all m/caide variables:
model_caide_linear_all <- lm(cognitive ~ i_bmi + i_systbp + physical_act + 
                               i_age + gender + i_education + apoe, 
                             data=iidp_caide_binary)

add_glance_source_note(tbl_regression(model_caide_linear_all)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 7. Linear Regression for CAIDE 
  with All Variables") %>%
  gt::gtsave(
    filename = "model_caide_linear_all.html"
  )

model_mcaide_linear_all <- lm(cognitive ~ i_bmi + i_systbp + physical_act 
                              + i_age + gender + i_education + apoe, 
                              data=iidp_mcaide_binary)

add_glance_source_note(tbl_regression(model_mcaide_linear_all)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 8. Linear Regression for MCAIDE 
  with All Variables") %>%
  gt::gtsave(
    filename = "model_mcaide_linear_all.html"
  )

# logistic/linear regression using the interaction variable caide_apoe:
# caide:

# remove missing values:
caide_binary <- subset(iidp_caide_binary, caide_apoe != "NA, e3/e3")
caide_binary_NoNAs <- subset(caide_binary, caide_apoe != "NA, e4+")

caide_binary_NoNAs %>% count(caide_apoe)

# Set reference:
caide_binary_NoNAs$caide_apoe <- factor(caide_binary_NoNAs$caide_apoe, 
                                       ordered = FALSE )
caide_binary_NoNAs$caide_apoe = relevel(caide_binary_NoNAs$caide_apoe, 
                                       ref = "Mid, e3/e3")

model_log_caide_it <- glm(diag ~ caide_apoe, data=caide_binary_NoNAs, 
                      family = binomial)

model_lin_caide_it <- lm(cognitive ~ caide_apoe, data=caide_binary_NoNAs)

# mcaide:
# remove missing values:
mcaide_binary <- subset(iidp_mcaide_binary, mcaide_apoe != "NA, e3/e3")
mcaide_binary_NoNAs <- subset(mcaide_binary, mcaide_apoe != "NA, e4+")

mcaide_binary_NoNAs %>% count(mcaide_apoe)

# Set reference:
mcaide_binary_NoNAs$mcaide_apoe <- factor(mcaide_binary_NoNAs$mcaide_apoe, 
                                        ordered = FALSE )
mcaide_binary_NoNAs$mcaide_apoe = relevel(mcaide_binary_NoNAs$mcaide_apoe, 
                                        ref = "Mid, e3/e3")

model_log_mcaide_it <- glm(diag ~ mcaide_apoe, data=mcaide_binary_NoNAs, 
                          family = binomial)

model_lin_mcaide_it <- lm(cognitive ~ mcaide_apoe, data=mcaide_binary_NoNAs)

add_glance_source_note(tbl_regression(model_lin_caide_it)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 1. Linear Regression for CAIDE-APOE 
                 Interaction Term") %>%
  gt::gtsave(
  filename = "caide_lin_interact_term.html"
  )

add_glance_source_note(tbl_regression(model_log_caide_it)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 2. Logistic Regression for CAIDE-APOE 
                 Interaction Term") %>%
  gt::gtsave(
    filename = "caide_log_interact_term.html"
  )

add_glance_source_note(tbl_regression(model_lin_mcaide_it)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 3. Linear Regression for MCAIDE-APOE 
                 Interaction Term") %>%
  gt::gtsave(
    filename = "mcaide_lin_interact_term.html"
  )

add_glance_source_note(tbl_regression(model_log_mcaide_it)) %>% 
  as_gt %>%
  gt::tab_header(title = "Table 4. Logistic Regression for MCAIDE-APOE 
                 Interaction Term") %>%
  gt::gtsave(
    filename = "mcaide_log_interact_term.html"
  )

