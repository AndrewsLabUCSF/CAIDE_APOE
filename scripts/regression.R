setwd("/Users/elinorvelasquez/Desktop/adni_nacc/")
library("tidyverse")
library("dplyr")
library("readxl")
library("readr")
library("ggplot2")
library("gtsummary")
library("broom")
library("MASS")
#library("cowplot")
#library("patchwork")

data_adni_nacc <- read_csv("~/Desktop/adni_nacc/adni_nacc.csv") %>%
  janitor::clean_names()

data_binary <- data_adni_nacc %>%
  mutate(
    diag = case_when(
      dx_bl == 'AD' ~ 1,
      dx_bl == 'ADRD' ~ 1,
      dx_bl == 'CN' ~ 0,
      dx_bl == 'EMCI' ~ 0,
      dx_bl == 'LMCI' ~ 1,
      dx_bl == "MCI" ~ 1,
      dx_bl == 'SMC' ~ 0
    )
  )

# Set reference
data_binary$apoe <- factor(data_binary$apoe, ordered = FALSE )
data_binary$apoe = relevel(data_binary$apoe, ref = "e3/e3")

# 1. logistic model: diagnosis ~ caide + apoe 
model_caide1_apoe <- glm(diag ~ z_caide1 + apoe, 
                            data=data_binary, family = binomial)
model_caide2_apoe <- glm(diag ~ z_caide2 + apoe, 
                         data=data_binary, family = binomial)
model_caide3_apoe <- glm(diag ~ z_caide3 + apoe, 
                         data=data_binary, family = binomial)

# 2. logistic model: diagnosis ~ mcaide + apoe 
model_mcaide1_apoe <- glm(diag ~ z_mcaide1 + apoe, 
                             data=data_binary, family = binomial)
model_mcaide2_apoe <- glm(diag ~ z_mcaide2 + apoe, 
                          data=data_binary, family = binomial)
model_mcaide3_apoe <- glm(diag ~ z_mcaide3 + apoe, 
                          data=data_binary, family = binomial)

# 3 logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education +
#                                   chol + sex + apoe
model_all <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                          data=data_binary, family = binomial)

##########################################################
# Glance
##########################################################
# 1. logistic model: diagnosis ~ caide + apoe  
add_glance_source_note(tbl_regression(model_caide1_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe.png")

add_glance_source_note(tbl_regression(model_caide2_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ CAIDE + APOE (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe.png")

add_glance_source_note(tbl_regression(model_caide3_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ CAIDE + APOE (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe.png")

# 2. logistic model: diagnosis ~ mcaide + apoe  
add_glance_source_note(tbl_regression(model_mcaide1_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ mCAIDE + APOE (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe)) %>% as_gt %>%
  gt::tab_header(title = "Diagnosis ~ mCAIDE + APOE (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe.png")

# 3. logistic model: diagnosis ~ i_bmi + i_systbp + i_age + i_education +
#                                   sex + apoe + cohort
add_glance_source_note(tbl_regression(model_all)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Diagnosis ~ i_bmi + i_systbp + i_age + i_education +
                    sex + apoe + cohort") %>%
  gt::gtsave(filename = "model_all.png")

####################################################
# Regression Stratified by Race 
####################################################
data_binary_nhw <-  dplyr::filter(data_binary, race == "NHW")

data_binary_black <-  dplyr::filter(data_binary, race == "Black")

data_binary_asian <-  dplyr::filter(data_binary, race == "Asian")

data_binary_his <-  dplyr::filter(data_binary, race == "Hispanic")

# 1. logistic model: diagnosis ~ caide + apoe 
model_caide1_apoe_nhw <- glm(diag ~ z_caide1 + apoe, 
                            data=data_binary_nhw, family = binomial)
model_caide2_apoe_nhw <- glm(diag ~ z_caide2 + apoe, 
                             data=data_binary_nhw, family = binomial)
model_caide3_apoe_nhw <- glm(diag ~ z_caide3 + apoe, 
                             data=data_binary_nhw, family = binomial)

# Glance for #1 model
add_glance_source_note(tbl_regression(model_caide1_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_nhw.png")

model_caide1_apoe_black <- glm(diag ~ z_caide1 + apoe, 
                              data=data_binary_black, family = binomial)
model_caide2_apoe_black <- glm(diag ~ z_caide2 + apoe, 
                               data=data_binary_black, family = binomial)
model_caide3_apoe_black <- glm(diag ~ z_caide3 + apoe, 
                               data=data_binary_black, family = binomial)

add_glance_source_note(tbl_regression(model_caide1_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_black.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_black.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_black.png")

model_caide1_apoe_asian <- glm(diag ~ z_caide1 + apoe, 
                               data=data_binary_asian, family = binomial)
model_caide2_apoe_asian <- glm(diag ~ z_caide2 + apoe, 
                               data=data_binary_asian, family = binomial)
model_caide3_apoe_asian <- glm(diag ~ z_caide3 + apoe, 
                               data=data_binary_asian, family = binomial)

add_glance_source_note(tbl_regression(model_caide1_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_asian.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_asian.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_asian.png")

model_caide1_apoe_his <- glm(diag ~ z_caide1 + apoe, 
                                data=data_binary_his, family = binomial)
model_caide2_apoe_his <- glm(diag ~ z_caide2 + apoe, 
                                data=data_binary_his, family = binomial)
model_caide3_apoe_his <- glm(diag ~ z_caide3 + apoe, 
                                data=data_binary_his, family = binomial)

add_glance_source_note(tbl_regression(model_caide1_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caide1_apoe_his.png")

add_glance_source_note(tbl_regression(model_caide2_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide2_apoe_his.png")

add_glance_source_note(tbl_regression(model_caide3_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ CAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caide3_apoe_his.png")

# 2. logistic model: diagnosis ~ mcaide + apoe 
model_mcaide1_apoe_nhw <- glm(diag ~ z_mcaide1 + apoe, 
                            data=data_binary_nhw, family = binomial)
model_mcaide2_apoe_nhw <- glm(diag ~ z_mcaide2 + apoe, 
                            data=data_binary_nhw, family = binomial)
model_mcaide3_apoe_nhw <- glm(diag ~ z_mcaide3 + apoe, 
                            data=data_binary_nhw, family = binomial)

add_glance_source_note(tbl_regression(model_mcaide1_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_nhw.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_nhw.png")

model_mcaide1_apoe_black <- glm(diag ~ z_mcaide1 + apoe, data=data_binary_black, 
                               family = binomial)
model_mcaide2_apoe_black <- glm(diag ~ z_mcaide2 + apoe, data=data_binary_black, 
                                family = binomial)
model_mcaide3_apoe_black <- glm(diag ~ z_mcaide3 + apoe, data=data_binary_black, 
                                family = binomial)

add_glance_source_note(tbl_regression(model_mcaide1_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_black.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_black.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_black.png")

model_mcaide1_apoe_asian <- glm(diag ~ z_mcaide1 + apoe, data=data_binary_asian, 
                               family = binomial)
model_mcaide2_apoe_asian <- glm(diag ~ z_mcaide2 + apoe, data=data_binary_asian, 
                                family = binomial)
model_mcaide3_apoe_asian <- glm(diag ~ z_mcaide3 + apoe, data=data_binary_asian, 
                                family = binomial)

add_glance_source_note(tbl_regression(model_mcaide1_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_asian.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_asian.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_asian.png")

model_mcaide1_apoe_his <- glm(diag ~ z_mcaide1 + apoe, data=data_binary_his, 
                              family = binomial)
model_mcaide2_apoe_his <- glm(diag ~ z_mcaide2 + apoe, data=data_binary_his, 
                              family = binomial)
model_mcaide3_apoe_his <- glm(diag ~ z_mcaide3 + apoe, data=data_binary_his, 
                              family = binomial)

add_glance_source_note(tbl_regression(model_mcaide1_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide1_apoe_his.png")

add_glance_source_note(tbl_regression(model_mcaide2_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide2_apoe_his.png")

add_glance_source_note(tbl_regression(model_mcaide3_apoe_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diagnosis ~ mCAIDE + APOE 
                 (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaide3_apoe_his.png")


# 3 logistic model: diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
#                            apoe + cohort + chol
model_all3_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol3, 
                        data=data_binary_nhw, family = binomial)

model_all2_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol2, 
                      data=data_binary_nhw, family = binomial)

model_all1_nhw <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort, 
                      data=data_binary_nhw, family = binomial)

add_glance_source_note(tbl_regression(model_all1_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort 
                        (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_nhw.png")

add_glance_source_note(tbl_regression(model_all2_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol2 
                        (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_nhw.png")

add_glance_source_note(tbl_regression(model_all3_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ i_bmi + i_systbp + i_age + i_education + gender + 
                        apoe + cohort + hypchol3 
                        (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_nhw.png")


model_all3_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol3, 
                          data=data_binary_black, family = binomial)

model_all2_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol2, 
                        data=data_binary_black, family = binomial)

model_all1_black <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                        data=data_binary_black, family = binomial)

add_glance_source_note(tbl_regression(model_all1_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_black.png")

add_glance_source_note(tbl_regression(model_all2_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_black.png")

add_glance_source_note(tbl_regression(model_all3_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_black.png")


model_all3_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol3, 
                        data=data_binary_asian, family = binomial)

model_all2_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol2, 
                        data=data_binary_asian, family = binomial)

model_all1_asian <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                        data=data_binary_asian, family = binomial)

add_glance_source_note(tbl_regression(model_all1_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_asian.png")

add_glance_source_note(tbl_regression(model_all2_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_asian.png")

add_glance_source_note(tbl_regression(model_all3_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_asian.png")


model_all3_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol3, 
                        data=data_binary_his, family = binomial)

model_all2_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort + hypchol2, 
                        data=data_binary_his, family = binomial)

model_all1_his <- glm(diag ~ i_bmi + i_systbp + i_age + i_education + 
                          gender + apoe + cohort, 
                        data=data_binary_his, family = binomial)

add_glance_source_note(tbl_regression(model_all1_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort (without Cholesterol)") %>%
  gt::gtsave(filename = "model_all1_his.png")

add_glance_source_note(tbl_regression(model_all2_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol2 
                   (with non-imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all2_his.png")

add_glance_source_note(tbl_regression(model_all3_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ i_bmi + i_systbp + i_age + i_education + 
                   gender + apoe + cohort + hypchol3 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_all3_his.png")

########## START HERE ############ Rewrite starting w line 560

#####################################################
# Model: diag ~ caide x apoe
######################################################
# Set reference for caide_apoe:

data_binary_nhw$caide_apoe1 <- factor(data_binary_nhw$caide_apoe1, 
                                             ordered = FALSE )
data_binary_nhw$caide_apoe1 = relevel(data_binary_nhw$caide_apoe1, 
                                             ref = "Mid, e3/e3")
data_binary_nhw$caide_apoe2 <- factor(data_binary_nhw$caide_apoe2, 
                                      ordered = FALSE )
data_binary_nhw$caide_apoe2 = relevel(data_binary_nhw$caide_apoe2, 
                                      ref = "Mid, e3/e3")
data_binary_nhw$caide_apoe3 <- factor(data_binary_nhw$caide_apoe3, 
                                      ordered = FALSE )
data_binary_nhw$caide_apoe3 = relevel(data_binary_nhw$caide_apoe3, 
                                      ref = "Mid, e3/e3")

data_binary_asian$caide_apoe1 <- factor(data_binary_asian$caide_apoe1, 
                                               ordered = FALSE )
data_binary_asian$caide_apoe1 = relevel(data_binary_asian$caide_apoe1, 
                                               ref = "Mid, e3/e3")
data_binary_asian$caide_apoe2 <- factor(data_binary_asian$caide_apoe2, 
                                        ordered = FALSE )
data_binary_asian$caide_apoe2 = relevel(data_binary_asian$caide_apoe2,                                         ref = "Mid, e3/e3")

data_binary_asian$caide_apoe3 <- factor(data_binary_asian$caide_apoe3, 
                                        ordered = FALSE )
data_binary_asian$caide_apoe3 = relevel(data_binary_asian$caide_apoe3, 
                                        ref = "Mid, e3/e3")

data_binary_black$caide_apoe1 <- factor(data_binary_black$caide_apoe1, 
                                               ordered = FALSE )
data_binary_black$caide_apoe1 = relevel(data_binary_black$caide_apoe1, 
                                               ref = "Mid, e3/e3")
data_binary_black$caide_apoe2 <- factor(data_binary_black$caide_apoe2, 
                                        ordered = FALSE )
data_binary_black$caide_apoe2 = relevel(data_binary_black$caide_apoe2, 
                                        ref = "Mid, e3/e3")
data_binary_black$caide_apoe3 <- factor(data_binary_black$caide_apoe3, 
                                        ordered = FALSE )
data_binary_black$caide_apoe3 = relevel(data_binary_black$caide_apoe3, 
                                        ref = "Mid, e3/e3")

data_binary_his$caide_apoe1 <- factor(data_binary_his$caide_apoe1, 
                                             ordered = FALSE )
data_binary_his$caide_apoe1 = relevel(data_binary_his$caide_apoe1, 
                                             ref = "Mid, e3/e3")
data_binary_his$caide_apoe2 <- factor(data_binary_his$caide_apoe2, 
                                      ordered = FALSE )
data_binary_his$caide_apoe2 = relevel(data_binary_his$caide_apoe2, 
                                      ref = "Mid, e3/e3")
data_binary_his$caide_apoe3 <- factor(data_binary_his$caide_apoe3, 
                                      ordered = FALSE )
data_binary_his$caide_apoe3 = relevel(data_binary_his$caide_apoe3, 
                                      ref = "Mid, e3/e3")

# 5. logistic model: diagnosis ~ caide x apoe 
model_caideXapoe1_nhw <- glm(diag ~ caide_apoe1, 
                              data=data_binary_nhw, family = binomial)
model_caideXapoe2_nhw <- glm(diag ~ caide_apoe2, 
                            data=data_binary_nhw, family = binomial)
model_caideXapoe3_nhw <- glm(diag ~ caide_apoe3, 
                            data=data_binary_nhw, family = binomial)


model_caideXapoe1_black <- glm(diag ~ caide_apoe1, 
                                data=data_binary_black, family = binomial)
model_caideXapoe2_black <- glm(diag ~ caide_apoe2, 
                              data=data_binary_black, family = binomial)
model_caideXapoe3_black <- glm(diag ~ caide_apoe3, 
                              data=data_binary_black, family = binomial)


model_caideXapoe1_asian <- glm(diag ~ caide_apoe1, 
                                data=data_binary_asian, family = binomial)
model_caideXapoe2_asian <- glm(diag ~ caide_apoe2, 
                              data=data_binary_asian, family = binomial)
model_caideXapoe3_asian <- glm(diag ~ caide_apoe3, 
                              data=data_binary_asian, family = binomial)

model_caideXapoe1_his <- glm(diag ~ caide_apoe1, 
                              data=data_binary_his, family = binomial)
model_caideXapoe2_his <- glm(diag ~ caide_apoe2, 
                            data=data_binary_his, family = binomial)
model_caideXapoe3_his <- glm(diag ~ caide_apoe3, 
                            data=data_binary_his, family = binomial)


# 6. logistic model: diagnosis ~ mcaide x apoe 
model_mcaideXapoe1_nhw <- glm(diag ~ mcaide_apoe1, 
                             data=data_binary_nhw, family = binomial)
model_mcaideXapoe2_nhw <- glm(diag ~ mcaide_apoe2, 
                             data=data_binary_nhw, family = binomial)
model_mcaideXapoe3_nhw <- glm(diag ~ mcaide_apoe3, 
                             data=data_binary_nhw, family = binomial)


model_mcaideXapoe1_black <- glm(diag ~ mcaide_apoe1, 
                               data=data_binary_black, family = binomial)
model_mcaideXapoe2_black <- glm(diag ~ mcaide_apoe2, 
                               data=data_binary_black, family = binomial)
model_mcaideXapoe3_black <- glm(diag ~ mcaide_apoe3, 
                               data=data_binary_black, family = binomial)

model_mcaideXapoe1_asian <- glm(diag ~ mcaide_apoe1, 
                               data=data_binary_asian, family = binomial)
model_mcaideXapoe2_asian <- glm(diag ~ mcaide_apoe2, 
                               data=data_binary_asian, family = binomial)
model_mcaideXapoe3_asian <- glm(diag ~ mcaide_apoe3, 
                               data=data_binary_asian, family = binomial)

model_mcaideXapoe1_his <- glm(diag ~ mcaide_apoe1, 
                             data=data_binary_his, family = binomial)
model_mcaideXapoe2_his <- glm(diag ~ mcaide_apoe2, 
                             data=data_binary_his, family = binomial)
model_mcaideXapoe3_his <- glm(diag ~ mcaide_apoe3, 
                             data=data_binary_his, family = binomial)

add_glance_source_note(tbl_regression(model_caideXapoe1_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ caideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe1_nhw.png")

add_glance_source_note(tbl_regression(model_caideXapoe2_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ caideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe2_nhw.png")

add_glance_source_note(tbl_regression(model_caideXapoe3_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ caideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe3_nhw.png")

add_glance_source_note(tbl_regression(model_caideXapoe1_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ caideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe1_black.png")

add_glance_source_note(tbl_regression(model_caideXapoe2_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ caideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe2_black.png")

add_glance_source_note(tbl_regression(model_caideXapoe3_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ caideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe3_black.png")

add_glance_source_note(tbl_regression(model_caideXapoe1_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ caideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe1_asian.png")

add_glance_source_note(tbl_regression(model_caideXapoe2_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ caideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe2_asian.png")

add_glance_source_note(tbl_regression(model_caideXapoe3_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ caideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe3_asian.png")

add_glance_source_note(tbl_regression(model_caideXapoe1_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ caideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe1_his.png")

add_glance_source_note(tbl_regression(model_caideXapoe2_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ caideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe2_his.png")

add_glance_source_note(tbl_regression(model_caideXapoe3_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ caideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_caideXapoe3_his.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe1_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ mcaideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe1_nhw.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe2_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ mcaideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe2_nhw.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe3_nhw)) %>% as_gt %>%
  gt::tab_header(title = 
                   "NHW: Diag ~ mcaideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe3_nhw.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe1_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ mcaideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe1_black.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe2_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ mcaideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe2_black.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe3_black)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Black: Diag ~ mcaideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe3_black.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe1_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ mcaideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe1_asian.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe2_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ mcaideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe2_asian.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe3_asian)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Asian: Diag ~ mcaideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe3_asian.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe1_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ mcaideXapoe 
                   (without Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe1_his.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe2_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ mcaideXapoe 
                   (with NAs Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe2_his.png")

add_glance_source_note(tbl_regression(model_mcaideXapoe3_his)) %>% as_gt %>%
  gt::tab_header(title = 
                   "Hispanic: Diag ~ mcaideXapoe 
                   (with imputed Cholesterol)") %>%
  gt::gtsave(filename = "model_mcaideXapoe3_his.png")






