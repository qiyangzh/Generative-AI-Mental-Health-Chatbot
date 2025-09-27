########################################################################################################
#AI Mental Health Interventions, JMIR
########################################################################################################
# Authors: Qiyang Zhang
# Contact: qiyang39@nus.edu.sg
# Created: 2024/11/14

########################################################################################################
# Initial Set-up
########################################################################################################
# Clear workspace
rm(list=ls(all=TRUE))

# Load packages
test<-require(googledrive)   #all gs_XXX() functions for reading data from Google
if (test == FALSE) {
  install.packages("googledrive")
  require(googledrive)
}
test<-require(googlesheets4)   #all gs_XXX() functions for reading data from Google
if (test == FALSE) {
  install.packages("googlesheets4")
  require(googlesheets4)
}
test<-require(plyr)   #rename()
if (test == FALSE) {
  install.packages("plyr")
  require(plyr)
}
test<-require(metafor)   #escalc(); rma();
if (test == FALSE) {
  install.packages("metafor")
  require(metafor)
}
test<-require(robumeta)
if (test == FALSE) {
  install.packages("robumeta")
  require(robumeta)
}
test<-require(weightr) #selection modeling
if (test == FALSE) {
  install.packages("weightr")
  require(weightr)
}
test<-require(clubSandwich) #coeftest
if (test == FALSE) {
  install.packages("clubSandwich")
  require(clubSandwich)
}
test<-require(tableone)   #CreateTableOne()
if (test == FALSE) {
  install.packages("tableone")
  require(tableone)
}
test<-require(flextable)   
if (test == FALSE) {
  install.packages("flextable")
  require(flextable)
}
test<-require(officer)   
if (test == FALSE) {
  install.packages("officer")
  require(officer)
}
test<-require(tidyverse)   
if (test == FALSE) {
  install.packages("tidyverse")
  require(tidyverse)
}
test<-require(ggrepel)   
if (test == FALSE) {
  install.packages("ggrepel")
  require(ggrepel)
}
test<-require(readxl)   
if (test == FALSE) {
  install.packages("readxl")
  require(readxl)
}
rm(test)

########################################################################################################
# Load data
########################################################################################################
# set up to load from Google
drive_auth(email = "zhangqiyang0329@gmail.com")
id <- drive_find(pattern = "AI_Chatbot_Well-being", type = "spreadsheet")$id[1]

# load findings and studies
gs4_auth(email = "zhangqiyang0329@gmail.com")
findings <- read_sheet(id, sheet = "Findings", col_types = "c")
studies <- read_sheet(id, sheet = "Studies", col_types = "c")   # includes separate effect sizes for each finding from a study

# setwd("~/Desktop/Research2025/AI_WellBeing")
# # set up to load from Google
# findings <- read_excel("AI_Chatbot_Well-being.xlsx", sheet = "Findings")
# studies <- read_excel("AI_Chatbot_Well-being.xlsx", sheet = "Studies")

rm(id)

########################################################################################################
# Clean data
########################################################################################################
# remove any empty rows & columns
studies <- subset(studies, is.na(studies$Study)==FALSE)
findings <- subset(findings, is.na(findings$Study)==FALSE)

studies <- subset(studies, studies$Drop==2 | is.na(studies$Drop)==TRUE)
findings <- subset(findings, is.na(findings$Drop) == TRUE)

# merge dataframes
full <- merge(studies, findings, by = c("Study"), all = TRUE, suffixes = c(".s", ".f"))
full <- subset(full, is.na(full$Authors.s)!=TRUE)

# format to correct variable types
nums <- c("Treatment.N.original", "Control.N.original",
          "T_Mean_Pre", "T_SD_Pre", "C_Mean_Pre", 
          "C_SD_Pre", "T_Mean_Post", "T_SD_Post", 
          "C_Mean_Post", "C_SD_Post", "Clustered",
          "Age.mean","Duration.weeks", "Follow-up",
          "Targeted", "Students", "Clinical", "Personalized",
          "Sample size", "Female", "HumanAssistance")

full[nums] <- lapply(full[nums], as.numeric)
rm(nums)

###############################################################
#Create unique identifiers (ES, study, program)
###############################################################
full$ESId <- as.numeric(rownames(full))
full$StudyID <- as.numeric(as.factor(full$Study))
summary(full$StudyID)
########################################################################################################
# Prep data
########################################################################################################
##### Calculate ESs #####
# calculate pretest ES, SMD is standardized mean difference
full <- escalc(measure = "SMD", m1i = T_Mean_Pre, sd1i = T_SD_Pre, n1i = Treatment.N.original,
               m2i = C_Mean_Pre, sd2i = C_SD_Pre, n2i = Control.N.original, data = full)
full$vi <- NULL
full <- plyr::rename(full, c("yi" = "ES_Pre"))

# calculate posttest ES
full <- escalc(measure = "SMD", m1i = T_Mean_Post, sd1i = T_SD_Post, n1i = Treatment.N.original,
               m2i = C_Mean_Post, sd2i = C_SD_Post, n2i = Control.N.original, data = full)
full$vi <- NULL
full <- plyr::rename(full, c("yi" = "ES_Post"))

# calculate DID (post ES - pre ES)
full$ES_DID <- full$ES_Post - full$ES_Pre

# put various ES together.  Options:
## 1) used reported ES (so it should be in the Effect.Size column, nothing to do)
## 2) Effect.Size is NA, and DID not missing, replace with that
full$ES_DID[which(is.na(full$ES_DID)==TRUE & is.na(full$ES)==FALSE)] <- full$ES[which(is.na(full$ES_DID)==TRUE & is.na(full$ES)==FALSE)]
full$Effect.Size <- as.numeric(full$ES_DID)
full$Effect.Size <- full$Effect.Size*-1

###############################################################
#Calculate meta-analytic variables: Sample sizes
###############################################################
#create full sample/total clusters variables
full$Sample <- full$Treatment.N.original + full$Control.N.original

################################################################
# Calculate meta-analytic variables: Variances (Lipsey & Wilson, 2000, Eq. 3.23)
################################################################
#calculate standard errors
full$se<-sqrt(((full$Treatment.N.original+full$Control.N.original)/(full$Treatment.N.original*full$Control.N.original))+((full$Effect.Size*full$Effect.Size)/(2*(full$Treatment.N.original+full$Control.N.original))))

#calculate variance
full$var<-full$se*full$se
#####################################
#Compare generative and rule-based
#####################################
# V_list <- impute_covariance_matrix(vi=full$var, cluster=full$StudyID, r=0.8)
# 
# MVnull <- rma.mv(yi=Effect.Size,
#                  V=V_list,
#                  random=~1 | StudyID/ESId,
#                  test="t",
#                  data=full,
#                  method="REML")
# MVnull
# #View(full[c("Study","Effect.Size", "T_Mean_Pre", "Outcomes", "Targeted", "Clinical", "ActiveorPassive", "Continent")])
# 
# #t-test of each covariate#
# MVnull.coef <- coef_test(MVnull, cluster=full$StudyID, vcov="CR2")
# MVnull.coef
# 
# terms <- c("Response.generation.approach")
# formula <- reformulate(termlabels = c(terms))
# 
# MVfull <- rma.mv(yi=Effect.Size,
#                  V=V_list,
#                  mods=formula,
#                  random=~1 | StudyID/ESId,
#                  test="t",
#                  data=full,
#                  method="REML")
# MVfull
# #t-test of each covariate#
# MVfull.coef <- coef_test(MVfull, cluster=full$StudyID, vcov="CR2")
# MVfull.coef

#########################
########################################
#meta-regression
########################################
#subset to generative only
full <- subset(full, is.na(full$Drop.s)==TRUE)

#Centering, when there is missing value, this won't work
full$FemalePercent <- full$Female/full$Sample.size
full$FiftyPercentFemale <- 0
full$FiftyPercentFemale[which(full$`FemalePercent` > 0.50)] <- 1
full$Students.c <- full$Students - mean(full$Students)
full$Clustered.c <- full$Clustered - mean(full$Clustered)
full$Follow.up.c <- full$Follow.up - mean(full$Follow.up)
full$Clinical.c <- full$Clinical - mean(full$Clinical)
full$FiftyPercentFemale.c <- full$FiftyPercentFemale - mean(full$FiftyPercentFemale)
full$Personalized.c <- full$Personalized - mean(full$Personalized)
#full$Self.guided.c <- full$Self.guided - mean(full$Self.guided)
full$Targeted.c <- full$Targeted - mean(full$Targeted)

full$ESId <- as.numeric(rownames(full))
full$StudyID <- as.numeric(as.factor(full$Study))
summary(full$StudyID)
#sensitivity analysis, vary correlation r
#Null Model
V_list <- impute_covariance_matrix(vi=full$var, cluster=full$StudyID, r=0.8)
#test="t"
MVnull <- metafor::rma.mv(yi=Effect.Size,
                 V=V_list,
                 random=~1 | StudyID/ESId,
                 dfs  = "contain",
                 data=full,
                 method="REML",
                 test = "t")
MVnull

#t-test of each covariate# # small-sample df approximation
MVnull.coef <- coef_test(MVnull, cluster=full$StudyID, vcov="CR2", test = "Satterthwaite")
MVnull.coef

### Output prediction interval ###
dat_clean <- data.frame(yi = full$Effect.Size, se_g = full$se)
dat_clean <- na.omit(dat_clean)

yi_clean <- dat_clean$yi
se_g_clean <- dat_clean$se_g
install.packages("pimeta")
library(pimeta)

pima_result <- pima(yi_clean, se_g_clean, method = "HK")  # Using the Hartung-Knapp method

print(pima_result)

# yi_uni, vi_uni should be one effect and its sampling variance per study
by_study_outcome <- full %>%
  dplyr::group_by(Study, Outcomes) %>%
  dplyr::summarise(
    Effect.Size = mean(Effect.Size, na.rm = TRUE),
    var = mean(var, na.rm = TRUE),
    n_rows = n()  # how many rows were averaged
  )
fit_HKSJ_SJ <- metafor::rma(
  yi     = Effect.Size,
  vi     = var,
  data = by_study_outcome,
  method = "SJ",        # Sidik–Jonkman tau^2
  test   = "knha"       # Hartung–Knapp with HA adjustment
)
summary(fit_HKSJ_SJ)
pred_HKSJ <- predict(fit_HKSJ_SJ)
pred_HKSJ
#####
#subgroup
######
# Ensure the outcome variable is correctly coded
by_study_outcome$Outcomes <- factor(by_study_outcome$Outcomes,
                                   levels = c("Depression", "Anxiety", "Stress", "Negative Affect/Mood"))

# Subgroup: depression
fit_dep <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = subset(by_study_outcome, Outcomes == "Depression"),
  method = "SJ",
  test   = "knha"
)
summary(fit_dep)

# Subgroup: anxiety
fit_anx <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = subset(by_study_outcome, Outcomes == "Anxiety"),
  method = "SJ",
  test   = "knha"
)
summary(fit_anx)

# Subgroup: NA
fit_negative <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = subset(by_study_outcome, Outcomes == "Negative Affect/Mood"),
  method = "SJ",
  test   = "knha"
)
summary(fit_negative)

# Subgroup: Stress
fit_stress <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = subset(by_study_outcome, Outcomes == "Stress"),
  method = "SJ",
  test   = "knha"
)
summary(fit_stress)
# Depression
pred_dep <- predict(fit_dep)
pred_dep$pi.lb; pred_dep$pi.ub

# Anxiety
pred_anx <- predict(fit_anx)
pred_anx$pi.lb; pred_anx$pi.ub

# Negative affect/mood
pred_neg <- predict(fit_negative)
pred_neg$pi.lb; pred_neg$pi.ub

# Stress
pred_stress <- predict(fit_stress)
pred_stress$pi.lb; pred_stress$pi.ub
########################
#Moderator analysis, single moderator per model, Test whether subgroups differ using a moderator in a single model
#####################
#Outcome as a moderator
by_study_outcome <- by_study_outcome %>%
  mutate(
    Outcomes = factor(Outcomes,
                      levels = c("Depression", "Anxiety", "Stress", "Negative Affect/Mood"))
  )

# Set Depression as the reference level
by_study_outcome$Outcomes <- relevel(by_study_outcome$Outcomes, ref = "Depression")

# Fit one-moderator model (univariate; HKSJ + SJ)
fit_outcomes <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = by_study_outcome,
  mods   = ~ Outcomes,         # coefficients are differences vs Depression
  method = "SJ",               # Sidik–Jonkman tau^2
  test   = "knha"              # Hartung–Knapp–Sidik–Jonkman inference
)

summary(fit_outcomes)

# Robust (CR2) SEs clustered on Study to address residual dependence
# Align cluster IDs to the rows actually used in the fit
mf <- model.frame(fit_outcomes)
# After fitting fit_outcomes
# Make sure Study exists in by_study_outcome
stopifnot("Study" %in% names(by_study_outcome))

# Subset the Study vector to the rows used in the model
cluster_vec <- by_study_outcome$Study[fit_outcomes$not.na]

# Now run CR2 robust inference
rob_outcomes <- coef_test(fit_outcomes,
                          cluster = cluster_vec,
                          vcov = "CR2",
                          test = "Satterthwaite")

rob_outcomes

#########
#"ActiveorPassive"
#########
by_study_outcome <- by_study_outcome %>%
  ungroup() %>%
  left_join(
    full %>% distinct(Study, ActiveorPassive),
    by = "Study"
  )

# Recode ActiveorPassive: 0/1 -> "Passive"/"Active" (adjust if your values are already strings)
by_study_outcome <- by_study_outcome %>%
  mutate(
    ActiveorPassive = factor(ActiveorPassive,
                             levels = c("passive", "active"),         # use c("Passive","Active") if already strings
                             labels = c("passive", "active"))
  )

# Set "Passive" as the reference level
by_study_outcome$ActiveorPassive <- relevel(by_study_outcome$ActiveorPassive, ref = "passive")

# Fit one-moderator model (univariate; HKSJ + SJ)
fit_AP <- rma(
  yi     = Effect.Size,
  vi     = var,
  data   = by_study_outcome,
  mods   = ~ ActiveorPassive,   # coefficient is Active vs Passive
  method = "SJ",                # Sidik–Jonkman tau^2
  test   = "knha"               # Hartung–Knapp–Sidik–Jonkman inference
)

summary(fit_AP)

# Robust (CR2) SEs clustered on Study to address residual dependence
# Subset Study to the rows used in the model to align cluster IDs
stopifnot("Study" %in% names(by_study_outcome))
cluster_vec <- by_study_outcome$Study[fit_AP$not.na]

rob_AP <- coef_test(fit_AP,
                    cluster = cluster_vec,
                    vcov    = "CR2",
                    test    = "Satterthwaite")

rob_AP

#######
#Age
########
# Build by_study_outcome if not already done
by_study_outcome <- full %>%
  group_by(Study, Outcomes) %>%
  summarise(
    Effect.Size = mean(Effect.Size, na.rm = TRUE),
    var         = mean(var, na.rm = TRUE),
    n_rows      = n(),
    .groups     = "drop"
  )

# Join Age from full (ensure one Age per Study)
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, Age), by = "Study")

# Normalize Age labels to exactly: "early adulthood", "middle adulthood", "late adulthood"
by_study_outcome <- by_study_outcome %>%
  mutate(
    Age = tolower(Age),
    Age = case_when(
      Age %in% c("late adulthood","late adult","late")    ~ "late adulthood",
      Age %in% c("middle adulthood","middle adult","mid") ~ "middle adulthood",
      Age %in% c("early adulthood","early adult","early") ~ "early adulthood",
      TRUE ~ Age
    ),
    Age3 = factor(Age, levels = c("early adulthood","middle adulthood","late adulthood"))
  )

# Drop rows with missing moderator or essentials
dat_age <- subset(by_study_outcome, !is.na(Age3) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_age) > 0)

# Set reference = early adulthood and fit HKSJ model (SJ tau^2)
dat_age$Age3 <- relevel(dat_age$Age3, ref = "early adulthood")

fit_age <- rma(
  yi     = dat_age$Effect.Size,
  vi     = dat_age$var,
  mods   = ~ Age3,       # coefficients: middle−early, late−early
  method = "SJ",
  test   = "knha",
  data   = dat_age
)
summary(fit_age)

# CR2 robust SEs clustered by Study, aligned to rows used in the fit
cluster_vec <- dat_age$Study[fit_age$not.na]
rob_age <- coef_test(fit_age, cluster = cluster_vec, vcov = "CR2", test = "Satterthwaite")
rob_age

############
#Clinical
############
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, Clinical), by = "Study")

# 3) Recode Clinical to a factor with explicit levels and reference
# If Clinical is numeric (0/1) or character "0"/"1", coerce then label
by_study_outcome <- by_study_outcome %>%
  mutate(
    Clinical = suppressWarnings(as.numeric(Clinical)),
    Clinical = factor(Clinical, levels = c(0, 1), labels = c("Nonclinical", "Clinical"))
  )
# Set Nonclinical as reference (already first level); relevel if needed:
by_study_outcome$Clinical <- relevel(by_study_outcome$Clinical, ref = "Nonclinical")

# 4) Drop rows with missing moderator or essentials
dat_cli <- subset(by_study_outcome, !is.na(Clinical) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_cli) > 0)

# 5) Fit one-moderator model (HKSJ + SJ)
fit_clin <- rma(
  yi     = dat_cli$Effect.Size,
  vi     = dat_cli$var,
  mods   = ~ Clinical,   # coefficient is Clinical vs Nonclinical
  method = "SJ",         # Sidik–Jonkman tau^2
  test   = "knha",       # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_cli
)
summary(fit_clin)

# 6) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_cli$Study[fit_clin$not.na]

rob_clin <- coef_test(fit_clin,
                      cluster = cluster_vec,
                      vcov    = "CR2",
                      test    = "Satterthwaite")
rob_clin

#############
#Personalized
##############
# 2) Join Personalized (0/1) from full (ensure one value per Study)
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, Personalized), by = "Study")

# 3) Recode Personalized to a factor with explicit levels and reference
# Handles numeric 0/1 or character "0"/"1" by coercing to numeric first
by_study_outcome <- by_study_outcome %>%
  mutate(
    Personalized = suppressWarnings(as.numeric(Personalized)),
    Personalized = factor(Personalized, levels = c(0, 1), labels = c("No", "Yes"))
  )
# Set "No" as reference (already first level); relevel explicitly if you like:
by_study_outcome$Personalized <- relevel(by_study_outcome$Personalized, ref = "No")

# 4) Drop rows with missing moderator or essentials
dat_per <- subset(by_study_outcome, !is.na(Personalized) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_per) > 0)

# 5) Fit one-moderator HKSJ model (SJ tau^2)
fit_per <- rma(
  yi     = dat_per$Effect.Size,
  vi     = dat_per$var,
  mods   = ~ Personalized,   # coefficient is Yes vs No
  method = "SJ",             # Sidik–Jonkman tau^2
  test   = "knha",           # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_per
)
summary(fit_per)

# 6) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_per$Study[fit_per$not.na]

rob_per <- coef_test(fit_per,
                     cluster = cluster_vec,
                     vcov    = "CR2",
                     test    = "Satterthwaite")
rob_per

################
#WEIRD
#################
# 2) Join WEIRD (0/1) from full (ensure one value per Study)
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, WEIRD), by = "Study")

# 3) Recode WEIRD to a factor with explicit levels and reference
# Handles numeric 0/1 or character "0"/"1" by coercing to numeric first
by_study_outcome <- by_study_outcome %>%
  mutate(
    WEIRD = suppressWarnings(as.numeric(WEIRD)),
    WEIRD = factor(WEIRD, levels = c(0, 1), labels = c("Non-WEIRD", "WEIRD"))
  )
# Set "Non-WEIRD" as reference (already first level); relevel explicitly if desired:
by_study_outcome$WEIRD <- relevel(by_study_outcome$WEIRD, ref = "Non-WEIRD")

# 4) Drop rows with missing moderator or essentials
dat_weird <- subset(by_study_outcome, !is.na(WEIRD) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_weird) > 0)

# 5) Fit one-moderator HKSJ model (SJ tau^2)
fit_weird <- rma(
  yi     = dat_weird$Effect.Size,
  vi     = dat_weird$var,
  mods   = ~ WEIRD,      # coefficient is WEIRD vs Non-WEIRD
  method = "SJ",         # Sidik–Jonkman tau^2
  test   = "knha",       # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_weird
)
summary(fit_weird)

# 6) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_weird$Study[fit_weird$not.na]

rob_weird <- coef_test(fit_weird,
                       cluster = cluster_vec,
                       vcov    = "CR2",
                       test    = "Satterthwaite")
rob_weird

#################
#Human Assisted
################
# 1) Join HumanAssistance (0/1) from full into by_study_outcome
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, HumanAssistance), by = "Study")

# 2) Recode HumanAssistance to a factor with explicit levels and reference
by_study_outcome <- by_study_outcome %>%
  mutate(
    HumanAssistance = suppressWarnings(as.numeric(HumanAssistance)),
    HumanAssistance = factor(HumanAssistance, levels = c(0, 1), labels = c("No", "Yes"))
  )
# Set "No" as the reference level
by_study_outcome$HumanAssistance <- relevel(by_study_outcome$HumanAssistance, ref = "No")

# 3) Drop rows with missing moderator or essentials
dat_ha <- subset(by_study_outcome, !is.na(HumanAssistance) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_ha) > 0)

# 4) Fit one-moderator model (HKSJ + SJ)
fit_HA <- rma(
  yi     = dat_ha$Effect.Size,
  vi     = dat_ha$var,
  mods   = ~ HumanAssistance,   # coefficient is Yes vs No
  method = "SJ",                # Sidik–Jonkman tau^2
  test   = "knha",              # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_ha
)
summary(fit_HA)

# 5) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_ha$Study[fit_HA$not.na]
rob_HA <- coef_test(fit_HA,
                    cluster = cluster_vec,
                    vcov    = "CR2",
                    test    = "Satterthwaite")
rob_HA

###########
#FiftyPercent
############
# 1) Join FiftyPercentFemale (0/1) from full into by_study_outcome
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, FiftyPercentFemale), by = "Study")

# 2) Recode FiftyPercentFemale to a factor with explicit levels and reference
by_study_outcome <- by_study_outcome %>%
  mutate(
    FiftyPercentFemale = suppressWarnings(as.numeric(FiftyPercentFemale)),
    FiftyPercentFemale = factor(FiftyPercentFemale, levels = c(0, 1), labels = c("No", "Yes"))
  )
# Set "No" as the reference level
by_study_outcome$FiftyPercentFemale <- relevel(by_study_outcome$FiftyPercentFemale, ref = "No")

# 3) Drop rows with missing moderator or essentials
dat_ff <- subset(by_study_outcome, !is.na(FiftyPercentFemale) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_ff) > 0)

# 4) Fit one-moderator model (HKSJ + SJ)
fit_FF <- rma(
  yi     = dat_ff$Effect.Size,
  vi     = dat_ff$var,
  mods   = ~ FiftyPercentFemale,  # coefficient is Yes vs No
  method = "SJ",                  # Sidik–Jonkman tau^2
  test   = "knha",                # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_ff
)
summary(fit_FF)

# 5) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_ff$Study[fit_FF$not.na]
rob_FF <- coef_test(fit_FF,
                    cluster = cluster_vec,
                    vcov    = "CR2",
                    test    = "Satterthwaite")
rob_FF
#############
#Social.function
#############
# 1) Join Social.function (two categories) from full into by_study_outcome
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, Social.function), by = "Study")

# 2) Recode Social.function to a factor with explicit levels and reference
# Normalize common variants, then set levels exactly as requested
by_study_outcome <- by_study_outcome %>%
  mutate(
    Social.function = tolower(trimws(Social.function)),
    Social.function = case_when(
      Social.function %in% c("task-oriented","task oriented","task")   ~ "task-oriented",
      Social.function %in% c("social-oriented","social oriented","social") ~ "social-oriented",
      TRUE ~ NA_character_
    ),
    Social.function = factor(Social.function, levels = c("task-oriented", "social-oriented"))
  )

# Set "task-oriented" as the reference level
by_study_outcome$Social.function <- relevel(by_study_outcome$Social.function, ref = "task-oriented")

# 3) Drop rows with missing moderator or essentials
dat_sf <- subset(by_study_outcome, !is.na(Social.function) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_sf) > 0)

# 4) Fit one-moderator model (HKSJ + SJ)
fit_SF <- rma(
  yi     = dat_sf$Effect.Size,
  vi     = dat_sf$var,
  mods   = ~ Social.function,   # coefficient is social-oriented vs task-oriented
  method = "SJ",                # Sidik–Jonkman tau^2
  test   = "knha",              # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_sf
)
summary(fit_SF)

# 5) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_sf$Study[fit_SF$not.na]
rob_SF <- coef_test(fit_SF,
                    cluster = cluster_vec,
                    vcov    = "CR2",
                    test    = "Satterthwaite")
rob_SF
#############
#Duration.weeks
#duration
###############
# 1) Join Duration.weeks (numeric) from full into by_study_outcome
by_study_outcome <- by_study_outcome %>%
  left_join(full %>% distinct(Study, Duration.weeks), by = "Study") %>%
  mutate(
    Duration_weeks = suppressWarnings(as.numeric(Duration.weeks)),
    # Center/scale for interpretability and stability (slope per 1 SD)
    Duration_c     = scale(Duration_weeks, center = TRUE, scale = TRUE)[,1]
    # If you prefer slope per 4 weeks instead, use:
    # Duration_4wk   = (Duration_weeks - mean(Duration_weeks, na.rm = TRUE)) / 4
  )

# 2) Drop rows with missing moderator or essentials
dat_dur <- subset(by_study_outcome, !is.na(Duration_c) & !is.na(Effect.Size) & !is.na(var))
stopifnot(nrow(dat_dur) > 0)

# 3) Fit one-moderator model (HKSJ + SJ) with continuous duration
fit_dur <- rma(
  yi     = dat_dur$Effect.Size,
  vi     = dat_dur$var,
  mods   = ~ Duration_c,   # slope per 1 SD increase in Duration.weeks
  method = "SJ",           # Sidik–Jonkman tau^2
  test   = "knha",         # Hartung–Knapp–Sidik–Jonkman inference
  data   = dat_dur
)
summary(fit_dur)

# 4) CR2 robust SEs clustered by Study; align cluster IDs to model rows
cluster_vec <- dat_dur$Study[fit_dur$not.na]
rob_dur <- coef_test(fit_dur,
                     cluster = cluster_vec,
                     vcov    = "CR2",
                     test    = "Satterthwaite")
rob_dur
#Method: "ActiveorPassive"
#Unit:"Clinical", "Age", "FiftyPercentFemale"
#Treatment:"Duration.weeks", "Personalized", "Self.guided", "Modality","Social.function"
#Outcome:"Outcomes"
#Setting: "WEIRD"
terms <- c("ActiveorPassive", "WEIRD",
           "Clinical", "Age", "FiftyPercentFemale",
           "Duration.weeks", "Personalized", "HumanAssistance",
           "Modality","Social.function", "Outcomes")
#interact <- c("Social.function*Age")
#formula <- reformulate(termlabels = c(terms, interact))
formula <- reformulate(termlabels = c(terms))
formula

MVfull <- rma.mv(yi=Effect.Size,
                 V=V_list,
                 mods=formula,
                 random=~1 | StudyID/ESId,
                 test="t",
                 data=full,
                 method="REML")
MVfull
#t-test of each covariate#
MVfull.coef <- coef_test(MVfull, cluster=full$StudyID, vcov="CR2")
MVfull.coef

###########################
#forest plot
###########################
study_averages <- full %>%
  dplyr::group_by(StudyID) %>%
  dplyr::summarise(avg_effect_size = mean(Effect.Size, na.rm = TRUE),
            across(everything(), ~ first(.)))

MVnull <- robu(formula = Effect.Size ~ 1, studynum = StudyID, data = study_averages, var.eff.size = var)
m.gen <- meta::metagen(TE = Effect.Size,
                 seTE = se,
                 studlab = Study,
                 data = study_averages,
                 sm = "SMD",
                 common = FALSE,
                 random = TRUE,
                 method.tau = "REML",
                 method.random.ci = "HK")
summary(m.gen)
png(file = "forestplot.png", width = 2800, height = 3000, res = 300)

meta::forest(m.gen,
             sortvar = TE,
             prediction = TRUE,
             print.tau2 = FALSE,
             leftlabs = c("Author", "g", "SE"))
dev.off()
#################################################################################
# Marginal Means
#################################################################################
# re-run model for each moderator to get marginal means for each #

# set up table to store results
means <- data.frame(moderator = character(0), group = character(0), beta = numeric(0), SE = numeric(0), 
                    tstat = numeric(0), df = numeric(0), p_Satt = numeric(0))

mods <- c("as.factor(Personalized)", "as.factor(Modality)", "as.factor(Social.function)", 
          "as.factor(Clinical)", "as.factor(ActiveorPassive)", "as.factor(WEIRD)",
          "as.factor(Age)", "as.factor(FiftyPercentFemale)", "as.factor(HumanAssistance)",
          "as.factor(Outcomes)")

for(i in 1:length(mods)){
  # i <- 1
  formula <- reformulate(termlabels = c(mods[i], terms, "-1"))   # Worth knowing - if you duplicate terms, it keeps the first one
  mod_means <- rma.mv(yi=Effect.Size, #effect size
                      V = V_list, #variance (tHIS IS WHAt CHANGES FROM HEmodel)
                      mods = formula, #ADD COVS HERE
                      random = ~1 | StudyID/ESId, #nesting structure
                      test= "t", #use t-tests
                      data=full, #define data
                      method="REML") #estimate variances using REML
  coef_mod_means <- as.data.frame(coef_test(mod_means,#estimation model above
                                            cluster=full$StudyID, #define cluster IDs
                                            vcov = "CR2")) #estimation method (CR2 is best)
  # limit to relevant rows (the means you are interested in)
  coef_mod_means$moderator <- gsub(x = mods[i], pattern = "as.factor", replacement = "")
  coef_mod_means$group <- rownames(coef_mod_means)
  rownames(coef_mod_means) <- c()
  coef_mod_means <- subset(coef_mod_means, substr(start = 1, stop = nchar(mods[i]), x = coef_mod_means$group)== mods[i])
  coef_mod_means$group <- substr(x = coef_mod_means$group, start = nchar(mods[i])+1, stop = nchar(coef_mod_means$group))
  means <- dplyr::bind_rows(means, coef_mod_means)
}
means
#################################################################################
# Heterogeneity
#################################################################################
# 95% prediction intervals
print(PI_upper <- MVfull$b[1] + (1.96*sqrt(MVfull$sigma2[1] + MVfull$sigma2[2])))
print(PI_lower <- MVfull$b[1] - (1.96*sqrt(MVfull$sigma2[1] + MVfull$sigma2[2])))

#################################################################################
#Create Descriptives Table
#################################################################################
# identify variables for descriptive tables (study-level and outcome-level)
#"Follow.up.Duration.weeks","Age.mean","FemalePercent",Duration.weeks",
vars_study <- c("WEIRD",
                "Clinical", "Age", "FiftyPercentFemale",
                "Personalized", "HumanAssistance", 
                "Modality","Social.function")
vars_outcome <- c("Outcomes", "Clustered", "ActiveorPassive", "Follow.up")

# To make this work, you will need a df that is at the study-level for study-level 
# variables (such as research design) you may have already created this (see above, with study-level ESs), but if you didn't, here is an easy way:
# 1) make df with *only* the study-level variables of interest and studyIDs in it
study_level_full <- full[c("StudyID", "WEIRD",
                           "Clinical", "Age", "FiftyPercentFemale",
                           "Personalized", "HumanAssistance", 
                           "Modality","Social.function")]
# 2) remove duplicated rows
study_level_full <- unique(study_level_full)
# 3) make sure it is the correct number of rows (should be same number of studies you have)
length(study_level_full$StudyID)==length(unique(study_level_full$StudyID))
# don't skip step 3 - depending on your data structure, some moderators can be
# study-level in one review, but outcome-level in another

# create the table "chunks"
table_study_df <- as.data.frame(print(CreateTableOne(vars = vars_study, data = study_level_full, 
                                                     includeNA = TRUE, 
                                                     factorVars = c("WEIRD",
                                                                    "Clinical", "Age", "FiftyPercentFemale",
                                                                    "Personalized", "HumanAssistance", 
                                                                    "Modality","Social.function")), 
                                      showAllLevels = TRUE))
table_outcome_df <- as.data.frame(print(CreateTableOne(vars = vars_outcome, data = full, includeNA = TRUE,
                                                       factorVars = c("Outcomes", "Clustered", "ActiveorPassive", "Follow.up")), 
                                        showAllLevels = TRUE))
rm(vars_study, vars_outcome)

################################
# Descriptives Table Formatting
################################
table_study_df$Category <- row.names(table_study_df)
rownames(table_study_df) <- c()
table_study_df <- table_study_df[c("Category", "level", "Overall")]
table_study_df$Category[which(substr(table_study_df$Category, 1, 1)=="X")] <- NA
table_study_df$Category <- gsub(pattern = "\\..mean..SD..", replacement = "", x = table_study_df$Category)
table_study_df$Overall <- gsub(pattern = "\\( ", replacement = "\\(", x = table_study_df$Overall)
table_study_df$Overall <- gsub(pattern = "\\) ", replacement = "\\)", x = table_study_df$Overall)
table_study_df$Category <- gsub(pattern = "\\.", replacement = "", x = table_study_df$Category)
table_study_df$Category[which(table_study_df$Category=="n")] <- "Total Studies"
table_study_df$level[which(table_study_df$level=="1")] <- "Yes"
table_study_df$level[which(table_study_df$level=="0")] <- "No                                                                              "
# fill in blank columns (to improve merged cells later)
for(i in 1:length(table_study_df$Category)) {
  if(is.na(table_study_df$Category[i])) {
    table_study_df$Category[i] <- table_study_df$Category[i-1]
  }
}

table_outcome_df$Category <- row.names(table_outcome_df)
rownames(table_outcome_df) <- c()
table_outcome_df <- table_outcome_df[c("Category", "level", "Overall")]
table_outcome_df$Category[which(substr(table_outcome_df$Category, 1, 1)=="X")] <- NA
table_outcome_df$Overall <- gsub(pattern = "\\( ", replacement = "\\(", x = table_outcome_df$Overall)
table_outcome_df$Overall <- gsub(pattern = "\\) ", replacement = "\\)", x = table_outcome_df$Overall)
table_outcome_df$Category <- gsub(pattern = "\\.", replacement = "", x = table_outcome_df$Category)
table_outcome_df$Category[which(table_outcome_df$Category=="n")] <- "Total Effect Sizes"
# fill in blank columns (to improve merged cells later)
for(i in 1:length(table_outcome_df$Category)) {
  if(is.na(table_outcome_df$Category[i])) {
    table_outcome_df$Category[i] <- table_outcome_df$Category[i-1]
  }
}

########################
#Output officer
########################
myreport<-read_docx()
# Descriptive Table
myreport <- body_add_par(x = myreport, value = "Table 4: Descriptive Statistics", style = "Normal")
descriptives_study <- flextable(head(table_study_df, n=nrow(table_study_df)))
descriptives_study <- add_header_lines(descriptives_study, values = c("Study Level"), top = FALSE)
descriptives_study <- theme_vanilla(descriptives_study)
descriptives_study <- merge_v(descriptives_study, j = c("Category"))
myreport <- body_add_flextable(x = myreport, descriptives_study)

descriptives_outcome <- flextable(head(table_outcome_df, n=nrow(table_outcome_df)))
descriptives_outcome <- delete_part(descriptives_outcome, part = "header")
descriptives_outcome <- add_header_lines(descriptives_outcome, values = c("Outcome Level"))
descriptives_outcome <- theme_vanilla(descriptives_outcome)
descriptives_outcome <- merge_v(descriptives_outcome, j = c("Category"))
myreport <- body_add_flextable(x = myreport, descriptives_outcome)
myreport <- body_add_par(x = myreport, value = "", style = "Normal")

########################
# MetaRegression Table
########################
MVnull.coef
str(MVnull.coef)
MVnull.coef$coef <- row.names(as.data.frame(MVnull.coef))
row.names(MVnull.coef) <- c()
MVnull.coef <- MVnull.coef[c("Coef", "beta", "SE", "tstat", "df_Satt", "p_Satt")]
MVnull.coef
str(MVnull.coef)

MVfull.coef$coef <- row.names(as.data.frame(MVfull.coef))
row.names(MVfull.coef) <- c()
MVfull.coef <- MVfull.coef[c("Coef", "beta", "SE", "tstat", "df_Satt", "p_Satt")]

# MetaRegression Table
model_null <- flextable(head(MVnull.coef, n=nrow(MVnull.coef)))
colkeys <- c("beta", "SE", "tstat", "df_Satt")
model_null <- colformat_double(model_null,  j = colkeys, digits = 2)
model_null <- colformat_double(model_null,  j = c("p_Satt"), digits = 3)
#model_null <- autofit(model_null)
model_null <- add_header_lines(model_null, values = c("Null Model"), top = FALSE)
model_null <- theme_vanilla(model_null)

myreport <- body_add_par(x = myreport, value = "Table 5: Model Results", style = "Normal")
myreport <- body_add_flextable(x = myreport, model_null)
#myreport <- body_add_par(x = myreport, value = "", style = "Normal")

model_full <- flextable(head(MVfull.coef, n=nrow(MVfull.coef)))
model_full <- colformat_double(model_full,  j = c("beta"), digits = 2)
model_full <- colformat_double(model_full,  j = c("p_Satt"), digits = 3)
#model_full <- autofit(model_full)
model_full <- delete_part(model_full, part = "header")
model_full <- add_header_lines(model_full, values = c("Meta-Regression"))
model_full <- theme_vanilla(model_full)

myreport <- body_add_flextable(x = myreport, model_full)
myreport <- body_add_par(x = myreport, value = "", style = "Normal")

# Marginal Means Table
marginalmeans <- flextable(head(means, n=nrow(means)))
colkeys <- c("moderator", "group", "SE", "tstat", "df")
marginalmeans <- colformat_double(marginalmeans,  j = colkeys, digits = 2)
marginalmeans <- colformat_double(marginalmeans,  j = c("p_Satt"), digits = 3)
rm(colkeys)
marginalmeans <- theme_vanilla(marginalmeans)
marginalmeans <- merge_v(marginalmeans, j = c("moderator"))
myreport <- body_add_par(x = myreport, value = "Table: Marginal Means", style = "Normal")
myreport <- body_add_flextable(x = myreport, marginalmeans)

# Write to word doc
file = paste("TableResults.docx", sep = "")
print(myreport, file)


#publication bias , selection modeling
full_y <- full$Effect.Size
full_v <- full$var
weightfunct(full_y, full_v)
weightfunct(full_y, full_v, steps = c(.025, .50, 1))

# source files for creating each visualization type
#source("/Users/apple/Desktop/SchoolbasedMentalHealth/Scripts/viz_MARC.R")

median(full$Sample.size)
sum(full$Sample.size)
# set sample sizes based on median sample size (331) and existing % weights for bar plot
#            w_j = meta-analytic weight for study j (before rescaling)
#            w_j_perc = percent weight allocated to study j

full <- full %>% 
  mutate(w_j = 1/(se^2)) %>%
  mutate(w_j_perc = w_j/sum(w_j)) %>%
  mutate(N_j = floor(w_j_perc*331))

#needed b/c tibbles throw error in viz_forest
full <- full %>% 
  ungroup()
full <- as.data.frame(full)

sd(full$Age.mean, na.rm=TRUE)
sd(full$Sample.size, na.rm=TRUE)
sd(full$Duration.weeks, na.rm=TRUE)

MVfull.modeling <- rma(yi=Effect.Size,
                       vi=var,
                       test="t",
                       data=full,
                       slab = Study,
                       method="REML")

#funnel plot
metafor::funnel(MVfull.modeling) 
reticulate::py_install("rispy")
sel <- selmodel(MVfull.modeling, type = "stepfun", steps = 0.025)
plot(sel, ylim=c(0,5))

jpeg(file="countour_funnel_color.jpeg")

# contour-enhanced funnel plot
metafor::funnel(MVfull.modeling, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=FALSE)
dev.off()

##############
#heatmap for outcomes
###########
library(ggplot2)
library(dplyr)
library(forcats)

# Ensure variables are factors
full <- full %>%
  mutate(
    Outcomes = as.factor(Outcomes),
    Study = as.factor(Study)
  )
by_study_outcome <- full %>%
  dplyr::group_by(Study, Outcomes) %>%
  dplyr::summarise(
    `Published Year` = first(as.integer(`Published.Year.f`)),
    `AI.Chatbot.Name` = first(`AI.Chatbot.Name`),
    Effect.Size = mean(Effect.Size, na.rm = TRUE),
    var = mean(var, na.rm = TRUE),
    n_rows = n()  # how many rows were averaged
  )
# Build ordered Study levels from the data frame (don’t overwrite it)
study_levels <- by_study_outcome %>%
  distinct(Study, `Published Year`) %>%      # one row per Study-year
  arrange(dplyr::desc(`Published Year`), Study) %>%       # older first, tie-break by name
  pull(Study)

# Apply the order to the Study column (choose rev(...) if you want oldest on top)
by_study_outcome <- by_study_outcome %>%
  mutate(Study = factor(Study, levels = study_levels))

ggplot(by_study_outcome, aes(x = Outcomes, y = Study, fill = Effect.Size)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F2F2F2", high = "#C53030", midpoint = 0,
    name = "Effect Size"
  ) +
  labs(
    x = "Mental Health Problems",
    y = "Study",
    title = "Heatmap of Effect Sizes by Outcome Category and Study"
  ) +
  geom_text(aes(label = sprintf("%.2f", Effect.Size)), size = 3) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

############
#traffic light
############
# Read the first sheet (or specify sheet = "Sheet1")
Rob <- read_excel("Rob.xlsx")

# Peek at your data
str(Rob)
names(Rob)
lapply(Rob, function(x) head(unique(x), 5))

library(dplyr)

standardize_judgement <- function(x) {
  x <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    x %in% c("y") ~ "Low",
    x %in% c("u") ~ "Some concerns",
    x %in% c("n") ~ "High"
  )
}

# Apply to all columns except Study (assumes first column is Study)
Rob_mapped <- Rob %>%
  dplyr::mutate(across(2:(ncol(Rob) - 1), standardize_judgement))

# Optional: check unique values after mapping
# lapply(Rob_mapped, function(col) sort(unique(col)))

# 2) Create traffic light plot using robvis in Generic mode
# If you have an Overall column already as last column, it will be included automatically below.
# If not, you can skip it or create one.

# Select only existing columns; ensure first column is Study (rename if needed)
if (names(Rob_mapped)[1] != "Study") names(Rob_mapped)[1] <- "Study"

# install.packages("robvis")
# library(robvis)
# 
# install.packages("remotes")
# remotes::install_github("mcguinlu/robvis")
# # Assuming rob_generic_df is your mapped JBI data
# # First column must be Study; middle columns are JBI items; last column (if any) is Overall.
# # Update robvis and dependencies
# install.packages("robvis", dependencies = TRUE)
# 
# # Restart R, then:
# library(robvis)
# packageVersion("robvis")          # confirm it updated
# getAnywhere("rob_traffic_light")  # check the function signature
# library(robvis)
# library(ggplot2)

# 2) Define how to categorize Overall Appraisal Score (set thresholds as appropriate)
map_overall <- function(x) {
  # Example thresholds: adjust to your rubric
  dplyr::case_when(
    x >= 11 ~ "Low",
    x >= 7 ~ "Some concerns",
    TRUE    ~ "High"
  )
}

# 3) Build plotting df with a mapped overall, but do NOT change your original last column
last_name <- names(Rob_mapped)[ncol(Rob_mapped)]  # "Overall Appraisal Score"
rob_generic_df <- Rob_mapped %>%
  mutate(Overall_mapped = map_overall(.data[[last_name]])) %>%
  select(-all_of(last_name), everything()) %>%     # drop numeric overall from plotting df
  relocate(Overall_mapped, .after = dplyr::last_col()) %>%
  dplyr::rename(Overall = Overall_mapped)
rob_generic_df <- rob_generic_df %>% dplyr::select(-`Overall Appraisal Score`)
# 4) Plot
p_traffic <- robvis::rob_traffic_light(rob_generic_df, tool = "Generic", psize = 12) +
  ggplot2::scale_fill_manual(values = c(
    "Low" = "#02C3BD",
    "Some concerns" = "#F5BE41",
    "High" = "#E45756"
  ), drop = FALSE, na.translate = FALSE)

n_studies <- nrow(rob_generic_df)
ggplot2::ggsave("traffic_light_JBI_RCT.png", p_traffic,
                width = 10, height = max(6, 0.35*n_studies), dpi = 300)

##########
#heatmap for social function
#########
# Ensure variables are factors (set Social.function order: task-oriented left, social-oriented right)
full <- full %>%
  mutate(
    Social.function = as.factor(Social.function),
    # If you want a specific ordering, uncomment and ensure labels match your data exactly:
    # Social.function = fct_relevel(Social.function, "task-oriented", "social-oriented"),
    Study = as.factor(Study)
  )

# Summarise by Study x Social.function
by_study_sf <- full %>%
  dplyr::group_by(Study, Social.function) %>%
  dplyr::summarise(
    `Published Year`   = first(as.integer(`Published.Year.f`)),
    `AI.Chatbot.Name`  = first(`AI.Chatbot.Name`),
    Effect.Size        = mean(Effect.Size, na.rm = TRUE),
    var                = mean(var, na.rm = TRUE),
    n_rows             = n(),                # how many rows were averaged
    .groups = "drop"
  )

# Build ordered Study levels (newest first; use rev(...) if you want oldest on top)
study_levels <- by_study_sf %>%
  dplyr::distinct(Study, `Published Year`) %>%
  dplyr::arrange(dplyr::desc(`Published Year`), Study) %>%
  dplyr::pull(Study)

# Apply the order to Study; also ensure Social.function order (task-oriented left, social-oriented right)
by_study_sf <- by_study_sf %>%
  mutate(
    Study = factor(Study, levels = study_levels),
    Social.function = fct_relevel(Social.function, "task-oriented", "social-oriented")
  )

# Plot
ggplot(by_study_sf, aes(x = Social.function, y = Study, fill = Effect.Size)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F2F2F2", high = "#C53030", midpoint = 0,
    name = "Effect Size"
  ) +
  labs(
    x = "Social Function Category",
    y = "Study",
    title = "Heatmap of Effect Sizes by Social Function and Study"
  ) +
  geom_text(aes(label = sprintf("%.2f", Effect.Size)), size = 3) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

