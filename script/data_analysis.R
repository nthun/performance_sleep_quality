library(tidyverse)
library(janitor)
library(skimr)
library(lme4)

# Data read
df_raw <- read_csv2("data/PSQIcikkhez adatok.csv", 
                    locale = locale(encoding = "WINDOWS-1252"))

# Define variable labels. This is important for making the plots pretty
variables_needed <- 
  c(
    "RT Triplet learning" = "trip_learn_all_rt",
    "ACC Triplet learning" = "trip_learn_all_acc",                 
    "ACC Higher-order sequence learning" = "highrer_order_acc",
    "RT Higher-order sequence learning" = "highrer_order_rt",
    "ACC Statistical learning" = "stat_learn_acc",
    "RT Statistical learning" = "stat_learn_rt",
    "Average ACC" = "acc_avg",
    "ACC general skill learning" = "acc_gs",
    "RT average" = "rt_avg",
    "RT general skill learning" = "rt_gs_1min4",
    "Counting Span" = "cs_avg",
    "WCST – perseverative error" = "wcst_pers_error",
    "Sleep disturbance" = "sleep_disturb",
    "study"
  )

variable_groups <- 
  tribble(~outcome_var, ~group,
          "RT Triplet learning", "RT learning indices",
          "ACC Triplet learning", "ACC learning indices",                 
          "ACC Higher-order sequence learning", "ACC learning indices",
          "RT Higher-order sequence learning", "RT learning indices",
          "ACC Statistical learning", "ACC learning indices",
          "RT Statistical learning", "RT learning indices",
          "Average ACC", "General skill indices",
          "ACC general skill learning", "General skill indices",
          "RT average", "General skill indices",
          "RT general skill learning", "General skill indices",
          "Counting Span", "WM and EF indices",
          "WCST – perseverative error", "WM and EF indices"
  )

# Process data for ease of use
df <-
  df_raw %>% 
  janitor::clean_names() %>% # Remove spaces, accented characters, capitals
  rename(sex = nem, age = eletkor, education_yrs = iskolazottsag) %>% 
  mutate(sex = recode(sex, `1` = "Male", `2` = "Female", `0` = "Not specified"))

# Principal component analysis of sleep quality
# This is required to make one PC that can represent sleep disturbance

sleepd_pca <- 
  df %>% 
  select(ais_ossz, psqi_osszpont) %>% 
  psych::pca(nfactors = 1)

# Add the PC to the original data 
df <- 
  df %>% 
  mutate(sleep_disturb = sleepd_pca$scores[,1])

# EDA
df %>% 
  select(one_of(variables_needed)) %>% 
  skimr::skim()

## Calculate association of the 
# First, create LMEM models with the same predictor (sleep disturbance), and various outcome measures. The LMEM contains a random intercept term by study

all_models <- 
  df %>% 
  select(one_of(variables_needed)) %>% 
  gather(outcome_var, value, -sleep_disturb, -study) %>% 
  split(.$outcome_var) %>% 
  map(., ~lmer(scale(value) ~ scale(sleep_disturb) + (1|study), REML = FALSE,
               data = .x))

# From the models, we can get the  standardized betas for the slopes and p the values
# We recalculate the models using restricted likelihood for this, so the estimated parameters will be more accurate

all_models %>% 
  map(~update(.x, REML = TRUE)) %>% 
  map(~broom::tidy(.)) %>% 
  bind_rows(.id = "outcome_var") %>% 
  filter(!str_detect(term, "(Intercept)|.Residual")) %>% 
  # Add number of observations by hand and calculate the p statistic baed on t value
  mutate(n = c(rep(231, 11), 225),
         p = 2*pt(q = -abs(statistic), df = n - 1))

## Calculate Bayes Factors (BF01)
# BF was calculated based on this tutoerial: https://rpubs.com/lindeloev/bayes_factors
# Calculate BIC for NULL models
# Maximum Likelihood estimation is used so the models can be compared

null_BIC <-
  df %>%
  select(one_of(variables_needed)) %>% 
  gather(outcome_var, value, -sleep_disturb, -study) %>% 
  split(.$outcome_var) %>%
  map_df(., ~lmer(scale(value) ~ 1 + (1|study), REML = FALSE, 
                  data = .x) %>% 
           broom::glance(.) %>% 
           select(nullBIC = BIC),
         .id = "outcome_var") 

# Calculate BF01 for all models
all_models %>% 
  map_df(~broom::glance(.),
         .id = "outcome_var") %>%
  left_join(null_BIC, by = "outcome_var") %>%
  mutate(BF_BIC = exp((BIC - nullBIC)/2))


# Visualize the association for each outcome measure
df %>%
  # Use pretty variable names
  select(variables_needed) %>%
  # Put data into long format so the outcome vriable can be faceted
  gather(outcome_var, value, -`Sleep disturbance`, -study) %>% 
  # Add variable group names
  left_join(variable_groups, by = "outcome_var") %>%
  ggplot() +
    aes(x = `Sleep disturbance`, y = value, color = as.factor(study)) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_grid(group~outcome_var, scales = "free", drop = TRUE) +
    scale_color_brewer(palette = "RdYlBu") +
    labs(y = NULL, color = "Study") +
    theme_bw()

