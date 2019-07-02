install.packages(c("tidyverse", "janitor","psych","lme4"))


library(MASS)
library(tidyverse)
library(lme4)
library(cowplot)
library(broom.mixed)

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

bartlett_result <- 
  df %>% 
  select(ais_ossz, psqi_osszpont) %>% 
  cor() %>% 
  psych::cortest.bartlett(n = 250)

# Show PCA result
sleepd_pca

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

# all_models <- 
#   df %>% 
#   select(one_of(variables_needed)) %>% 
#   gather(outcome_var, value, -sleep_disturb, -study) %>% 
#   split(.$outcome_var) %>% 
#   map(., ~lmer(scale(value) ~ scale(sleep_disturb) + (1|study), REML = FALSE,
#                data = .x))

# Roboust lmer
all_models <- 
  df %>% 
  select(one_of(variables_needed)) %>% 
  gather(outcome_var, value, -sleep_disturb, -study) %>% 
  split(.$outcome_var) %>% 
  map(., ~robustlmm::rlmer(scale(value) ~ scale(sleep_disturb) + (scale(sleep_disturb)|study), REML = FALSE,
               data = .x))

all_models[[1]] %>% summary()


# From the models, we can get the  standardized betas for the slopes and p the values
# We recalculate the models using restricted likelihood for this, so the estimated parameters will be more accurate

all_models %>% 
  map(~update(.x, REML = TRUE)) %>% 
  map(~tidy(.)) %>% 
  bind_rows(.id = "outcome_var") %>% 
  filter(term == "scale(sleep_disturb)") %>% 
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
  mutate(BF01 = exp((BIC - nullBIC)/2))


# Visualize the association for each outcome measure

plot_df <-
  df %>%
  # Use pretty variable names
  select(variables_needed) %>% 
  # Put data into long format so the outcome vriable can be faceted
  gather(outcome_var, value, -`Sleep disturbance`, -study) %>% 
  # Add variable group names
  left_join(variable_groups, by = "outcome_var")

# Create a plotting function for all groups separately
plot_groups <- function(group_name){
  plot_df %>%
    filter(group == group_name) %>% 
    ggplot() +
    aes(x = `Sleep disturbance`, y = value, color = as.factor(study)) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_brewer(palette = "Paired") +
    facet_grid(group~outcome_var, scales = "free_y", switch = "y") + 
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
}

# Put all rows of the plot into a list
group_plots <- 
  map(unique(plot_df$group), plot_groups)

# Plot all panels separately --------------------------------------------------------

# Create a plotting function for all groups separately
plot_single <- function(outcome_name){
  plot_df %>%
    filter(outcome_var == outcome_name) %>% 
    ggplot() +
    aes(x = `Sleep disturbance`, y = value, color = as.factor(study)) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    scale_color_brewer(palette = "Paired") +
    facet_wrap(~outcome_var, scales = "free_y") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    cowplot::panel_border()
}

# Put all rows of the plot into a list
single_plots <- 
  map(unique(plot_df$outcome_var), plot_single)

# Create one plot to get the legend
plot_1 <-
  plot_df %>%
  filter(outcome_var == "RT Triplet learning") %>% 
  ggplot() +
  aes(x = `Sleep disturbance`, y = value, color = as.factor(study)) +
  geom_point(alpha = .5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_brewer(palette = "Paired") +
  labs(color = "Study") +
  theme(legend.direction = "horizontal",
        legend.justification = "right",
        legend.box.just = "top", 
        legend.box.margin = margin())

# Get legend
plot_legend <- cowplot::get_legend(plot_1)

plot_grid(
  NULL, NULL, NULL, NULL,
  single_plots[[1]], single_plots[[4]], single_plots[[6]], NULL,
  NULL, NULL, NULL, NULL,
  single_plots[[2]], single_plots[[3]], single_plots[[5]], NULL,
  NULL, NULL, NULL, NULL,
  single_plots[[7]], single_plots[[8]], single_plots[[9]], single_plots[[10]],
  NULL, NULL, NULL, NULL,
  single_plots[[11]], single_plots[[12]], NULL, plot_legend, 
  labels = c(
    "", "", "","",
    "RT learning indices", "", "", "",
    "", "", "","",
    "ACC learning indices", "", "", "",
    "", "", "","",
    "General skill indices", "", "", "",
    "", "", "","",
    "WM and EF indices"),
  hjust = -.05,
  vjust = .05,
  nrow = 8,
  rel_heights = c(.15, 1, .15, 1, .15, 1, .15, 1)
  ) 

# After adjustments
temp <-
plot_grid(
  NULL, NULL, NULL, NULL,
  single_plots[[1]], single_plots[[2]], single_plots[[7]], single_plots[[9]],
  NULL, NULL, NULL, NULL,
  single_plots[[4]], single_plots[[3]], single_plots[[8]], single_plots[[10]],
  NULL, NULL, NULL, NULL,
  single_plots[[6]], single_plots[[5]], single_plots[[11]], single_plots[[12]],
  plot_legend, NULL, NULL, NULL,
  labels = c(
    "\n\n\nRT learning indices", "\n\n\nACC learning indices", "\n\n\nGeneral skill indices", "",
    "", "", "","",
    "", "", "","",
    "", "", "", "",
    "", "", "\n\n\nWM and EF indices","",
    "", "", "","",
    "", "", "",""),
  hjust = -.15,
  vjust = .25,
  nrow = 7,
  rel_heights = c(.2, 1, .26, 1, .2, 1, .2)
) 

temp + 
  geom_rect(xmin = 0.003, ymin = 0.05, xmax = .245, ymax = .995, alpha = 0, color = "black", size = .5) +
  geom_rect(xmin = 0.25, ymin = 0.05, xmax = .50, ymax = .995, alpha = 0, color = "black", size = .5) +
  geom_rect(xmin = 0.505, ymin = 0.365, xmax = .997, ymax = .995, alpha = 0, color = "black", size = .5) +
  geom_rect(xmin = 0.505, ymin = 0.05, xmax = .997, ymax = .355, alpha = 0, color = "black", size = .5) +
  NULL

ggdraw()

# geom_rect(xmin = 0.005, ymin = 0.05, xmax = .50, ymax = .99, alpha = 0, color = "black", linetype = "dashed", size = 1)



# Version 3 -------------------------------------------------------------------------

plot_grid(
  NULL, NULL, plot_legend,
  single_plots[[1]], single_plots[[4]], single_plots[[6]], 
  NULL, NULL, NULL, 
  single_plots[[2]], single_plots[[3]], single_plots[[5]], 
  NULL, NULL, NULL, 
  single_plots[[9]], single_plots[[10]],single_plots[[11]], 
  NULL, NULL, NULL, 
  single_plots[[7]], single_plots[[8]], single_plots[[12]],
  labels = c(
    "", "", "",
    "A) RT learning indices", "", "",
    "", "", "",
    "B) ACC learning indices", "", "", 
    "", "", "",
    "C) General skill indices", "", "D) WM and EF indices",
    "", "", "",
    "", "", ""),
  hjust = -.05,
  vjust = .05,
  nrow = 8,
  rel_heights = c(.20, 1, .15, 1, .15, 1, .15, 1)
) + 
  geom_segment(y = .74, yend = .74, x = .02, xend = .98, size = 1.1) +
  geom_segment(y = .49, yend = .49, x = .02, xend = .98, size = 1.1) +
  geom_segment(y = .02, yend = .49, x = .668, xend = .668, size = 1.1) +
  NULL



# Check residuals -------------------------------------------------------------------
model_resid <-
  all_models %>% 
  map(~update(.x, REML = TRUE) %>% 
      broom::augment(.) %>% 
      select(.resid)) %>%
      bind_rows(.id = "outcome_var")


model_resid %>% 
  ggplot() +
  aes(x = .resid) +
  geom_histogram() +
  facet_wrap(~outcome_var)

model_resid %>% 
  ggplot() +
  aes(sample = .resid) +
  geom_qq_line() +
  geom_qq() +
  facet_wrap(~outcome_var)


model_resid %>% 
  group_by(outcome_var) %>% 
  nest() %>% 
  mutate(shapiro = map(data, ~shapiro.test(.x$.resid) %>% 
                               broom::tidy())) %>% 
  unnest(shapiro) %>% 
  View()


variables_arrange <-
  tibble(variables_needed,
       outcome_var = names(variables_needed)) %>% 
  slice(-(13:14)) %>% 
  left_join(variable_groups) %>% 
  arrange(variables_needed)


write_tsv(variables_arrange, "clipboard")

df %>% 
  select(ais_ossz, psqi_osszpont, sleep_disturb) %>% 
  rename("Sleep disturbance" = sleep_disturb,
         "AIS" = ais_ossz, 
         "PSQI" = psqi_osszpont) %>% 
  gather(measure, value) %>% 
  ggplot() +
  aes(x = value) +
  geom_histogram() +
  facet_wrap(~measure, scales = "free")
  
  
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
  
  variables_recode <-
    names(variables_needed) %>% 
    set_names(variables_needed)
  
  models %>% 
    mutate(outcome_var = recode(outcome_var, !!!variables_recode)) %>% 
    select(-group, -effect, -term) %>% 
    left_join(variable_groups, by = "outcome_var") %>% 
    arrange(group) %>% 
    select(group, outcome_var, everything()) %>% 
    write_excel_csv("model_summary.csv")
  
  df %>% 
    select(one_of(variables_needed)) %>% 
    gather(outcome_var, value, -sleep_disturb, -study) %>% 
    split(.$outcome_var) %>% 
    map(., ~lmer(scale(value) ~ scale(sleep_disturb) + (1|study), REML = FALSE,
                 data = .x)) %>% 
    map_df(~broom::glance(.),
           .id = "outcome_var") %>%
    left_join(null_BIC, by = "outcome_var") %>%
    mutate(BF01 = exp((BIC - nullBIC)/2)) %>% 
    select(outcome_var, BF01) %>% 
    mutate(outcome_var = recode(outcome_var, !!!variables_recode)) %>% 
    left_join(variable_groups, by = "outcome_var") %>% 
    arrange(group) %>% 
    select(group, outcome_var, BF01) %>% 
    write_excel_csv("model_bf.csv")
  
  
  glimpse(df)

  df %>% 
    select(study, groningen_ossz, nb_ossz) %>% View()
  
  df %>%
    rename(diary_psqi_5_f = naplo_psqi_5_f) %>% 
    select(-diary_psqi_4) %>% 
    mutate(diary = rowSums(select(., contains("diary")))) %>% 
    select(diary, groningen_ossz, variables_needed, -study) %>% 
    cor(use = "pairwise.complete.obs") %>% 
    View()
    
    View()
    
    nrow()
    summary()
  
  
  
  df %>% 
    drop_na(groningen_ossz) %>% 
    nrow()
  
  df %>% 
    filter(study == 2) %>% 
    glimpse()
    

# GSQ models ------------------------------------------------------------------------

  study2_models <- 
    df %>% 
    filter(study == 2) %>% 
    select(one_of(variables_needed), groningen_ossz, -study, -sleep_disturb) %>% 
    gather(outcome_var, value, -groningen_ossz) %>% 
    split(.$outcome_var) %>% 
    map(., ~lm(scale(rank(value)) ~ scale(groningen_ossz),
                  data = .x))
  
  sudy2_models_df <-
    study2_models %>% 
    map(~tidy(.)) %>% 
    bind_rows(.id = "outcome_var") %>% 
    filter(term == "scale(groningen_ossz)") %>% 
    # Add number of observations by hand and calculate the two-sided p statistic 
    # based on t value and residual df
    mutate(df = rep(103, 12)-1)
  
  sudy2_models_df %>% 
    mutate(outcome_var = recode(outcome_var, !!!variables_recode)) %>% 
    left_join(variable_groups, by = "outcome_var") %>% 
    arrange(group) %>% 
    select(group, outcome_var:statistic, df, p.value) %>% 
    knitr::kable(digits = 4)
  
  
  null_study2_BIC <-
    df %>%
    select(one_of(variables_needed), groningen_ossz,  -study, -sleep_disturb) %>% 
    gather(outcome_var, value, -groningen_ossz) %>% 
    split(.$outcome_var) %>%
    map_df(., ~rlm(scale(rank(value)) ~ 1,
                    data = .x) %>% 
             glance(.) %>% 
             select(nullBIC = BIC),
           .id = "outcome_var") 
  
  study2_models %>% 
    map_df(~broom::glance(.),
           .id = "outcome_var") %>%
    left_join(null_study2_BIC, by = "outcome_var") %>%
    mutate(BF01 = exp((BIC - nullBIC)/2)) %>% 
    select(outcome_var, BF01) %>% 
    mutate(outcome_var = recode(outcome_var, !!!variables_recode)) %>% 
    left_join(variable_groups, by = "outcome_var") %>% 
    arrange(group) %>% 
    select(group, outcome_var, BF01)
  
  
# This will be good!  
study2_long <-
  df %>% 
  filter(study == 2) %>%
  select(one_of(variables_needed), -one_of(sleep_variables), -sleep_disturb, -study,
         groningen_ossz, diary_psqi_4_f) %>% 
  gather(outcome_var, out_value, -groningen_ossz, -diary_psqi_4_f, na.rm = TRUE) %>% 
  gather(predictor_var, pred_value, groningen_ossz, diary_psqi_4_f) %>% 
  drop_na(pred_value, out_value) %>% 
  group_nest(outcome_var, predictor_var)

# Get null BICs
study2_nullBIC <-
  study2_long %>% 
  mutate(model = map(data, ~rlm(scale(out_value) ~ 1,
                                data = .x) %>% 
                       glance(.) %>% 
                       select(nullBIC = BIC))) %>% 
  unnest(model) %>% 
  select(-data)
  
  
study2_models <-
  study2_long %>% 
  mutate(model = map(data, ~rlm(scale(out_value) ~ scale(pred_value),
                               data = .x)),
         n = map_int(data, ~nrow(.x)))

# Show predictors
study2_models %>% 
  mutate(model = map(model, ~tidy(.x))) %>% 
  unnest(model) %>% 
  filter(str_detect(term, "scale")) %>% 
  mutate_at(vars(outcome_var, predictor_var),
            ~recode(., !!!variables_recode)) %>% 
  left_join(variable_groups, ., by = "outcome_var") %>% 
  select(group, everything(), -term) %>% 
  arrange(predictor_var, group, outcome_var)
  
# Calculate BF
study2_models %>% 
  mutate(model = map(model, ~glance(.x))) %>% 
  unnest(model) %>% 
  left_join(study2_nullBIC, by = c("outcome_var", "predictor_var")) %>% 
  mutate(BF01 = exp((BIC - nullBIC)/2) %>% round(4)) %>% 
  mutate_at(vars(outcome_var, predictor_var),
            ~recode(., !!!variables_recode)) %>% 
  left_join(variable_groups, ., by = "outcome_var") %>% 
  select(group, outcome_var, predictor_var, n, BIC, nullBIC, BF01) %>%
  arrange(predictor_var, group, outcome_var)



# PLOT STUDY 2 ----------------------------------------------------------------------
# GSQI

plot2_df <-
  study2_long %>% 
  mutate_at(vars(outcome_var, predictor_var),
            ~recode(., !!!variables_recode)) %>% 
  unnest()

# Create a plotting function for all groups separately
plot_double <- function(outcome_name){
  plot2_df %>%
    filter(outcome_var == outcome_name) %>% 
    mutate_at(vars(outcome_var, predictor_var),
              ~recode(., !!!variables_recode)) %>% 
    ggplot() +
    aes(x = pred_value, y = out_value, color = predictor_var) +
    geom_point(alpha = .5) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~outcome_var, scales = "free_y") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    cowplot::panel_border()
}

# Put all rows of the plot into a list
double_plots <-
  plot2_df %>%
  group_nest(outcome_var, keep = TRUE) %>% 
  mutate(plot = map(data, ~plot_double(.x$outcome_var))) %>% 
  pull(plot)

double_plots[[2]]

# Create one plot to get the legend
plot_2 <-
  plot2_df %>%
  filter(outcome_var == "RT Triplet learning") %>% 
  ggplot() +
  aes(x = pred_value, y = out_value, color = predictor_var) +
  geom_point(alpha = .5) +
  geom_smooth(method = "lm", se = FALSE) +
  # scale_color_viridis_d() +
  labs(color = "Predictor") +
  theme(legend.direction = "horizontal",
        legend.justification = "right",
        legend.box.just = "top", 
        legend.box.margin = margin())

# Get legend
plot_legend <- cowplot::get_legend(plot_2)

# Assemble the plot. 
# I realize that the code is a bit ugly, but the outcome seems acceptable
plot_grid(
  NULL, plot_legend, NULL, 
  double_plots[[1]], double_plots[[4]], double_plots[[6]], 
  NULL, NULL, NULL, 
  double_plots[[2]], double_plots[[3]], double_plots[[5]], 
  NULL, NULL, NULL, 
  double_plots[[9]], double_plots[[10]],double_plots[[11]], 
  NULL, NULL, NULL, 
  double_plots[[7]], double_plots[[8]], double_plots[[12]],
  labels = c(
    "", "", "",
    "A) RT learning indices", "", "",
    "", "", "",
    "B) ACC learning indices", "", "", 
    "", "", "",
    "C) General skill indices", "", "D) WM and EF indices",
    "", "", "",
    "", "", ""),
  hjust = -.05,
  vjust = .05,
  nrow = 8,
  rel_heights = c(.20, 1, .15, 1, .15, 1, .15, 1)
) + 
  geom_segment(y = .74, yend = .74, x = .02, xend = .98, size = 1.1) +
  geom_segment(y = .49, yend = .49, x = .02, xend = .98, size = 1.1) +
  geom_segment(y = .02, yend = .49, x = .668, xend = .668, size = 1.1) +
  NULL

nrow(all_models[[1]]@frame)

  