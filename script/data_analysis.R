install.packages(c("tidyverse", "janitor","psych","lme4"))


library(tidyverse)
library(lme4)
library(cowplot)

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
  mutate(BF01 = exp((BIC - nullBIC)/2))


# Visualize the association for each outcome measure
library(patchwork)

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
  theme(legend.direction = "vertical",
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





