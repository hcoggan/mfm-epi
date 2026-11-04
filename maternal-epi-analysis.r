library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(testit)
library(httpgd)
library(broom)
library(stringr)
library(patchwork)
library(forcats)
library(pROC)
library(ggsignif)
library(boot)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggdendro)
library(car)
library(weights)
library(survey)
library(metafor)
library(pwrss)
library(svyVGAM)

library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)
library(MatchIt)
library(data.table)
library(exact2x2)
library(cobalt)
library(zipcodeR)
library(xgboost)
library(glmnet)
library(cowplot)
library(yardstick)
library(ggridges)
library(comorbidity)
library(readr)
library(fuzzyjoin)
library(WeightIt)

#Set random seed and working directory.
set.seed(210226)
load_filepath <- "/Users/helenacoggan/mnt/chip-lacava-e2/Groups/ch265493/mfm-penn/epi-work/mar-res-main/"
base_filepath <- "/Users/helenacoggan/mnt/chip-lacava-e2/Groups/ch265493/mfm-penn/epi-work/results/"

#-----ANALYSIS SETTINGS------
pH_threshold <- 7.1
base_excess_threshold <- Inf #Functionally no upper limit to the base excess
ref_population <- "White" #, 'Non-Black', 'Non-Black Except Unknown'
trim <- 0.0 #, 0.1, 0.05
truncate <- 0.05 #0.0, 0.05 Truncate at the PS-bound limits rather than trimming.
power_calc <- 0.8
doubly_robust <- FALSE
savename <- 'base'


#Implement these settings.
risk_metric <- 'Risk Ratio'
ps_bounds <- c(trim, 1-trim)
ps_bounds_truncate <- c(truncate, 1-truncate)


#Create directory if necessary.
dir.create(paste0(base_filepath, savename), showWarnings = FALSE)
setwd(paste0(base_filepath, savename))


###-----LOAD DATA-----

births <- read.csv(paste0(load_filepath, "preprocessed-births-with-outcomes.csv")) %>%
    rename(previous_birth_indicated = categorised_num_births) %>%
    mutate(acidemia = ifelse((p_h_cord <= pH_threshold) & (base_excess_cord <= base_excess_threshold), 1, 0)) %>% #Note- total labor length 80% missing
    mutate(birth_history = ifelse(
      previous_birth_indicated == 0, 'Nulliparous',
      ifelse(previous_c_section_indicated == 0, 'Previous births, no CD', 'Previous CD')
    )) %>% mutate(baby_mrn = as.character(baby_mrn))




#Combine early births to avoid sample size errors. After this, "preterm" means <37 weeks.
births$pregnancy_term <- ifelse(births$pregnancy_term == 'very preterm', 'preterm', births$pregnancy_term)

assert(nrow(births) == length(unique(births$baby_mrn))) #Passes

###----DEFINE CONFOUNDERS------
categorical_confounders <- c('year_of_delivery', #Include year but exclude seasonal trends (following Attanasio)
    'maternal_age', 'pregnancy_term')

binary_confounders <- c('previous_birth_indicated', 'previous_c_section_indicated')

#Create a dictionary.
confounders <- c(rep('categorical', length(categorical_confounders)), rep('binary', length(binary_confounders)))
names(confounders) <- c(categorical_confounders, binary_confounders)


###----HELPERS FOR DESCRIPTIVE TABLES------

# Chi-squared test when all EXPECTED counts are >= 5; otherwise Fisher's exact test.
# For tables larger than 2x2 the Fisher p-value is Monte Carlo simulated (reproducible
# via the seed above), since exact RxC computation can be very slow or run out of workspace.
test_independence <- function(x, y) {
  tab <- table(x, y)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (all(expected >= 5)) {
    chisq.test(tab)$p.value
  } else {
    fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value
  }
}

# N (%) of each covariate category, overall and within each level of group_var,
# plus a test of independence between the covariate and group_var.
describe_by <- function(data, group_var, covs) {
  out <- list()
  for (conf in covs) {
    overall <- data %>%
      count(Category = .data[[conf]], name = "N") %>%
      mutate(Perc = N / sum(N), group = "Overall")

    by_group <- data %>%
      count(group = .data[[group_var]], Category = .data[[conf]], name = "N") %>%
      group_by(group) %>%
      mutate(Perc = N / sum(N)) %>%   # fraction of births in group Y that are in category X
      ungroup()

    out[[conf]] <- bind_rows(
        overall  %>% mutate(Category = as.character(Category)),
        by_group %>% mutate(Category = as.character(Category), group = as.character(group))
      ) %>%
      mutate(Characteristic = conf,
             label = paste0(N, " (", round(100 * Perc, digits = 1), ")"),
             pval  = test_independence(data[[group_var]], data[[conf]]))
  }

  bind_rows(out) %>%
    tidyr::pivot_wider(id_cols = c(Characteristic, Category, pval),
                       names_from = group, values_from = label) %>%
    relocate(pval, .after = last_col())
}


#----CREATE COHORT TABLE (CHARACTERISTICS BY DELIVERY MODE)-------

biv_df <- describe_by(births, "delivery_method",
                      c('race', 'acidemia', categorical_confounders, binary_confounders))

characteristic_labels <- c(
  race                         = 'Race',
  maternal_age                 = 'Maternal Age',
  pregnancy_term               = 'Gestational Age (weeks)',
  previous_birth_indicated     = 'Previous Births',
  previous_c_section_indicated = 'Previous CD',
  year_of_delivery             = 'Year of Delivery',
  acidemia                     = 'Fetal Acidemia',
  delivery_method              = 'Delivery Method'
  #birth_history                = 'Obstetric Record'
)

# Names are the values in the data; values are the display labels.
category_labels <- c(
  "Asian" = "Asian", "Black" = "Black", "Other" = "Other", "Unknown" = "Unknown", "White" = "White",
  "Cesarean Delivery (FHT)" = "Cesarean Delivery (FHT)" ,
                    "Cesarean Delivery (Labor Arrest)" = "Cesarean Delivery (Labor Arrest)",
                    "Cesarean Delivery (Planned)"      = "Cesarean Delivery (Planned)" ,
                    "Vaginal" = "Vaginal",
  "0" = "None", "1" = "Any",
  "2017" = "2017", "2018" = "2018", "2019" = "2019", "2020" = "2020",
  "19 and under" = "19 and under", "20-29" = "20-29", "30-35" = "30-34", "35-39" = "35-39", "40+" = "40+",
  "preterm" = "<37", "early term" = "37-38", "full term" = "39-40", "late term" = "41+"
  #"Nulliparous" = "Nulliparous", "Previous births, no CD" = "Previous births, no CD", "Previous CD" = "Previous CD"
)

unlabelled <- setdiff(na.omit(unique(biv_df$Category)), names(category_labels))
if (length(unlabelled) > 0) {
  warning("Categories with no display label (they will show as NA): ",
          paste(unlabelled, collapse = ", "))
}

biv_df <- biv_df %>%
  mutate(Characteristic = factor(Characteristic, levels = names(characteristic_labels),
                                 labels = characteristic_labels),
         Category       = factor(Category, levels = names(category_labels),
                                 labels = category_labels)) %>%
  arrange(Characteristic, Category)

write.csv(biv_df, "association-with-delivery-method.csv", row.names = FALSE)


#----BIVARIATE ANALYSIS (DELIVERY METHOD AND ACIDEMIA BY RACE, AMONGST BLACK/REF PATIENTS)-------
biv_df_race <- describe_by(
  if  (ref_population == 'White'){
    births %>% filter(race %in% c("Black", "White"))
  } else if (ref_population == 'Non-Black') {
    births  %>% mutate(race=ifelse(race=='Black', 'Black', 'Non-Black'))                     #Include everyone
  } else if (ref_population == 'Non-Black Except Unknown') {
    births %>% filter(!(race == 'Unknown')) %>% mutate(race=ifelse(race=='Black', 'Black', 'Non-Black'))
  } else {
    stop(paste0("Unknown reference population: ", ref_population))
  }, "race",
                           c('delivery_method', categorical_confounders, binary_confounders, 'acidemia'))


unlabelled <- setdiff(na.omit(unique(biv_df_race$Category)), names(category_labels))
if (length(unlabelled) > 0) {
  warning("Categories with no display label (they will show as NA): ",
          paste(unlabelled, collapse = ", "))
}

biv_df_race <- biv_df_race %>%
  mutate(Characteristic = factor(Characteristic, levels = names(characteristic_labels),
                                 labels = characteristic_labels),
         Category       = factor(Category, levels = names(category_labels),
                                 labels = category_labels)) %>%
  arrange(Characteristic, Category)


write.csv(biv_df_race, "covariates-by-race.csv", row.names = FALSE)


################################################################################
# SETTINGS
################################################################################

covariate_vars <- names(confounders)   # defined above

ps_formula     <- as.formula(paste("race_black ~", paste(covariate_vars, collapse = " + ")))
outcome_covars <- paste(c("race_black", covariate_vars), collapse = " + ")

dm_short_names <- c("Cesarean Delivery (FHT)"          = "fht",
                    "Cesarean Delivery (Labor Arrest)" = "labor_arrest",
                    "Cesarean Delivery (Planned)"      = "planned_cd",
                    "Vaginal"                          = 'vaginal')

var_labels <- c(
  prop.score                   = "Propensity score",
  maternal_age                 = "Maternal age",
  pregnancy_term               = "Gestational age",
  previous_birth_indicated     = "Any previous birth",
  previous_c_section_indicated = "Any previous C-section",
  year_of_delivery             = "Year of delivery"
)

################################################################################
# HELPERS
################################################################################

relabel <- function(x) {
  for (v in names(var_labels)) {
    x <- gsub(paste0("\\b", v, "_"), paste0(var_labels[[v]], ": "), x)  # factor levels
    x <- gsub(paste0("\\b", v, "\\b"), var_labels[[v]], x)               # bare names
  }
  x
}

# Restrict to Black/White births (by default), convert covariates to factors, code treatment.
prep_data <- function(data) {
  if  (ref_population == 'White'){
    data <- data %>% filter(race %in% c("Black", "White"))
  } else if (ref_population == 'Non-Black') {
    data <- data                        #Include everyone
  } else if (ref_population == 'Non-Black Except Unknown') {
    data <- data %>% filter(!(race == 'Unknown'))
  } else {
    stop(paste0("Unknown reference population: ", ref_population))
  }
  for (v in covariate_vars) data[[v]] <- factor(data[[v]])
  data$race_black <- factor(ifelse(data$race == "Black", 1, 0), levels = c(0, 1))
  data$preterm <- factor(ifelse(data$pregnancy_term == 'preterm', 1, 0))
  data
}

# Truncate weights.
truncate <- function(vec, as_weights=FALSE) {
  lower_limit <- if (as_weights) 0 else ps_bounds_truncate[1]
  upper_limit <- if (as_weights) 1/ps_bounds_truncate[1] else ps_bounds_truncate[2]
  returned_vec <- ifelse(
    vec < lower_limit, lower_limit, ifelse(vec > upper_limit, upper_limit, vec)
  )
  return(returned_vec)

}

# Fit CBPS on the full sample, trim to the overlap region, refit on the trimmed sample.
#Also truncate weights where requested
fit_trimmed_weights <- function(data, label, subset_var=NULL) {
  fit_full <- weightit(ps_formula, data = data, method = "cbps", estimand = "ATE", by=subset_var)
  keep <- fit_full$ps >= ps_bounds[1] & fit_full$ps <= ps_bounds[2]
  births_retained <- sum(keep) / nrow(data)
 
  cat("\n====", label, "====\n")
  cat(sprintf("Units retained after trimming PS to [%.1f, %.1f]:\n", ps_bounds[1], ps_bounds[2]))
  print(table(Race = ifelse(data$race_black == 1, "Black", ref_population), Retained = keep))
  cat(sprintf("Fraction of births retained: %.1f%%\n", births_retained * 100))

  # droplevels() removes factor levels left empty by trimming, which would
  # otherwise create all-zero dummy columns in the CBPS and outcome models.
  data_trim <- droplevels(data[keep, ])
  fit_trim <- weightit(ps_formula, data = data_trim, method = "cbps", estimand = "ATE", by=subset_var)

  # Truncate propensity scores and weights if requested.
  fit_full$ps   <- truncate(fit_full$ps)
  fit_trim$ps   <- truncate(fit_trim$ps)
  fit_full$weights <- truncate(fit_full$weights, as_weights=TRUE)
  fit_trim$weights <- truncate(fit_trim$weights, as_weights=TRUE)

  # Add weights
  data_trim$iptw <- fit_trim$weights
 
  # Weights on the full cohort: IPTW for retained units, 0 for trimmed units.
  # Lets cobalt compare the full unweighted cohort against the trimmed weighted one.
  w_full <- numeric(nrow(data))
  w_full[keep] <- fit_trim$weights
 


  list(data = data_trim, fit = fit_trim, fit_full = fit_full, keep = keep,
       data_full = data, w_full = w_full, frac_retained = births_retained)
}



# Balance table and love plot: unadjusted = full cohort, adjusted = trimmed + IPTW.
# Takes the list returned by fit_trimmed_weights().
show_balance <- function(a, file, title, subset_var=NULL, subset_val=NULL) {

  #Pull out subsets if requested
  if(is.null(subset_var)) {
    subset_vec <- NULL
    app <- ""
  } else {
    subset_vec <- ifelse(a$data_full[[subset_var]] == subset_val, TRUE, FALSE)
    app <- paste0("_", subset_var, "_", subset_val)
  }

  bal <- bal.tab(ps_formula,
                 data       = a$data_full,
                 weights    = a$w_full,
                 method     = "weighting",
                 estimand   = "ATE",
                 stats      = "m",
                 thresholds = c(m = 0.1),
                 s.d.denom  = "pooled",   # pooled SD of the full unweighted cohort, used for both
                 binary     = "std",
                 un         = TRUE,
                 subset     = subset_vec)
  print(bal)
 
  old_names <- rownames(bal$Balance)
  new_names <- setNames(relabel(old_names), old_names)
 
  p <- love.plot(
    bal,
    stats        = "mean.diffs",
    abs          = TRUE,
    thresholds   = c(m = 0.1),
    var.order    = "unadjusted",
    sample.names = c("Unweighted (full cohort)", "IPTW (final cohort)"),
    colors       = c("#E69F00", "#0072B2"),
    shapes       = c("circle", "triangle"),
    title        = title,
    var.names    = new_names
  ) + theme(legend.position = "top")
 
  ggsave(paste0(file, app, ".pdf"), p, width = 7, height = 5)
  invisible(bal)
}

# Estimate for the race term: odds ratio for logistic models, rate ratio for poisson, mean difference for linear ones.
get_est <- function(fit, term = "race_black1", exponentiate = TRUE) {
  tf <- if (exponentiate) exp else identity
  ci <- confint(fit)[term, ]
  c(estimate = tf(unname(coef(fit)[term])),
    lower    = tf(unname(ci[1])),
    upper    = tf(unname(ci[2])),
    pval     = summary(fit)$coefficients[term, 4])
}

# Effect size for a binary event comparing Black vs White (by default), at the observed proportions
# and at 80% power. We assume R^2=0 in all of this for simplicity.
effect_for <- function(event, is_black, label, power = 0.8,
                      outcome_model = c("binomial", "poisson", "gaussian"),
                      alpha = 0.05) {
  outcome_model <- match.arg(outcome_model)

  is_black <- as.logical(is_black)
  treated_frac <- mean(is_black)
  control_prob <- mean(as.numeric(event[!is_black])) 
  n            <- length(event)
  sd.predictor <- sqrt(treated_frac*(1-treated_frac))
  sd.outcome   <- sd(as.numeric(event))



  args <- list(alpha = alpha, n = n,
               distribution = list(dist = "binomial", size = 1, prob = treated_frac),
               verbose = 0)

  #Binomial model (for acidemia)
  if (outcome_model == "binomial") {
    get_power     <- function(or) {
      pt <- do.call(power.z.logistic,
                      c(list(base.prob = control_prob, odds.ratio = or), args))
      pt$power - power} 
    
    #Get the size of the effect (assuming OR < 1/>1) we can detect at 80% power.
    low_effect <- uniroot(get_power, c(0.1, 1-1e-6), tol=1e-8)$root
    high_effect <- uniroot(get_power, c(1+1e-6, 100), tol=1e-8)$root
      
    metric <- "Odds ratio"

  #Poisson model
  } else if (outcome_model == 'poisson') {
    
    get_power     <- function(rr) {
      pt <- do.call(power.z.poisson,
                      c(list(base.rate = control_prob, rate.ratio = rr), args))
      pt$power - power} 
    
    #Get the size of the effect (assuming RR < 1/>1) we can detect at 80% power.
    low_effect <- uniroot(get_power, c(0.1, 1-1e-6), tol=1e-8)$root
    high_effect <- uniroot(get_power, c(1+1e-6, 100), tol=1e-8)$root
      
    metric <- "Rate ratio"

  #Gaussian model (for continuous differences)
  } else {
    get_power     <- function(beta) {
      pt <- do.call(power.t.regression, 
                      c(list(beta=beta,
                          r.squared = 0,
                          sd.predictor=sd.predictor, #Specified in the documentation to be necessary for binary predictors
                          sd.outcome = sd.outcome,
                          alpha = alpha, n = n, verbose = 0)))
      pt$power - power} 
    
    #Get the size of the effect (assuming difference < 0/>0) we can detect at 80% power.
    low_effect <- uniroot(get_power, c(-100, -1e-6), tol=1e-8)$root
    high_effect <- uniroot(get_power, c(1e-6, 100), tol=1e-8)$root
      
    metric <- "Difference"

  }



  data.frame(Outcome = label, SampleSize = n, Power = power,
             Metric = metric, LowEffectSize = low_effect, HighEffectSize = high_effect)
}
power_tests <- c()


################################################################################
# ANALYSIS 1: IPT-WEIGHTED POISSON REGRESSION OF RACE AND DELIVERY METHOD
################################################################################

a1 <- fit_trimmed_weights(prep_data(births), "Analysis 1: Delivery method")
df <- a1$data
df_full <- a1$data_full

show_balance(a1, "love-plot", "Covariate balance before and after weighting")


###---EFFECT SIZE TESTS (on the TRIMMED sample)---
######---Tests whether we can detect differences in rates of e.g. vaginal, CD for NRFHT births--
######---This is effectively for the one-vs-all model rather than full multinomial model but a rough guide--
for (dm in names(dm_short_names)) {
  df[[dm_short_names[[dm]]]] <- as.integer(df$delivery_method == dm)
  df_full[[dm_short_names[[dm]]]] <- as.integer(df_full$delivery_method == dm)
  power_tests[[dm]] <- effect_for(df$delivery_method == dm, df$race_black == 1, dm, outcome_model = "poisson", power = power_calc)
}



###---MULTINOMIAL MODEL------
###---Fit individual Poisson models first then subtract them----
###---Following procedure described in https://grodri.github.io/glms/notes/c6.pdf-----

#Pivot the table so we have 1 row per birth and delivery method.
pivot <- function(df, cluster = 'baby_mrn', covariates = names(confounders), use_weights=FALSE) {

  total_long_df <- data.frame()

  if (use_weights) {
    df_ <- df 
  } else {
    df_ <- df %>% mutate(iptw=1)
  }

  #One version of the table per model.
  for (dm in names(dm_short_names)) {
    long_df <- df_ %>%  mutate(outcome=ifelse(delivery_method==dm, 1, 0)) %>% 
      select(all_of(c('outcome', cluster, covariates, 'iptw', 'race_black'))) %>% 
      mutate(model_dm = dm)  #Highlight the dm relevant to the model; 'outcome' is 1 if the birth has that specific delivery method
    

    total_long_df <- rbind(total_long_df, long_df) 
  }
  #Add dummy columns for model.
  for (dm in names(dm_short_names)) {
    total_long_df[[paste0("model_dm_", dm_short_names[[dm]])]] <-
      as.integer(total_long_df$model_dm == dm)
  }
  return(total_long_df)

}

adj   <- svydesign(ids = ~ baby_mrn, weights = ~iptw, data = pivot(df, use_weights=TRUE))
unadj <- svydesign(ids = ~ baby_mrn, weights = ~iptw, data = pivot(df_full)) #Unadjusted ORs will be fit on the full cohort.

multinom_results <- c()

###----Fit single Poisson model-----
form <- as.formula("outcome ~ race_black + model_dm_fht + model_dm_planned_cd + model_dm_labor_arrest +
          model_dm_fht:race_black + model_dm_planned_cd:race_black + model_dm_labor_arrest:race_black ") 
          #'Case-specific parameter' can be reduced to 'race-specific' as probability must be identical within race groups
adj_fit <- svyglm(
  form,
  design = adj,
  family = quasipoisson(link = "log")
)
unadj_fit <- svyglm(
  form,
  design = unadj,
  family = quasipoisson(link = "log")
)

#Extract risk ratios for Black women relative to white women for each method relative to vaginal births
for (dm in names(dm_short_names)) {
  if (!(dm=='Vaginal')) {
    term_name <- paste0('race_black1:model_dm_', dm_short_names[[dm]])
    model_results <- data.frame(
      method = dm,
      model  = c("Adjusted", "Unadjusted"),
      rbind(get_est(adj_fit, term=term_name), get_est(unadj_fit, term=term_name)),
      row.names = NULL,
      frac_retained = round(a1$frac_retained*100, digits=1) #Fraction of births retained (after trimming)
    )

    multinom_results <- rbind(multinom_results, model_results)
  }
}

###---CORRECT AND SAVE---
multinom_results <- multinom_results %>%
  group_by(model) %>%
  mutate(pval = p.adjust(pval, method = "fdr")) %>%   # Benjamini-Hochberg
  ungroup()
write.csv(multinom_results, "multinomial-lr-of-delivery-method.csv", row.names = FALSE)

print(multinom_results)

###---PLOT---
method_order <- c("Vaginal",
                  "Cesarean Delivery (Planned)",
                  "Cesarean Delivery (Labor Arrest)",
                  "Cesarean Delivery (FHT)")

y_base <- setNames(rev(seq_along(method_order)), method_order)
offset <- 0.17

plot_df <- multinom_results %>%
  mutate(
    model = factor(model, levels = c("Unadjusted", "Adjusted")),
    y     = y_base[method] + ifelse(model == "Adjusted", offset, -offset),
    stars = case_when(pval < 0.001 ~ "***",
                      pval < 0.01  ~ "**",
                      pval < 0.05  ~ "*",
                      TRUE         ~ ""),
    label = paste0(sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper), stars)
  )

log_rng <- log(range(c(plot_df$lower, plot_df$upper, 1)))
brks    <- exp(seq(floor(log_rng[1] * 2) / 2, ceiling(log_rng[2] * 2) / 2, by = 0.5))
x_lim   <- range(brks)
text_x  <- x_lim[2] * 1.3
ref_y   <- y_base[["Vaginal"]]
head_y  <- max(y_base) + 0.6

cols <- c(Adjusted = "#1F4E79", Unadjusted = "#9A9A9A")

p <- ggplot(plot_df, aes(x = estimate, y = y, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 0.8) +
  geom_point(aes(shape = model), size = 3) +
  geom_text(aes(x = text_x, label = label), hjust = 0, size = 3.6,
            show.legend = FALSE) +
  annotate("point", x = 1, y = ref_y, shape = 18, size = 3.5, colour = "grey30") +
  annotate("text", x = text_x, y = ref_y, label = "Reference", hjust = 0,
           size = 3.6, colour = "grey30") +
  annotate("text", x = text_x, y = head_y, label = "Estimate [95% CI]",
           hjust = 0, size = 3.8, fontface = "bold") +
  scale_x_log10(breaks = brks,
                labels = function(b) sapply(round(b, 2), format),
                limits = x_lim,
                oob    = scales::oob_keep) +
  scale_y_continuous(breaks = y_base, labels = names(y_base),
                     expand = expansion(add = 0.5)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = c(Adjusted = 15, Unadjusted = 17)) +
  coord_cartesian(clip = "off") +
  labs(x = paste0(risk_metric, ": Black race relative to ", ref_population, " (log scale)"),
       y = NULL, colour = NULL, shape = NULL) +
       #caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-adjusted)") +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "top",
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin        = margin(10, 160, 10, 10),
        plot.caption       = element_text(colour = "grey40", size = 10))


ggsave("delivery-method-race-forest.pdf", p, width = 9.5, height = 5.5)

#################################################################################################
# ANALYSIS 1.5: EFFECT MODIFICATION OF RACE ON BIRTH HISTORY AND GESTATIONAL AGE FOR CD FOR NRFHT
#################################################################################################

###---Does race modify the effect of preterm birth, previous CD or previous birth on the likelihood of CD for NRFHT?----

stratification_vars <- c('preterm', 'previous_birth_indicated', 'previous_c_section_indicated')

###---Perform analysis of interaction-----
task_1_rrs <- data.frame()
task_2_rrs <- data.frame()
task_3_effects <- data.frame()
for (var in stratification_vars) {

  #Refit weights so that balance is achieved within strata
  a_sub <- fit_trimmed_weights(prep_data(births), "Analysis 2: Effect modification of race on ", subset_var = var)

  #Plot balance in both strata.
  for (val in c(1, 0)) {
    show_balance(a_sub, "love_group_strat", 'Covariate balancing within subsets', subset_var = var, subset_val = val) #Balance is near perfect in all strata
  }

  #Define data.
  df_sub <- pivot(a_sub$data, covariates=c("preterm", names(confounders)), use_weights = TRUE)
  df_sub$stratum <- factor(paste(ifelse(df_sub$race_black == 1, "Black", "White"),
          ifelse(df_sub[[var]] == 1, "pos", "neg"), sep = "_"),
    levels = c("White_neg", "White_pos", "Black_neg", "Black_pos"))  # first = reference-- e.g. White term, White nulliparous

  adj  <- svydesign(ids = ~ baby_mrn, weights = ~iptw, data = df_sub)

  #####TASK 1--- get all RRs relative to a single reference group (e.g. white term women)----.
  subform <- as.formula(paste0("outcome ~ stratum + model_dm_fht + model_dm_planned_cd + model_dm_labor_arrest +
          model_dm_fht:stratum + model_dm_planned_cd:stratum + model_dm_labor_arrest:stratum")) 

  #Fit single-outcome model.
  adj_fit <- svyglm(
    subform, 
    design = adj,
    family = quasipoisson(link = "log")
  )

  #Extract risk ratios for each outcome relative to vaginal birth-- for each stratum relative to white term women. 
  ref_group_rrs <- c()
  for (dm in names(dm_short_names)) {
    if (!(dm=='Vaginal')) {
      for (stratum in  c("White_pos", "Black_neg", "Black_pos")) {
        term_name <- paste0('stratum', stratum, ':model_dm_', dm_short_names[[dm]])
        model_results <- data.frame(
            method = dm,
            exposure = stratum,
            stratum_var = var,
            model  = "Adjusted", #We don't care about unadjusted results here
            rbind(get_est(adj_fit, term_name)),
            row.names = NULL,
            frac_retained = round(a_sub$frac_retained*100, digits=1) #Fraction of births retained (after trimming)
          )
        ref_group_rrs <- rbind(ref_group_rrs, model_results)
      }
      }
  }
  task_1_rrs <- rbind(task_1_rrs, ref_group_rrs)

  #####TASK 2--- get stratum-specific RRs for race alone----

  #Fit single-outcome model.
  pos_fit <- svyglm(
    form, #Defined in the previous analysis
    design = subset(adj, stratum %in% c("Black_pos", "White_pos")),
    family = quasipoisson(link = "log")
  )
  neg_fit <- svyglm(
    form,
    design = subset(adj, stratum %in% c("Black_neg", "White_neg")),
    family = quasipoisson(link = "log")
  )

  #Extract risk ratios for each outcome relative to vaginal birth-- for each stratum relative to white term women. 
  stratum_specific_rrs <- c()
  for (dm in names(dm_short_names)) {
    if (!(dm=='Vaginal')) {
      term_name <- paste0('race_black1:model_dm_', dm_short_names[[dm]])
      model_results <- data.frame(
          method = dm,
          stratum = c("pos", "neg"),
          stratum_var = var,
          model  = "Adjusted", #We don't care about unadjusted results here
          rbind(get_est(pos_fit, term_name), get_est(neg_fit, term_name)),
          row.names = NULL,
          frac_retained = round(a_sub$frac_retained*100, digits=1) #Fraction of births retained (after trimming)
        )
      stratum_specific_rrs <- rbind(stratum_specific_rrs, model_results)
      }
    }
    task_2_rrs <- rbind(task_2_rrs, stratum_specific_rrs)

    ###----TASK 3: Assess interaction on additive and multiplicative scale-----
    ###----Similarly to task 1, bootstrap R11, R10, ... from the balanced cohort----
    boot_stat <- function(births_df, idx, method = "fht") {
      d <- births_df[idx, ]
      d$baby_mrn <- paste0(d$baby_mrn, "_", seq_len(nrow(d)))      # unique ids for repeated births

      invisible(capture.output(w <- fit_trimmed_weights(d, "boot", subset_var = var)))  # re-estimate CBPS + trimming
      df_sub <- pivot(w$data, covariates=c(names(confounders), var), use_weights = TRUE)

      df_sub$stratum <- factor(paste(ifelse(df_sub$race_black == 1, "Black", "White"),
          ifelse(df_sub[[var]] == 1, "pos", "neg"), sep = "_"),
        levels = c("White_neg", "White_pos", "Black_neg", "Black_pos"))  # first = reference-- e.g. White term, White nulliparous

      fit <- glm(subform, data = df_sub, family = quasipoisson(link = "log"), weights = iptw)
      b <- coef(fit)

      #Get coefficients.
      log_rr <- function(s) {      
        unname(b[paste0("stratum", s, ":model_dm_", method)])
      }
      rr <- exp(c(log_rr("Black_pos"), log_rr("Black_neg"), log_rr("White_pos")))
      c(r11 = rr[1], r10 = rr[2], r01 = rr[3])
    }

    births_prepped <- prep_data(births)      # one row per birth, before weighting
    boot_results <- boot(births_prepped, boot_stat, R = 1000, method = "fht", parallel = "multicore", ncpus = 8)
    interaction_effects <- as.data.frame(boot_results$t)
    names(interaction_effects) <- names(boot_results$t0)
    interaction_effects <- interaction_effects %>% mutate(
      add = r11 + 1 - r10 - r01,
      mult = r11/(r10*r01)
    )

    summarise_interaction_effects <- function(est, reps, null) {
      reps <- reps[is.finite(reps)]                       # drop NA/Inf replicates
      ci   <- quantile(reps, c(0.025, 0.975))
      p    <- min(1, 2 * min(mean(reps <= null), mean(reps >= null)))
      c(est = unname(est), lc = unname(ci[1]), uc = unname(ci[2]),
        pval = p, n_used = length(reps))
    }

    #Get point estimates.
    t0 <- boot_results$t0
    est_add  <- t0["r11"] + 1 - t0["r10"] - t0["r01"]
    est_mult <- t0["r11"] / (t0["r10"] * t0["r01"])

    effects <- data.frame(
      stratification   = var,
      interaction_type = c("additive", "multiplicative"),
      rbind(summarise_interaction_effects(est_add,  interaction_effects$add,  null = 0),
            summarise_interaction_effects(est_mult, interaction_effects$mult, null = 1))
    )
    
    task_3_effects <- rbind(task_3_effects, effects)
}

#---CORRECT AND SAVE ANALYSES SEPARATELY----
task_1_rrs <- task_1_rrs %>% filter(method=="Cesarean Delivery (FHT)") %>%
  mutate(pval=p.adjust(pval, method='fdr'))
task_2_rrs <- task_2_rrs %>% filter(method=="Cesarean Delivery (FHT)") %>%
  mutate(pval=p.adjust(pval, method='fdr'))
task_3_effects <- task_3_effects %>%
  mutate(pval=p.adjust(pval, method='fdr'))
write.csv(task_1_rrs, "rr-between-strata-for-fht.csv")
write.csv(task_2_rrs, "stratum-specific-rr-for-fht.csv")
write.csv(task_3_effects, "interaction-effects-for-fht.csv")

#----PLOTTING CODE------

# Display labels for each stratification variable: panel title, and what "pos"/"neg" mean.
strat_labels <- list(
  preterm                      = c(title = "Gestational age", pos = "Preterm",     neg = "Term"),
  previous_birth_indicated     = c(title = "Previous births", pos = "Parous",      neg = "Nulliparous"),
  previous_c_section_indicated = c(title = "Previous CD",     pos = "Previous CD", neg = "No previous CD")
)
panel_title <- function(v)      unname(sapply(v, function(x) strat_labels[[x]][["title"]]))
level_label <- function(v, lvl) unname(mapply(function(x, l) strat_labels[[x]][[l]], v, lvl))
 
add_stars <- function(p) case_when(is.na(p)   ~ "",
                                   p < 0.001  ~ "***",
                                   p < 0.01   ~ "**",
                                   p < 0.05   ~ "*",
                                   TRUE       ~ "")
 
# Generic forest plot. `d` needs: label, panel, estimate, lower, upper, pval, is_ref.
# Rows are drawn top-to-bottom in the order they appear in `d`.
forest_plot <- function(d, xlab, null = 1, log_x = TRUE, caption = NULL) {
  d <- d %>% mutate(
    stars = ifelse(is_ref, "", add_stars(pval)),
    txt   = ifelse(is_ref, "Reference", sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)),
    label = factor(label, levels = rev(unique(label))),
    panel = factor(panel, levels = unique(panel)),
    col   = ifelse(is_ref, "grey30", "#1F4E79"),
    shp   = ifelse(is_ref, 18, 15)
  )
 
  rng <- range(c(d$lower, d$upper, null), na.rm = TRUE)
  if (log_x) {
    brks    <- exp(seq(floor(log(rng[1]) * 2) / 2, ceiling(log(rng[2]) * 2) / 2, by = 0.5))
    x_lim   <- range(brks)
    text_x  <- x_lim[2] * 1.3
    x_scale <- scale_x_log10(breaks = brks, labels = function(b) sapply(round(b, 2), format))
  } else {
    x_lim   <- rng + c(-1, 1) * 0.05 * diff(rng)
    text_x  <- x_lim[2] + 0.08 * diff(x_lim)
    x_scale <- scale_x_continuous()
  }
 
  ggplot(d, aes(x = estimate, y = label)) +
    geom_vline(xintercept = null, linetype = "dashed", colour = "grey50") +
    geom_linerange(aes(xmin = lower, xmax = upper, colour = col), linewidth = 0.8, na.rm = TRUE) +
    geom_point(aes(colour = col, shape = shp), size = 3) +
    geom_text(aes(x = upper, label = stars), hjust = -0.25, vjust = 0.75,
              size = 4, na.rm = TRUE) +
    geom_text(aes(x = text_x, label = txt), hjust = 0, size = 3.4) +
    scale_colour_identity() +
    scale_shape_identity() +
    x_scale +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    coord_cartesian(xlim = x_lim, clip = "off") +
    labs(x = xlab, y = NULL, caption = caption) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text         = element_text(face = "bold", hjust = 0),
          plot.margin        = margin(10, 150, 10, 10),
          plot.caption       = element_text(colour = "grey40", size = 9))
}
 
caption_wald <- "" #"* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-adjusted). CIs are unadjusted."
stratum_order <- c("White_neg", "White_pos", "Black_neg", "Black_pos")
 
 
###---TASK 1: Each race x stratum group vs the reference group (White, 'neg' level)---
refs <- data.frame(stratum_var = stratification_vars, exposure = "White_neg",
                   estimate = 1, lower = NA_real_, upper = NA_real_, pval = NA_real_,
                   is_ref = TRUE)
 
p1_df <- task_1_rrs %>%
  transmute(stratum_var, exposure, estimate, lower, upper, pval, is_ref = FALSE) %>%
  bind_rows(refs, .) %>%
  mutate(race  = sub("_.*", "", exposure),
         lvl   = sub(".*_", "", exposure),
         label = paste0(race, ", ", level_label(stratum_var, lvl)),
         panel = panel_title(stratum_var)) %>%
  arrange(match(stratum_var, stratification_vars), match(exposure, stratum_order))
 
p1 <- forest_plot(p1_df,
                  xlab    = "Relative risk ratio vs reference group: CD for NRFHT relative to vaginal birth (log scale)",
                  caption = caption_wald)
ggsave("fht-rr-between-strata.pdf", p1, width = 10, height = 8)
 
 
###---TASK 2: Black vs reference race within each stratum---
p2_df <- task_2_rrs %>%
  transmute(stratum_var,
            label  = level_label(stratum_var, stratum),
            panel  = panel_title(stratum_var),
            estimate, lower, upper, pval,
            is_ref = FALSE,
            ord    = match(stratum, c("neg", "pos"))) %>%
  arrange(match(stratum_var, stratification_vars), ord)
 
p2 <- forest_plot(p2_df,
                  xlab    = paste0("Relative risk ratio, Black vs ", ref_population,
                                   ": CD for NRFHT relative to vaginal birth (log scale)"),
                  caption = caption_wald)
ggsave("fht-stratum-specific-rr.pdf", p2, width = 10, height = 6)
 
 
###---TASK 3: Interaction measures (bootstrap)---
p3_df <- task_3_effects %>%
  transmute(interaction_type,
            label    = panel_title(stratification),
            panel    = ifelse(interaction_type == "additive", "Additive scale", "Multiplicative scale"),
            estimate = est, lower = lc, upper = uc, pval,
            is_ref   = FALSE)
 
p3_add  <- forest_plot(filter(p3_df, interaction_type == "additive"),
                       xlab = "Additive modification for the effect of Black race vs White",
                       null = 0, log_x = FALSE)
p3_mult <- forest_plot(filter(p3_df, interaction_type == "multiplicative"),
                       xlab = "Multiplicative modification for the effect of Black race vs White (log scale)",
                       null = 1, log_x = TRUE)
 
p3 <- p3_add + p3_mult +
  plot_annotation(#caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-adjusted). Bootstrap percentile CIs (R = 1000), unadjusted.",
                  theme = theme(plot.caption = element_text(colour = "grey40", size = 9)))
ggsave("fht-interaction-effects.pdf", p3, width = 14, height = 4.5)
 

################################################################################
# ANALYSIS 2: SUBGROUP ANALYSIS OF ACIDEMIA IN NRFHT BIRTHS
################################################################################

a2 <- fit_trimmed_weights(
  prep_data(births) %>% filter(delivery_method == "Cesarean Delivery (FHT)"),
  "Analysis 2: Acidemia in NRFHT births"
)
df2 <- a2$data

show_balance(a2, "love-plot-nrfht", "Covariate balance before and after weighting")

###---EFFECT SIZE (on the trimmed analytic sample)---
###---ACIDEMIA-------
power_tests[["acidemia_nrfht"]] <- effect_for(df2$acidemia == 1, df2$race_black == 1,
                                       "Acidemia (for NRFHT births)", outcome_model = 'binomial',  #Odds ratios are retained for acidemia
                                       power = power_calc)
###----PH DIFFERENCES-------
power_tests[["ph_nrfht"]] <- effect_for(df2$p_h_cord, df2$race_black == 1,
                                       "pH (for NRFHT births)", outcome_model = 'gaussian',  
                                       power = power_calc)

###---MODELS---
adj   <- svydesign(ids = ~1, weights = ~iptw, data = df2)
unadj <- svydesign(ids = ~1, weights = ~1,    data = a2$data_full)

results_list <- list()

for (y in c("acidemia", "depressed_apgar", "p_h_cord")) {
  is_continuous <- y == "p_h_cord"
  fam <- if (is_continuous) gaussian() else quasibinomial(link = "logit") #Fit binomial model here
  form <- if(doubly_robust) as.formula(paste(y, "~", outcome_covars)) else as.formula(paste(y, "~ race_black"))

  adj_fit   <- svyglm(form, design = adj,   family = fam)
  unadj_fit <- svyglm(as.formula(paste(y, "~ race_black")),      design = unadj, family = fam)

  results_list[[y]] <- data.frame(
    outcome = y,
    scale   = if (is_continuous) "Mean difference" else 'Odds Ratio',
    model   = c("Adjusted", "Unadjusted"),
    rbind(get_est(adj_fit,   exponentiate = !is_continuous),
          get_est(unadj_fit, exponentiate = !is_continuous)),
    row.names = NULL,
    frac_retained = round(a2$frac_retained*100, digits=1) #Fraction of births retained
  )
}

model_results <- do.call(rbind, results_list)
rownames(model_results) <- NULL

###--- SAVE---
# model_results <- model_results %>%
#   group_by(model) %>%
#   mutate(pval = p.adjust(pval, method = "fdr")) %>%   # Benjamini-Hochberg
#   ungroup()
write.csv(model_results, "acidemia-nrfht-analysis.csv", row.names = FALSE)
print(model_results)


################################################################################
# ANALYSIS 3: SUBGROUP ANALYSIS OF ACIDEMIA IN VAGINAL BIRTHS
################################################################################

a3 <- fit_trimmed_weights(
  prep_data(births) %>% filter(delivery_method == "Vaginal"),
  "Analysis 3: Acidemia in vaginal births"
)
df3 <- a3$data

show_balance(a3, "love-plot-vaginal", "Covariate balance before and after weighting")

###---EFFECT SIZE (on the trimmed analytic sample)---
power_tests[["acidemia_vaginal"]] <- effect_for(df3$acidemia == 1, df3$race_black == 1,
                                       "Acidemia (for vaginal births)", outcome_model = 'binomial',  #Odds ratios are retained for acidemia
                                       power = power_calc)
###----PH DIFFERENCES-------
power_tests[["ph_vaginal"]] <- effect_for(df3$p_h_cord, df3$race_black == 1,
                                       "pH (for vaginal births)", outcome_model = 'gaussian',  
                                       power = power_calc)


###---MODELS---
adj   <- svydesign(ids = ~1, weights = ~iptw, data = df3)
unadj <- svydesign(ids = ~1, weights = ~1,    data = a3$data_full)

results_list <- list()

for (y in c("acidemia", "depressed_apgar", "p_h_cord")) {
  is_continuous <- y == "p_h_cord"
  fam <- if (is_continuous) gaussian() else quasibinomial(link = "logit") #Fit binomial model here
  form <- if(doubly_robust) as.formula(paste(y, "~", outcome_covars)) else as.formula(paste(y, "~ race_black"))

  adj_fit   <- svyglm(form, design = adj,   family = fam)
  unadj_fit <- svyglm(as.formula(paste(y, "~ race_black")),      design = unadj, family = fam)

  results_list[[y]] <- data.frame(
    outcome = y,
    scale   = if (is_continuous) "Mean difference" else 'Odds Ratio',
    model   = c("Adjusted", "Unadjusted"),
    rbind(get_est(adj_fit,   exponentiate = !is_continuous),
          get_est(unadj_fit, exponentiate = !is_continuous)),
    row.names = NULL,
    frac_retained = round(a3$frac_retained*100, digits=1) #Fraction of births retained
  )
}

model_results <- do.call(rbind, results_list)
rownames(model_results) <- NULL

# ###---CORRECT AND SAVE---
# model_results <- model_results %>%
#   group_by(model) %>%
#   mutate(pval = p.adjust(pval, method = "fdr")) %>%   # Benjamini-Hochberg
#   ungroup()
write.csv(model_results, "acidemia-vaginal-analysis.csv", row.names = FALSE)
print(model_results)


################################################################################
# PLOT ACIDEMIA RESULTS
################################################################################

plot_acidemia_results <- function(outcome_name ='acidemia', savename="acidemia-forest-plot.pdf",
  label='Odds ratio', high_precision=FALSE, intercept=1) {
  results <- rbind(
    read.csv("acidemia-vaginal-analysis.csv") %>% mutate(DeliveryMethod = "Vaginal"),
    read.csv("acidemia-nrfht-analysis.csv")   %>% mutate(DeliveryMethod = "Cesarean Delivery (FHT)")
  ) %>%
    filter(outcome == outcome_name) %>%
    mutate(
      DeliveryMethod = factor(DeliveryMethod,
                              levels = c("Vaginal", "Cesarean Delivery (FHT)")),
      model = factor(model, levels = c("Unadjusted", "Adjusted")))
  if (high_precision) {
    results <- results %>% mutate(label = sprintf("%.3f [%.3f, %.3f]", estimate, lower, upper))
  } else {
     results <- results %>% mutate(label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper))
  } 
  results <- results %>%
    mutate(pval = p.adjust(pval, method = 'fdr'),
           stars = case_when(pval < 0.001 ~ "***",
                             pval < 0.01  ~ "**",
                             pval < 0.05  ~ "*",
                             TRUE         ~ ""),
           label = paste0(label, stars))

  dodge <- position_dodge(width = 0.5)
  print(results)
  
  p_acidemia <- ggplot(results,
                      aes(x = estimate, y = DeliveryMethod, colour = model, shape = model)) +
    geom_vline(xintercept = intercept, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lower, xmax = upper),
                  width = 0.15, orientation = "y", position = dodge) +
    geom_point(size = 2.5, position = dodge) +
    geom_text(aes(label = label), position = dodge, vjust = -1.1,
              size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c(Unadjusted = "grey55", Adjusted = "black")) +
    labs(
      title = "",
      x = paste0(label, ", Black vs ", ref_population),
      y = NULL, colour = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  if (intercept == 1) {
    p_acidemia <- p_acidemia + scale_x_log10() 
  }
  
  ggsave(savename, p_acidemia, width = 7, height = 4)
}

plot_acidemia_results(outcome_name='acidemia', savename="acidemia-forest-plot.pdf",
  label='Odds ratio (pH < 7.1) ')
plot_acidemia_results(outcome_name='p_h_cord', savename="pH-forest-plot.pdf",
  label='Mean pH difference', high_precision = TRUE, intercept = 0)

################################################################################
# PLOT pH DISTRIBUTION
################################################################################
pull_data <- function(a, delivery_method_label) {
  unweighted <- a$data_full %>%
    select(baby_mrn, p_h_cord, race, delivery_method) %>%
    mutate(Cohort = "Unweighted", iptw = 1)
  weighted <- a$data %>%
    select(baby_mrn, p_h_cord, iptw, race, delivery_method) %>%
    mutate(Cohort = "Weighted")
  rbind(unweighted, weighted)
}

plot_ph_distribution <- function(a2, a3, savename) {

  full_cohort <- rbind(pull_data(a2), pull_data(a3)) %>%
    mutate(
      race   = factor(race,   levels = c(ref_population, "Black")),
      Cohort = factor(Cohort, levels = c("Unweighted", "Weighted")),
      delivery_method = recode(delivery_method,                 # shorter facet labels
        "Cesarean Delivery (FHT)" = "CD for NRFHT",
        "Vaginal"                 = "Vaginal")
    )

  p <- ggplot(full_cohort, aes(x = p_h_cord, fill = race, colour = race, weight = iptw)) +
    geom_density(alpha = 0.25, linewidth = 0.6) +
    geom_vline(xintercept = pH_threshold, linetype = "dashed",   # mark acidemia cutoff
               colour = "grey40", linewidth = 0.5) +
    scale_fill_manual(values   = c("Black" = "#6BAED6", "White" = "#FD8D3C",
                                   "Non-Black" = "#FD8D3C",
                                   "Non-Black Except Unknown" = "#FD8D3C")) +
    scale_colour_manual(values = c("Black" = "#6BAED6", "White" = "#FD8D3C",
                                   "Non-Black" = "#FD8D3C",
                                   "Non-Black Except Unknown" = "#FD8D3C")) +
    facet_grid(Cohort ~ delivery_method) +                       # grid clearer than wrap
    labs(
      x      = "Umbilical cord gas pH",
      y      = "Probability density",
      fill   = "Race",
      colour = "Race"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.text         = element_text(face = "bold")
    )

  ggsave(savename, p, width = 8, height = 5)
}

plot_ph_distribution(a2, a3, savename = "ph_distribution_plot.pdf")
################################################################################
# SAVE POWER TESTS
################################################################################

write.csv(do.call(rbind, power_tests), "power_tests.csv", row.names = FALSE)