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
multinom_model <- "poisson" #, "binomial"
trim <- 0.1 #, 0.0, 0.05
power_calc <- 0.8
doubly_robust <- FALSE
savename <- 'base'


#Implement these settings.
risk_metric <- if (multinom_model == "binomial") {
  "Odds Ratio"
} else if (multinom_model == "poisson") {
  "Risk Ratio"
} else {
  stop("Unknown model: ", multinom_model)
}
ps_bounds <- c(trim, 1-trim)


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
    ))




#Combine early births to avoid sample size errors. After this, "preterm" means <37 weeks.
births$pregnancy_term <- ifelse(births$pregnancy_term == 'very preterm', 'preterm', births$pregnancy_term)

assert(nrow(births) == length(unique(births$baby_mrn))) #Passes

###----DEFINE CONFOUNDERS------
categorical_confounders <- c('year_of_delivery', #Include year but exclude seasonal trends (following Attanasio)
    'maternal_age', 'pregnancy_term', 'birth_history')

binary_confounders <- c()

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
  # previous_birth_indicated     = 'Previous Births',
  # previous_c_section_indicated = 'Previous CD',
  year_of_delivery             = 'Year of Delivery',
  acidemia                     = 'Fetal Acidemia',
  delivery_method              = 'Delivery Method',
  birth_history                = 'Obstetric Record'
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
  "preterm" = "<37", "early term" = "37-38", "full term" = "39-40", "late term" = "41+",
  "Nulliparous" = "Nulliparous", "Previous births, no CD" = "Previous births, no CD", "Previous CD" = "Previous CD"
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
    births  %>% mutate(race==ifelse(race=='Black', 'Black', 'Non-Black'))                     #Include everyone
  } else if (ref_population == 'Non-Black Except Unknown') {
    births %>% filter(!(race == 'Unknown')) %>% mutate(race==ifelse(race=='Black', 'Black', 'Non-Black'))
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
                    "Cesarean Delivery (Planned)"      = "planned_cd")

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
  data
}


# Fit CBPS on the full sample, trim to the overlap region, refit on the trimmed sample.
fit_trimmed_weights <- function(data, label) {
  fit_full <- weightit(ps_formula, data = data, method = "cbps", estimand = "ATE")
  keep <- fit_full$ps >= ps_bounds[1] & fit_full$ps <= ps_bounds[2]
  births_retained <- sum(keep) / nrow(data)
 
  cat("\n====", label, "====\n")
  cat(sprintf("Units retained after trimming PS to [%.1f, %.1f]:\n", ps_bounds[1], ps_bounds[2]))
  print(table(Race = ifelse(data$race_black == 1, "Black", ref_population), Retained = keep))
  cat(sprintf("Fraction of births retained: %.1f%%\n", births_retained * 100))
 
  # droplevels() removes factor levels left empty by trimming, which would
  # otherwise create all-zero dummy columns in the CBPS and outcome models.
  data_trim <- droplevels(data[keep, ])
 
  fit_trim <- weightit(ps_formula, data = data_trim, method = "cbps", estimand = "ATE")
  data_trim$iptw <- fit_trim$weights
 
  # Weights on the full cohort: IPTW for retained units, 0 for trimmed units.
  # Lets cobalt compare the full unweighted cohort against the trimmed weighted one.
  w_full <- numeric(nrow(data))
  w_full[keep] <- fit_trim$weights
 
  print(summary(fit_trim))
  list(data = data_trim, fit = fit_trim, fit_full = fit_full, keep = keep,
       data_full = data, w_full = w_full, frac_retained = births_retained)
}
 
# Balance table and love plot: unadjusted = full cohort, adjusted = trimmed + IPTW.
# Takes the list returned by fit_trimmed_weights().
show_balance <- function(a, file, title) {

  bal <- bal.tab(ps_formula,
                 data       = a$data_full,
                 weights    = a$w_full,
                 method     = "weighting",
                 estimand   = "ATE",
                 stats      = "m",
                 thresholds = c(m = 0.1),
                 s.d.denom  = "pooled",   # pooled SD of the full unweighted cohort, used for both
                 binary     = "std",
                 un         = TRUE)
  print(bal)
 
  old_names <- rownames(bal$Balance)
  new_names <- setNames(relabel(old_names), old_names)
 
  p <- love.plot(
    bal,
    stats        = "mean.diffs",
    abs          = TRUE,
    thresholds   = c(m = 0.1),
    var.order    = "unadjusted",
    sample.names = c("Unweighted (full cohort)", "IPTW (trimmed)"),
    colors       = c("#E69F00", "#0072B2"),
    shapes       = c("circle", "triangle"),
    title        = title,
    var.names    = new_names
  ) + theme(legend.position = "top")
 
  ggsave(file, p, width = 7, height = 5)
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
# and at 80% power.
effect_for <- function(event, is_black, label, power = 0.8,
                      outcome_model = c("binomial", "poisson"),
                      alpha = 0.05) {
  outcome_model <- match.arg(outcome_model)

  is_black <- as.logical(is_black)
  treated_frac <- mean(is_black)
  control_prob <- mean(event[!is_black])
  n            <- length(event)


  args <- list(alpha = alpha, n = n,
               distribution = list(dist = "binomial", size = 1, prob = treated_frac),
               verbose = 0)

  if (outcome_model == "binomial") {
    get_power     <- function(or) {
      pt <- do.call(power.z.logistic,
                      c(list(base.prob = control_prob, odds.ratio = or), args))
      pt$power - power} 
    
    #Get the size of the effect (assuming OR < 1/>1) we can detect at 80% power.
    low_effect <- uniroot(get_power, c(0.1, 1-1e-6), tol=1e-8)$root
    high_effect <- uniroot(get_power, c(1+1e-6, 100), tol=1e-8)$root
      
    metric <- "Odds ratio"
  } else {
    
    get_power     <- function(rr) {
      pt <- do.call(power.z.poisson,
                      c(list(base.rate = control_prob, rate.ratio = rr), args))
      pt$power - power} 
    
    #Get the size of the effect (assuming RR < 1/>1) we can detect at 80% power.
    low_effect <- uniroot(get_power, c(0.1, 1-1e-6), tol=1e-8)$root
    high_effect <- uniroot(get_power, c(1+1e-6, 100), tol=1e-8)$root
      
    metric <- "Rate ratio"
  }

  data.frame(Outcome = label, SampleSize = n, Power = power,
             Metric = metric, LowEffectSize = low_effect, HighEffectSize = high_effect)
}
power_tests <- c()


################################################################################
# ANALYSIS 1: IPT-WEIGHTED LOGISTIC REGRESSION OF RACE AND DELIVERY METHOD
################################################################################

a1 <- fit_trimmed_weights(prep_data(births), "Analysis 1: Delivery method")
df <- a1$data
df_full <- a1$data_full

show_balance(a1, "love-plot.pdf", "Covariate balance before and after weighting")



###---EFFECT SIZE TESTS (on the TRIMMED sample)---
for (dm in names(dm_short_names)) {
  df[[dm_short_names[[dm]]]] <- as.integer(df$delivery_method == dm)
  df_full[[dm_short_names[[dm]]]] <- as.integer(df_full$delivery_method == dm)

  subdf <- df %>% filter(delivery_method %in% c(dm, "Vaginal"))
  power_tests[[dm]] <- effect_for(subdf$delivery_method == dm, subdf$race_black == 1, dm, outcome_model = multinom_model, power = power_calc)
}

###---MODELS---
# Each delivery method is compared with vaginal birth among births with one of
# those two outcomes, which reproduces the multinomial (baseline-category) ORs.
adj   <- svydesign(ids = ~1, weights = ~iptw, data = df)
unadj <- svydesign(ids = ~1, weights = ~1,    data = df_full) #Unadjusted ORs will be fit on the full cohort.
model_family <- if (multinom_model == "binomial") {
  quasibinomial(link = "logit")
} else if (multinom_model == "poisson") {
  quasipoisson(link = "log")
} else {
  stop("Unknown model: ", multinom_model)
}

results_list <- list()


for (dm in names(dm_short_names)) {
  y <- dm_short_names[[dm]]
  form <- if(doubly_robust) as.formula(paste(y, "~", outcome_covars)) else as.formula(paste(y, "~ race_black"))



  adj_fit <- svyglm(
    form,
    design = subset(adj, delivery_method %in% c(dm, "Vaginal")),
    family = model_family
  )
  unadj_fit <- svyglm(
    form,
    design = subset(unadj, delivery_method %in% c(dm, "Vaginal")),
    family = model_family
  )

  results_list[[dm]] <- data.frame(
    method = dm,
    model  = c("Adjusted", "Unadjusted"),
    rbind(get_est(adj_fit), get_est(unadj_fit)),
    row.names = NULL,
    frac_retained = round(a1$frac_retained*100, digits=1) #Fraction of births retained (after trimming)
  )
}
model_results <- do.call(rbind, results_list)
rownames(model_results) <- NULL

###---CORRECT AND SAVE---
model_results <- model_results %>%
  group_by(model) %>%
  mutate(pval = p.adjust(pval, method = "fdr")) %>%   # Benjamini-Hochberg
  ungroup()
write.csv(model_results, "multinomial-lr-of-delivery-method.csv", row.names = FALSE)

###---PLOT---
method_order <- c("Vaginal",
                  "Cesarean Delivery (Planned)",
                  "Cesarean Delivery (Labor Arrest)",
                  "Cesarean Delivery (FHT)")

y_base <- setNames(rev(seq_along(method_order)), method_order)
offset <- 0.17

plot_df <- model_results %>%
  mutate(
    model = factor(model, levels = c("Unadjusted", "Adjusted")),
    y     = y_base[method] + ifelse(model == "Adjusted", offset, -offset),
    stars = case_when(pval < 0.001 ~ "***",
                      pval < 0.01  ~ "**",
                      pval < 0.05  ~ "*",
                      TRUE         ~ ""),
    label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)
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
  geom_text(aes(x = upper, label = stars), hjust = -0.25, vjust = 0.75,
            size = 4, show.legend = FALSE) +
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
       y = NULL, colour = NULL, shape = NULL,
       caption = "* p < 0.05, ** p < 0.01, *** p < 0.001 (FDR-adjusted)") +
  theme_minimal(base_size = 13) +
  theme(legend.position    = "top",
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin        = margin(10, 160, 10, 10),
        plot.caption       = element_text(colour = "grey40", size = 10))

ggsave("delivery-method-race-forest.pdf", p, width = 9.5, height = 5.5)

################################################################################
# ANALYSIS 2: SUBGROUP ANALYSIS OF ACIDEMIA IN NRFHT BIRTHS
################################################################################

a2 <- fit_trimmed_weights(
  prep_data(births) %>% filter(delivery_method == "Cesarean Delivery (FHT)"),
  "Analysis 2: Acidemia in NRFHT births"
)
df2 <- a2$data

show_balance(a2, "love-plot-nrfht.pdf", "Covariate balance before and after weighting")

###---EFFECT SIZE (on the trimmed analytic sample)---
power_tests[["acidemia"]] <- effect_for(df2$acidemia == 1, df2$race_black == 1,
                                       "Acidemia (for NRFHT births)", outcome_model = 'binomial',  #Odds ratios are retained for acidemia
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

show_balance(a3, "love-plot-vaginal.pdf", "Covariate balance before and after weighting")

###---EFFECT SIZE (on the trimmed analytic sample)---
power_tests[["acidemia"]] <- effect_for(df3$acidemia == 1, df3$race_black == 1,
                                       "Acidemia (for vaginal births)", outcome_model = 'binomial',  #Odds ratios are retained for acidemia
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


results <- rbind(
  read.csv("acidemia-vaginal-analysis.csv") %>% mutate(DeliveryMethod = "Vaginal"),
  read.csv("acidemia-nrfht-analysis.csv")   %>% mutate(DeliveryMethod = "Cesarean Delivery (FHT)")
) %>%
  filter(outcome == "acidemia") %>%     # outcome values are lowercase in the saved CSVs
  mutate(
    DeliveryMethod = factor(DeliveryMethod,
                            levels = c("Vaginal", "Cesarean Delivery (FHT)")),
    model = factor(model, levels = c("Unadjusted", "Adjusted")),
    label = sprintf("%.2f [%.2f, %.2f]", estimate, lower, upper)
  )
 
dodge <- position_dodge(width = 0.5)
 
p_acidemia <- ggplot(results,
                     aes(x = estimate, y = DeliveryMethod, colour = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                width = 0.15, orientation = "y", position = dodge) +
  geom_point(size = 2.5, position = dodge) +
  geom_text(aes(label = label), position = dodge, vjust = -1.1,
            size = 3, show.legend = FALSE) +
  scale_x_log10() +
  scale_colour_manual(values = c(Unadjusted = "grey55", Adjusted = "black")) +
  labs(
    title = "",
    x = "Odds ratio, Black vs White (log scale)",
    y = NULL, colour = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
 
ggsave("acidemia-forest-plot.pdf", p_acidemia, width = 7, height = 4)

 

################################################################################
# SAVE POWER TESTS
################################################################################

write.csv(do.call(rbind, power_tests), "power_tests.csv", row.names = FALSE)