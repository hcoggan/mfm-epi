
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

#Set random seed and working directory.
set.seed(210226)
setwd("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/ch265493/mfm-penn/epi-work/mar-res-main/")
pH_threshold <- 7.10
base_excess_threshold <- Inf #Functionally no upper limit to the base excess


###-----LOAD DATA-----

births <- read.csv("preprocessed-births-with-outcomes.csv") %>%
    rename(previous_birth_indicated=categorised_num_births)

assert(nrow(births)==length(unique(births$baby_mrn))) #Passes

###----DEFINE CONFOUNDERS------
categorical_confounders <- c('year_of_delivery', #Include year but exclude seasonal trends  (following Attanasio)
    'maternal_age', 'pregnancy_term')

binary_confounders <- c('previous_c_section_indicated', 'previous_birth_indicated')

#Create a dictionary.
confounders <- c(rep('categorical', length(categorical_confounders)), rep('binary', length(binary_confounders)))
names(confounders) <- c(categorical_confounders, binary_confounders)


#----CREATE COHORT TABLE (CHARACTERISTICS BY DELIVERY MODE)-------

#Break race, acidemia and all confounders out by delivery method.
#Print the raw numbers and tests of bivariate association.
biv_df <- data.frame()
pvals <- c()
covs <- c('race', 'acidemia', categorical_confounders, binary_confounders)

for (conf in covs) {

    #Break births by confounder and delivery method.
    overall <- births %>%
        mutate(TotalBirths=n()) %>%
        group_by(.data[[conf]]) %>%
        summarise(N=n(), Perc=N/first(TotalBirths)) %>% #TotalBirths is a single value so 'first' just prevents datatype errors
        mutate(delivery_method='Overall')
    
    by_dm <- births %>%
        group_by(delivery_method) %>%
        mutate(TotalBirths=n()) %>% #By delivery method.
        group_by(delivery_method, .data[[conf]]) %>%
        summarise(N=n(), Perc=N/first(TotalBirths)) #i.e. the fraction of births of type Y which are of, for example, race X

    
    #Test for independence (#Perform chi-squared test.)
    tab <- table(births$delivery_method, births[[conf]])
    pval <- chisq.test(tab)$p.value

    #Now put all of this together.
    df <- rbind(overall, by_dm) %>% 
        rename(Category=all_of(conf)) %>% 
        mutate(Characteristic=conf)

    df$label <- paste0(
        df$N, " (", round(100*df$Perc, digits=1), ")"
    )
    pvals <- c(pvals, pval)
    biv_df <- rbind(biv_df, df)

}

#Pivot, format and return.
biv_df <- biv_df %>% tidyr::pivot_wider(
    id_cols=c('Characteristic', 'Category'),
    names_from='delivery_method',
    values_from = 'label'

)

#Attach pvalues.
names(pvals) <- covs
biv_df$pval <- unname(unlist(pvals[as.character(biv_df$Characteristic)]))

biv_df$Characteristic <- factor(
    biv_df$Characteristic, 
    levels = c('race', 'maternal_age', 'pregnancy_term',
         'previous_birth_indicated', 'previous_c_section_indicated', 
         'year_of_delivery', 'acidemia'),
    labels = c('Race', 'Maternal Age', 'Gestational Age (weeks)',
    'Previous Births', 'Previous CD', 'Year of Delivery', 'Fetal Acidemia')
)

biv_df$Category <- factor(biv_df$Category,
    levels = c("Asian", "Black", "Other", "Unknown", "White",
         "0", "1", "2017", "2018", "2019", "2020",
          "19 and under", "20-29", "30-35", "35-39",
           "40+", "very preterm", "preterm", "early term", "full term", "late term"),
    labels = c("Asian", "Black", "Other", "Unknown", "White",
         "None", "Any", "2017", "2018", "2019", "2020",
          "19 and under", "20-29", "30-34", "35-39",
           "40+", "<32", "32-36", "37-38", "39-40", "41+"))

biv_df <- biv_df %>%
    arrange(Characteristic, Category)

write.csv(biv_df, "association-with-delivery-method.csv")

#----CREATE OVERVIEW TABLE (DELIVERY METHOD AND ACIDEMIA BY RACE)-------

#Examine bivariate association of delivery method and acidemia with race.
biv_df <- data.frame()
pvals <- c()
covs <- c('delivery_method', 'acidemia')

for (conf in covs) {

    #Break births by outcome.
    overall <- births %>%
        mutate(TotalBirths=n()) %>%
        group_by(.data[[conf]]) %>%
        summarise(N=n(), Perc=N/first(TotalBirths)) %>% #
        mutate(race='Overall')
    
    by_race <- births %>%
        group_by(race) %>%
        mutate(TotalBirths=n()) %>% #By race
        group_by(race, .data[[conf]]) %>%
        summarise(N=n(), Perc=N/first(TotalBirths)) #i.e. the fraction of births of race Y which are of, for example, delivery method Y.

    
    #Test for independence (#Perform chi-squared test.)
    tab <- table(births$race, births[[conf]])
    pval <- chisq.test(tab)$p.value

    #Now put all of this together.
    df <- rbind(overall, by_race) %>% 
        rename(Category=all_of(conf)) %>% 
        mutate(Characteristic=conf)

    df$label <- paste0(
        df$N, " (", round(100*df$Perc, digits=1), ")"
    )
    pvals <- c(pvals, pval)
    biv_df <- rbind(biv_df, df)

}

#Pivot, format and return.
biv_df <- biv_df %>% tidyr::pivot_wider(
    id_cols=c('Characteristic', 'Category'),
    names_from='race',
    values_from = 'label'

)

#Attach pvalues.
names(pvals) <- covs
biv_df$pval <- unname(unlist(pvals[as.character(biv_df$Characteristic)]))

write.csv(biv_df, "dm-and-acidemia-by-race.csv")

#--------POWER ANALYSIS--------------

#File down the dataframe to just Black and White births.
df <- births %>% filter(race %in% c("Black", "White"))

#####-----DIFFERENCE IN LIKELIHOOD OF DELIVERY (RELATIVE TO VAGINAL BIRTHS)-----
for (dm in unique(df$delivery_method)){

    if (!(dm=='Vaginal')) {

        subdf <- df %>% filter(delivery_method %in% c(dm, 'Vaginal')) #73% power to detect difference in CD (Labor Arrest) relative to Vaginal birth; 76% to detect difference in likelihood of planned CD; ~100% powered to detect difference in likelihood of 

        #Get necessary proportions.
        treated_frac <- mean(subdf$race=='Black')
        treated_prob <- mean((subdf$delivery_method==dm) & (subdf$race=='Black'))/mean(subdf$race=='Black')
        control_prob <- mean((subdf$delivery_method==dm) & (subdf$race=='White'))/mean(subdf$race=='White')

        print(dm) #61% for Labor Arrest, 69% for planned CD, ~100% for planned
        power_test <- power.z.logistic(
                        base.prob=control_prob,
                        prob=treated_prob,
                        alpha=0.05,
                        n=nrow(subdf),
                        distribution = list(dist='binomial', size=1, prob=treated_frac))
    }

}


#####-----DIFFERENCE IN LIKELIHOOD OF ACIDEMIA (FOR NRFHT BIRTHS)-----

    subdf <- df %>% filter(delivery_method=="Cesarean Delivery (FHT)") #

    #Get necessary proportions.
    treated_frac <- mean(subdf$race=='Black')
    treated_prob <- mean((subdf$acidemia==1) & (subdf$race=='Black'))/mean(subdf$race=='Black')
    control_prob <- mean((subdf$acidemia==1) & (subdf$race=='White'))/mean(subdf$race=='White')

    print("Acidemia (for NRHFHT births)") #12.2
    power_test <- power.z.logistic(
                    base.prob=control_prob,
                    prob=treated_prob,
                    alpha=0.05,
                    n=nrow(subdf),
                    distribution = list(dist='binomial', size=1, prob=treated_frac))


#--------ANALYSIS 1: MULTINOMIAL IPT-WEIGHTED LOGISTIC REGRESSION ANALYSIS OF RACE AND DELIVERY METHOD-------
df <- births %>% filter(race %in% c("Black", "White"))

#Transform all covariates into factors.
for (cov in names(confounders)) {
    df[[cov]] <- factor(df[[cov]])
}

#Now fit propensity score.
df$race_black <- factor(ifelse(df$race=='Black', 1, 0))

ps_model <- glm(race_black ~ maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery,
        data=df, family=binomial(link='logit'))
df$propensity_score <- predict(ps_model, newdata = df, type = "response")


print('------PS MODEL SUMMARY--------')
print(summary(ps_model))

#Plot distribution of propensity scores.
p <- ggplot(df, aes(x = propensity_score, fill = race, color = race)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values = c(Black = "#E69F00", White = "#0072B2")) +
  scale_color_manual(values = c(Black = "#E69F00", White = "#0072B2")) +
  scale_x_continuous(limits = c(NA, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Propensity score distribution by race",
    x = "Propensity score",
    y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("ps-by-race.pdf", p, width = 7, height = 4.5)

###---ASSIGN IPT WEIGHTS, TRUNCATE AND STABILISE-----
df$race_black <- as.numeric(as.character(df$race_black))
treated_frac <- mean(df$race_black)

df$iptw <- ifelse(df$race_black, treated_frac/df$propensity_score, (1-treated_frac)/(1-df$propensity_score))
lower_limit <- quantile(df$iptw, 0.00)[[1]]
upper_limit <- quantile(df$iptw, 0.98)[[1]]
df$iptw <- ifelse(df$iptw < lower_limit, lower_limit, ifelse(df$iptw > upper_limit, upper_limit, df$iptw))

###----SHOW BALANCE---------
bal <- bal.tab(
  race_black ~ maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery,  
  data = df,
  weights = c("iptw"),
  method = "weighting",
  estimand = "ATE",
  s.d.denom = "pooled",
  un = TRUE,             # also show unweighted balance
  binary = "std"         # standardize binary covariates too
)

print(bal)  # table of SMDs before/after

var_labels <- c(
  maternal_age                 = "Maternal age",
  pregnancy_term               = "Gestational age",
  previous_birth_indicated     = "Any previous birth",
  previous_c_section_indicated = "Any previous C-section",
  year_of_delivery             = "Year of delivery"
)

relabel <- function(x) {
  for (v in names(var_labels)) {
    x <- gsub(paste0("\\b", v, "_"), paste0(var_labels[[v]], ": "), x)  # factor levels
    x <- gsub(paste0("\\b", v, "\\b"), var_labels[[v]], x)               # bare names
  }
  x
}

old_names <- rownames(bal$Balance)
new_names <- setNames(relabel(old_names), old_names)

love.plot(
  bal,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  sample.names = c("Unweighted", "IPTW"),
  colors = c("#E69F00", "#0072B2"),
  shapes = c("circle", "triangle"),
  title = "Covariate balance before and after weighting",
  var.names = new_names
) +
  ggplot2::theme(legend.position = "top")

ggplot2::ggsave("love-plot.pdf", width = 7, height = 5)

###-----MULTINOMIAL LOGISTIC REGRESSION -------

#Add adjusted and unadjusted ORs (include original specification to control for remaining biases) [https://pmc.ncbi.nlm.nih.gov/articles/PMC10087564/#bimj2388-sec-0140].


#########----Specify adjusted and unadjusted models---------
dm_short_names <- c("Cesarean Delivery (FHT)"          = "fht",
                    "Cesarean Delivery (Labor Arrest)" = "labor_arrest",
                    "Cesarean Delivery (Planned)"      = "planned_cd")

for (dm in names(dm_short_names)) {
  df[[dm_short_names[[dm]]]] <- ifelse(df$delivery_method == dm, 1, 0)
}

adj   <- svydesign(ids = ~1, weights = ~iptw, data = df)
unadj <- svydesign(ids = ~1, weights = ~1,    data = df)

get_or <- function(fit, term = "race_black") {
  ci <- confint(fit)[term, ]
  c(OR    = unname(exp(coef(fit)[term])),
    lower = unname(exp(ci[1])),
    upper = unname(exp(ci[2])),
    pval  = summary(fit)$coefficients[term, 4])
}
####-----DEFINE COVARIATES-----------
covars <- "race_black + maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery"


###------RUN ANALYSIS--------------
results_list <- list()

for (dm in names(dm_short_names)) {
  y <- dm_short_names[[dm]]

  adj_fit <- svyglm(
    as.formula(paste(y, "~", covars)),
    design = subset(adj, delivery_method %in% c(dm, "Vaginal")),
    #design = adj,
    family = quasibinomial(link='logit')
  )
  unadj_fit <- svyglm(
    as.formula(paste(y, "~ race_black")),
    design = subset(unadj, delivery_method %in% c(dm, "Vaginal")),
    #design = unadj,
    family = quasibinomial(link='logit')
  )

  results_list[[dm]] <- data.frame(
    method = dm,
    model  = c("Adjusted", "Unadjusted"),
    rbind(get_or(adj_fit), get_or(unadj_fit)),
    row.names = NULL
  )
}

model_results <- do.call(rbind, results_list)
rownames(model_results) <- NULL

###-----CORRECT AND SAVE----------
model_results <- model_results %>% 
        group_by(model) %>%
        mutate(pval=p.adjust(pval, method='fdr')) #Benjamini-Hochberg correction
write.csv(model_results, "multinomial-lr-of-delivery-method.csv")

####------PLOT----------

# Rows from top to bottom. Vaginal is the reference category (no estimate).
method_order <- c("Vaginal",
                  "Cesarean Delivery (Planned)",
                  "Cesarean Delivery (Labor Arrest)",
                  "Cesarean Delivery (FHT)")

y_base <- setNames(rev(seq_along(method_order)), method_order)
offset <- 0.17   # vertical gap between adjusted (above) and unadjusted (below)

plot_df <- model_results %>%
  mutate(
    model = factor(model, levels = c("Unadjusted", "Adjusted")),
    y     = y_base[method] + ifelse(model == "Adjusted", offset, -offset),
    stars = case_when(pval < 0.001 ~ "***",
                      pval < 0.01  ~ "**",
                      pval < 0.05  ~ "*",
                      TRUE         ~ ""),
    label = sprintf("%.2f [%.2f, %.2f]", OR, lower, upper)
  )

# Log-scale breaks at half-unit steps on the log scale (0.61, 1, 1.65, 2.72, ...),
# widened automatically to cover all CIs.
log_rng <- log(range(c(plot_df$lower, plot_df$upper, 1)))
brks    <- exp(seq(floor(log_rng[1] * 2) / 2, ceiling(log_rng[2] * 2) / 2, by = 0.5))
x_lim   <- range(brks)
text_x  <- x_lim[2] * 1.3   # where the estimate column starts (right of the panel)
ref_y   <- y_base[["Vaginal"]]
head_y  <- max(y_base) + 0.6

cols <- c(Adjusted = "#1F4E79", Unadjusted = "#9A9A9A")

p <- ggplot(plot_df, aes(x = OR, y = y, colour = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 0.8) +
  geom_point(aes(shape = model), size = 3) +
  # significance stars just right of each CI
  geom_text(aes(x = upper, label = stars), hjust = -0.25, vjust = 0.75,
            size = 4, show.legend = FALSE) +
  # estimate column
  geom_text(aes(x = text_x, label = label), hjust = 0, size = 3.6,
            show.legend = FALSE) +
  # reference row
  annotate("point", x = 1, y = ref_y, shape = 18, size = 3.5, colour = "grey30") +  # diamond
  annotate("text", x = text_x, y = ref_y, label = "Reference", hjust = 0,
           size = 3.6, colour = "grey30") +
  # column header
  annotate("text", x = text_x, y = head_y, label = "Estimate [95% CI]",
           hjust = 0, size = 3.8, fontface = "bold") +
  scale_x_log10(breaks = brks,
                labels = function(b) sapply(round(b, 2), format),
                limits = x_lim,
                oob = scales::oob_keep) +
  scale_y_continuous(breaks = y_base, labels = names(y_base),
                     expand = expansion(add = 0.5)) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = c(Adjusted = 15, Unadjusted = 17)) +  # square, triangle
  coord_cartesian(clip = "off") +
  labs(x = "Odds ratio: Black race relative to White (log scale)",
       y = NULL, colour = NULL, shape = NULL,
       caption = "* p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_minimal(base_size = 13) +
  theme(legend.position     = "top",
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin        = margin(10, 160, 10, 10),
        plot.caption       = element_text(colour = "grey40", size = 10))

ggsave("delivery-method-race-forest.pdf", p, width = 9.5, height = 5.5)

#------ANALYSIS 2: SUBGROUP ANALYSIS OF ACIDEMIA IN NRFHT BIRTHS-----------------

df <- births %>% filter(race %in% c("Black", "White"), delivery_method=='Cesarean Delivery (FHT)')

#Transform all covariates into factors.
for (cov in names(confounders)) {
    df[[cov]] <- factor(df[[cov]])
}

#Now fit propensity score.
df$race_black <- factor(ifelse(df$race=='Black', 1, 0))

ps_model <- glm(race_black ~ maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery,
        data=df, family=binomial(link='logit'))
df$propensity_score <- predict(ps_model, newdata = df, type = "response")


print('------PS MODEL SUMMARY--------')
print(summary(ps_model))

#Plot distribution of propensity scores.
p <- ggplot(df, aes(x = propensity_score, fill = race, color = race)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values = c(Black = "#E69F00", White = "#0072B2")) +
  scale_color_manual(values = c(Black = "#E69F00", White = "#0072B2")) +
  scale_x_continuous(limits = c(NA, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    title = "Propensity score distribution by race",
    x = "Propensity score",
    y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave("ps-by-race-in-subgroup.pdf", p, width = 7, height = 4.5)

###---ASSIGN IPT WEIGHTS, TRUNCATE AND STABILISE-----
df$race_black <- as.numeric(as.character(df$race_black))
treated_frac <- mean(df$race_black)

df$iptw <- ifelse(df$race_black, treated_frac/df$propensity_score, (1-treated_frac)/(1-df$propensity_score))
lower_limit <- quantile(df$iptw, 0.0)[[1]]
upper_limit <- quantile(df$iptw, 0.98)[[1]]
df$iptw <- ifelse(df$iptw < lower_limit, lower_limit, ifelse(df$iptw > upper_limit, upper_limit, df$iptw))

###----SHOW BALANCE---------
bal <- bal.tab(
  race_black ~ maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery,  
  data = df,
  weights = c("iptw"),
  method = "weighting",
  estimand = "ATE",
  s.d.denom = "pooled",
  un = TRUE,             # also show unweighted balance
  binary = "std"         # standardize binary covariates too
)

print(bal)  # table of SMDs before/after

var_labels <- c(
  maternal_age                 = "Maternal age",
  pregnancy_term               = "Gestational age",
  previous_birth_indicated     = "Any previous birth",
  previous_c_section_indicated = "Any previous C-section",
  year_of_delivery             = "Year of delivery"
)

relabel <- function(x) {
  for (v in names(var_labels)) {
    x <- gsub(paste0("\\b", v, "_"), paste0(var_labels[[v]], ": "), x)  # factor levels
    x <- gsub(paste0("\\b", v, "\\b"), var_labels[[v]], x)               # bare names
  }
  x
}

old_names <- rownames(bal$Balance)
new_names <- setNames(relabel(old_names), old_names)

love.plot(
  bal,
  stats = "mean.diffs",
  abs = TRUE,
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  sample.names = c("Unweighted", "IPTW"),
  colors = c("#E69F00", "#0072B2"),
  shapes = c("circle", "triangle"),
  title = "Covariate balance before and after weighting",
  var.names = new_names
) +
  ggplot2::theme(legend.position = "top")

ggplot2::ggsave("love-plot-subgroup.pdf", width = 7, height = 5)

###-----MULTINOMIAL LOGISTIC REGRESSION -------

#Add adjusted and unadjusted ORs (include original specification to control for remaining biases) [https://pmc.ncbi.nlm.nih.gov/articles/PMC10087564/#bimj2388-sec-0140].


#########----Specify adjusted and unadjusted models---------
adj   <- svydesign(ids = ~1, weights = ~iptw, data = df)
unadj <- svydesign(ids = ~1, weights = ~1,    data = df)

####-----DEFINE COVARIATES-----------
covars <- "race_black + maternal_age + pregnancy_term + previous_birth_indicated + previous_c_section_indicated + year_of_delivery"


###------RUN ANALYSIS--------------
results_list <- list()

for (y in c('acidemia')) {

  adj_fit <- svyglm(
    as.formula(paste(y, "~", covars)),
    design = adj,
    family = quasibinomial(link='logit')
  )
  unadj_fit <- svyglm(
    as.formula(paste(y, "~ race_black")),
    design = unadj,
    family = quasibinomial(link='logit')
  )

  results_list[[y]] <- data.frame(
    outcome = y,
    model  = c("Adjusted", "Unadjusted"),
    rbind(get_or(adj_fit), get_or(unadj_fit)),
    row.names = NULL
  )
}

model_results <- do.call(rbind, results_list)
rownames(model_results) <- NULL

###-----CORRECT AND SAVE----------
model_results <- model_results %>% 
        group_by(model) %>%
        mutate(pval=p.adjust(pval, method='fdr')) #Benjamini-Hochberg correction
write.csv(model_results, "acidemia-subgroup-analysis.csv")

# print(model_results)