
library(lubridate)
library(jsonlite)
library(tidycmprsk)
library(testit)
library(broom)
library(stringr)
library(weights)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(janitor)
library(MatchIt)
library(data.table)
library(tidyverse)
library(xgboost)
library(pROC)

#Set random seed and working directory.
set.seed(210226)
setwd("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/ch265493/mfm-penn/epi-work/preprocessing/")
load_filepath <- "/Users/helenacoggan/Documents/OSXLAP12460/mfm/new-dataset/"


########-----LOAD DATASETS-----------------------------------

raw_deliveries  <- fread(paste0(load_filepath, "Deliveries Dataset 2026.11.01.csv")) 
deliveries  <- fread("preprocessed-deliveries-without-acidemia-outcomes.csv") %>%
    inner_join(select(raw_deliveries, c(PregnancyId, NewbornId, DeliveryMrn, BirthInstant)),
    by = c("PregnancyId", "NewbornId", "DeliveryMrn"))

# labs        <- read.csv(paste0(load_filepath, "Labs Dataset 2026.11.01.csv"))
# pregnancies <- fread(paste0(load_filepath, "Pregnancies Dataset 2026.11.01.csv"))
newborns    <- fread(paste0(load_filepath, "Newborns Dataset 2026.11.01.csv"))


#------GET DIFFERENCES BETWEEN CORD GAS MEASUREMENTS IN THE SAME PATIENT----
outcomes <- newborns %>% rename(DeliveryMrn=MomDeliveryMrn) 
delivery_times <- deliveries %>% select(BirthInstant, PregnancyId, DeliveryMrn, NewbornId)


#First get all cord gas measurements within 30 mins of birth.
cord_gas <- get_json_measurements(outcomes, delivery_times) %>% 
        inner_join(delivery_times, by=c("PregnancyId", "NewbornId", "DeliveryMrn")) %>%
        mutate(
            CollectionInstant=as_datetime(CollectionInstant),
            BirthInstant=as_datetime(BirthInstant),
            cord_gas_interval = as.numeric(CollectionInstant-BirthInstant, units='mins')) %>%
        filter(cord_gas_interval >= 0, cord_gas_interval <= 30)

print(head(cord_gas))
num_measurements <- cord_gas %>%
    group_by(PregnancyId, NewbornId, DeliveryMrn) %>%
    summarise(num_measurements=n())
print(table(num_measurements$num_measurements))
assert(1==0)

#---Everything below this point is pointless because the overwhelming majority of people only have 1 valid measurement (231 have both). 
#---Ah well.



# #Get all patients with both measurements.
# get_measurement_comparison <- function(cord_gas, type='Venous') {
#     venous_to_arterial_comparison <- cord_gas %>%
#         filter(!is.na(as.numeric(CordpH))) %>%

#         #Get the first measurements of each type.
#         group_by(PregnancyId, NewbornId, DeliveryMrn, Source) %>%
#         filter(cord_gas_interval == min(cord_gas_interval)) %>%

#         #Get patients with a venous and arterial measurement.
#         group_by(PregnancyId, NewbornId, DeliveryMrn) %>%
#         filter(any(Source==type) & any(Source=='Arterial')) %>%
#         select(PregnancyId, NewbornId, DeliveryMrn, Source, CordpH) %>%
#         ungroup() %>%

#         tidyr::pivot_wider(names_from='Source', values_from='CordpH') 
#     return(venous_to_arterial_comparison)
# }


# venous_to_arterial_comparison <- get_measurement_comparison(cord_gas)
# unspecified_to_arterial_comparison <- get_measurement_comparison(cord_gas, type='Cord (Unspecified)') %>%
#     rename(Unspecified = !!sym("Cord (Unspecified)"))

# print(nrow(venous_to_arterial_comparison))


# print(nrow(unspecified_to_arterial_comparison))
# assert(1==0)

# #-----TEST OUR ABILITY TO DEDUCE ARTERIAL pH FROM VENOUS-------

# out_dir <- "venous_to_arterial_plots"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# #####----- DATA PREP -----
# model_data <- venous_to_arterial_comparison %>%
#     ungroup() %>%
#     select(DeliveryMrn, Venous, Arterial)

# stopifnot(is.numeric(model_data$Venous), is.numeric(model_data$Arterial))

# # Fail loudly on any NA
# na_counts <- colSums(is.na(model_data))
# if (any(na_counts > 0)) {
#     stop("Unexpected missing values in model data:\n",
#          paste0("  ", names(na_counts)[na_counts > 0], ": ",
#                 na_counts[na_counts > 0], " NA(s)", collapse = "\n"),
#          call. = FALSE)
# }

# #####----- FOLDS (grouped by mother) -----
# K <- 5
# mrns <- unique(model_data$DeliveryMrn)
# mrn_fold <- setNames(sample(rep(1:K, length.out = length(mrns))), mrns)
# fold <- unname(mrn_fold[as.character(model_data$DeliveryMrn)])

# #####----- SPLINE HELPER -----
# # Natural cubic spline regression; df chosen by AIC on the data passed in
# fit_spline <- function(d, dfs = 3) {
#     fits <- lapply(dfs, function(df) lm(Arterial ~ ns(Venous, df = df), data = d))
#     fits[[which.min(vapply(fits, AIC, numeric(1)))]]
# }

# #####----- OUT-OF-FOLD PREDICTIONS -----
# y_all <- model_data$Arterial

# oof <- tibble(y = y_all, venous = model_data$Venous,
#               identity = model_data$Venous,   # baseline: arterial = venous
#               linear = NA_real_, spline = NA_real_)

# for (k in 1:K) {
#     tr <- fold != k; te <- fold == k

#     fit_lm <- lm(Arterial ~ Venous, data = model_data[tr, ])
#     oof$linear[te] <- predict(fit_lm, newdata = model_data[te, ])

#     # df chosen on the training fold only
#     fit_sp <- fit_spline(model_data[tr, ])
#     oof$spline[te] <- predict(fit_sp, newdata = model_data[te, ])
# }

# stopifnot(!anyNA(oof[c("identity", "linear", "spline")]))

# #####----- ACCURACY MEASURES -----
# eval_pred <- function(p, y) {
#     cal <- coef(lm(y ~ p))   # regress observed on predicted
#     tibble(
#         RMSE = sqrt(mean((y - p)^2)),
#         MAE  = mean(abs(y - p)),
#         R2   = 1 - sum((y - p)^2) / sum((y - mean(y))^2),
#         Mean_Error    = mean(p - y),     # bias; ideal = 0
#         Cal_Intercept = unname(cal[1]),  # ideal = 0
#         Cal_Slope     = unname(cal[2])   # ideal = 1
#     )
# }

# model_levels <- c("Baseline: arterial = venous", "Linear regression", "Natural spline")

# long <- oof %>%
#     pivot_longer(c(identity, linear, spline), names_to = "Model", values_to = "pred") %>%
#     mutate(Model = factor(recode(Model,
#         identity = "Baseline: arterial = venous",
#         linear   = "Linear regression",
#         spline   = "Natural spline"), levels = model_levels),
#         resid = y - pred)

# metrics <- long %>%
#     group_by(Model) %>%
#     reframe(eval_pred(pred, y))

# print(metrics, width = Inf)
# write.csv(metrics, file.path(out_dir, "metrics.csv"), row.names = FALSE)

# #####----- PLOTS -----
# save_plot <- function(p, name, w = 14, h = 5) {
#     print(p)
#     ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 300, bg = "white")
# }

# # 1. Predicted vs actual (out-of-fold)
# pred_vs_actual <- ggplot(long, aes(y, pred)) +
#     geom_point(alpha = 0.2, size = 0.8) +
#     geom_abline(linetype = "dashed", colour = "red") +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "blue") +
#     facet_wrap(~ Model) +
#     labs(x = "Actual arterial pH", y = "Predicted arterial pH (out-of-fold)",
#          title = "Predicted vs actual arterial pH",
#          subtitle = "Red = perfect prediction; blue = loess fit") +
#     theme_minimal()
# save_plot(pred_vs_actual, "predicted_vs_actual.pdf")

# # 2. Residual vs actual (out-of-fold)
# resid_vs_actual <- ggplot(long, aes(y, resid)) +
#     geom_point(alpha = 0.2, size = 0.8) +
#     geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "blue") +
#     facet_wrap(~ Model) +
#     labs(x = "Actual arterial pH", y = "Residual (actual - predicted)",
#          title = "Residuals vs actual arterial pH (out-of-fold)") +
#     theme_minimal()
# save_plot(resid_vs_actual, "residuals_vs_actual.pdf")

# # 3. Calibration plot: deciles of predicted value vs mean actual value
# n_bins <- 10
# cal_df <- long %>%
#     group_by(Model) %>%
#     mutate(bin = ntile(pred, n_bins)) %>%
#     group_by(Model, bin) %>%
#     summarise(mean_pred = mean(pred), mean_obs = mean(y), n = n(),
#               se = sd(y) / sqrt(n), .groups = "drop")

# cal_plot <- ggplot(cal_df, aes(mean_pred, mean_obs)) +
#     geom_abline(linetype = "dashed", colour = "red") +
#     geom_errorbar(aes(ymin = mean_obs - 1.96 * se, ymax = mean_obs + 1.96 * se),
#                   width = 0.002, alpha = 0.6) +
#     geom_line(colour = "blue") +
#     geom_point(aes(size = n)) +
#     facet_wrap(~ Model) +
#     labs(x = "Mean predicted arterial pH (decile of prediction)",
#          y = "Mean actual arterial pH",
#          title = "Calibration (out-of-fold, deciles of predicted pH)",
#          size = "Bin size") +
#     theme_minimal()
# save_plot(cal_plot, "calibration.pdf")

# #####----- FINAL FITS (all data) -----
# final_lm <- lm(Arterial ~ Venous, data = model_data)
# print(summary(final_lm))
# print(confint(final_lm))

# final_spline <- fit_spline(model_data)
# cat("Spline df chosen on full data:", length(coef(final_spline)) - 1, "\n")
# print(summary(final_spline))

# #------- TEST THRESHOLDS FOR ACIDEMIA ------
# out_dir <- "thresholds_test"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# acidemia_cutoff <- 7.1            # arterial pH <= this counts as acidemia
# targets         <- c(0.90, 0.95, 0.99)
# min_n           <- 20             # ignore thresholds with fewer babies than this below them
# step            <- 0.005

# # Inputs: each data frame has an Arterial column plus the measurement to be tested
# inputs <- list(
#     Venous      = list(measured = venous_to_arterial_comparison$Venous,
#                        arterial = venous_to_arterial_comparison$Arterial),
#     Unspecified = list(measured = unspecified_to_arterial_comparison$Unspecified,
#                        arterial = unspecified_to_arterial_comparison$Arterial)
# )

# # Fail loudly on NAs
# for (nm in names(inputs)) {
#     stopifnot(!anyNA(inputs[[nm]]$measured), !anyNA(inputs[[nm]]$arterial))
# }

# # Threshold grid covers the observed pH range of both measurement types
# all_ph     <- unlist(lapply(inputs, `[[`, "measured"))
# thresholds <- seq(floor(min(all_ph) / step) * step,
#                   ceiling(max(all_ph) / step) * step, by = step)

# wilson <- function(k, n, z = 1.96) {
#     p <- k / n
#     centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
#     half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
#     tibble(lower = centre - half, upper = centre + half)
# }

# # P(arterial acidemia | measured pH <= threshold), for each threshold
# threshold_curve <- function(measured, arterial, label) {
#     acidemic <- arterial <= acidemia_cutoff
#     n <- vapply(thresholds, function(t) sum(measured <= t), numeric(1))
#     k <- vapply(thresholds, function(t) sum(acidemic[measured <= t]), numeric(1))
#     ok <- n >= min_n

#     tibble(Type = label, Threshold = thresholds, n = n, k = k) %>%
#         filter(ok) %>%
#         mutate(Probability = k / n, wilson(k, n))
# }

# threshold_prob <- bind_rows(lapply(names(inputs), function(nm) {
#     threshold_curve(inputs[[nm]]$measured, inputs[[nm]]$arterial, nm)
# }))

# print(head(threshold_prob))
# write.csv(threshold_prob, file.path(out_dir, "acidemia_threshold_curves.csv"), row.names = FALSE)

# #------- REQUIRED THRESHOLDS FOR 90 / 95 / 99% -------
# # For each target: the highest threshold at which the probability is still >= target.
# # "point estimate" uses the observed fraction; "lower 95% bound" is more conservative
# # and uses the lower end of the Wilson interval.
# highest_meeting <- function(thr, value, target) {
#     ok <- value >= target
#     if (any(ok)) max(thr[ok]) else NA_real_
# }

# required <- threshold_prob %>%
#     group_by(Type) %>%
#     reframe(bind_rows(lapply(targets, function(tg) {
#         tibble(Target = tg,
#                `Point estimate`    = highest_meeting(Threshold, Probability, tg),
#                `Lower 95% bound`   = highest_meeting(Threshold, lower, tg))
#     })))

# print(required)
# write.csv(required, file.path(out_dir, "required_thresholds.csv"), row.names = FALSE)

# if (anyNA(required$`Point estimate`)) {
#     warning("Some targets were never reached for some measurement types (NA in table).")
# }

# #------- PLOT -------
# marker_df <- required %>%
#     select(Type, Target, Threshold = `Point estimate`) %>%
#     filter(!is.na(Threshold)) %>%
#     left_join(threshold_prob %>% select(Type, Threshold, Probability),
#               by = c("Type", "Threshold"))

# threshold_plot <- ggplot(threshold_prob, aes(Threshold, Probability, colour = Type, fill = Type)) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
#     geom_hline(yintercept = targets, linetype = "dashed", colour = "grey50") +
#     geom_line(linewidth = 1) +
#     geom_segment(data = marker_df,
#                  aes(x = Threshold, xend = Threshold, y = 0, yend = Probability),
#                  linetype = "dotted", show.legend = FALSE) +
#     geom_point(data = marker_df, size = 2.5, show.legend = FALSE) +
#     geom_text(data = marker_df, aes(label = sprintf("%.3f", Threshold), y = 0.02),
#               angle = 90, hjust = 0, vjust = -0.5, size = 3, show.legend = FALSE) +
#     scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, targets, 1),
#                        labels = scales::percent) +
#     labs(x = "pH threshold (measured pH <= threshold)",
#          y = sprintf("P(arterial pH <= %.1f)", acidemia_cutoff),
#          title = "Probability of arterial acidemia given a venous or unspecified pH below a threshold",
#          subtitle = sprintf("Dashed lines: %s. Markers: highest threshold still meeting each target. Ribbons: Wilson 95%% CI; thresholds with n < %d hidden.",
#                             paste0(targets * 100, "%", collapse = " / "), min_n)) +
#     theme_minimal() +
#     theme(legend.position = "bottom")

# print(threshold_plot)
# ggsave(file.path(out_dir, "acidemia_thresholds.png"), threshold_plot,
#        width = 9, height = 6, dpi = 300, bg = "white")