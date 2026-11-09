
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


#####-----CAN WE IMPUTE THE SOURCE OF A CORD PH MEASUREMENT FROM OTHER COVARIATES?-----
deliveries <- read.csv("preprocessed-deliveries-without-acidemia-outcomes.csv")


#####----- OUTPUT SETTINGS -----
out_dir <- "cordph_source_model_plots"   # change as needed
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#####----- DATA PREP -----
predictors <- c(covariates, "delivery_method")

model_data <- known_source %>%
    select(arterial, all_of(predictors)) %>%
    mutate(across(where(is.character), factor))

# Fail loudly if there are any missing values anywhere
na_counts <- colSums(is.na(model_data))
if (any(na_counts > 0)) {
    stop(
        "Unexpected missing values in model data:\n",
        paste0("  ", names(na_counts)[na_counts > 0], ": ",
               na_counts[na_counts > 0], " NA(s)", collapse = "\n"),
        "\nTotal rows with at least one NA: ", sum(!complete.cases(model_data)),
        call. = FALSE
    )
}

form  <- as.formula(paste("arterial ~", paste(predictors, collapse = " + ")))
X_all <- model.matrix(form, data = model_data)[, -1, drop = FALSE]
y_all <- model_data$arterial

# Belt and braces: model.matrix silently drops NA rows, so confirm alignment
stopifnot(nrow(X_all) == length(y_all), nrow(X_all) == nrow(model_data))

#####----- STRATIFIED K-FOLD -----
K <- 5
fold <- integer(nrow(model_data))
for (cls in 0:1) {
    idx <- which(y_all == cls)
    fold[idx] <- sample(rep(1:K, length.out = length(idx)))
}

#####----- OUT-OF-FOLD PREDICTIONS -----
xgb_params <- list(
    objective = "binary:logistic", eval_metric = "logloss",
    eta = 0.05, max_depth = 3, subsample = 0.8,
    colsample_bytree = 0.8, min_child_weight = 5
)

oof <- tibble(y = y_all, logit = NA_real_, xgb = NA_real_)

for (k in 1:K) {
    tr <- fold != k; te <- fold == k

    # Logistic regression
    fit_glm <- glm(form, data = model_data[tr, ], family = binomial)
    oof$logit[te] <- predict(fit_glm, newdata = model_data[te, ], type = "response")

    # XGBoost: choose nrounds by inner CV on the training fold only
    dtrain <- xgb.DMatrix(X_all[tr, , drop = FALSE], label = y_all[tr])
    cv <- xgb.cv(params = xgb_params, data = dtrain, nrounds = 1000, nfold = 5,
                 early_stopping_rounds = 25, verbose = 0, stratified = TRUE)
    best_n <- which.min(cv$evaluation_log$test_logloss_mean)

    fit_xgb <- xgb.train(params = xgb_params, data = dtrain, nrounds = best_n)
    oof$xgb[te] <- predict(fit_xgb, xgb.DMatrix(X_all[te, , drop = FALSE]))
}

# Every row should have received a prediction from both models
stopifnot(!anyNA(oof$logit), !anyNA(oof$xgb))

#####----- ACCURACY MEASURES -----
eval_model <- function(p, y, threshold = 0.5) {
    p_c   <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    pred  <- as.integer(p_c >= threshold)
    cal   <- coef(glm(y ~ qlogis(p_c), family = binomial))
    tibble(
        AUC         = as.numeric(pROC::auc(y, p_c, quiet = TRUE)),
        Accuracy    = mean(pred == y),
        Sensitivity = mean(pred[y == 1] == 1),
        Specificity = mean(pred[y == 0] == 0),
        Brier       = mean((p_c - y)^2),
        LogLoss     = -mean(y * log(p_c) + (1 - y) * log(1 - p_c)),
        Cal_Intercept = unname(cal[1]),   # ideal = 0
        Cal_Slope     = unname(cal[2])    # ideal = 1
    )
}

metrics <- bind_rows(
    Logistic = eval_model(oof$logit, oof$y),
    XGBoost  = eval_model(oof$xgb,   oof$y),
    .id = "Model"
) %>%
    add_row(Model = "Baseline (majority class)",
            Accuracy = max(mean(y_all), 1 - mean(y_all)))

print(metrics, width = Inf)

#####----- CALIBRATION PLOT -----
n_bins <- 10
cal_df <- oof %>%
    pivot_longer(c(logit, xgb), names_to = "Model", values_to = "p") %>%
    mutate(Model = recode(Model, logit = "Logistic", xgb = "XGBoost")) %>%
    group_by(Model) %>%
    mutate(bin = ntile(p, n_bins)) %>%
    group_by(Model, bin) %>%
    summarise(mean_pred = mean(p), obs = mean(y), n = n(),
              se = sqrt(obs * (1 - obs) / n), .groups = "drop")

calibration_plot <- ggplot(cal_df, aes(mean_pred, obs, colour = Model)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = pmax(obs - 1.96 * se, 0), ymax = pmin(obs + 1.96 * se, 1)),
                  width = 0.01, alpha = 0.5) +
    geom_line() + geom_point(aes(size = n)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Mean predicted P(arterial)", y = "Observed proportion arterial",
         title = "Calibration (out-of-fold, deciles of predicted risk)",
         size = "Bin size") +
    theme_minimal()

print(calibration_plot)
ggsave(file.path(out_dir, "calibration_plot.pdf"), calibration_plot,
       width = 7, height = 6, dpi = 300, bg = "white")

#####----- ROC CURVES -----
roc_list <- list(Logistic = roc(oof$y, oof$logit, quiet = TRUE),
                 XGBoost  = roc(oof$y, oof$xgb,   quiet = TRUE))

roc_plot <- ggroc(roc_list) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    labs(title = "ROC (out-of-fold)", colour = "Model") +
    theme_minimal()

print(roc_plot)
ggsave(file.path(out_dir, "roc_curves.pdf"), roc_plot,
       width = 7, height = 6, dpi = 300, bg = "white")

#####----- FINAL XGBOOST FIT FOR INTERPRETATION -----
dall <- xgb.DMatrix(X_all, label = y_all)
cv_all <- xgb.cv(params = xgb_params, data = dall, nrounds = 1000, nfold = 5,
                 early_stopping_rounds = 25, verbose = 0, stratified = TRUE)
final_xgb <- xgb.train(params = xgb_params, data = dall,
                       nrounds = which.min(cv_all$evaluation_log$test_logloss_mean))

importance <- xgb.importance(model = final_xgb)
print(importance)

importance_plot <- importance %>%
    dplyr::slice_max(Gain, n = 20) %>%
    ggplot(aes(x = reorder(Feature, Gain), y = Gain)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "Gain", title = "XGBoost feature importance (top 20)") +
    theme_minimal()

print(importance_plot)
ggsave(file.path(out_dir, "xgb_importance.pdf"), importance_plot,
       width = 7, height = 6, dpi = 300, bg = "white")


