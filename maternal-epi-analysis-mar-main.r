
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
setwd("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/epi-work/mar-res-main/")
pH_threshold <- 7.10
base_excess_threshold <- Inf #Functionally no upper limit to the base excess
#This contains the main analysis; changing these thresholds and working directory results in the sensitivity-analysis definition 



#####------SUPPLEMENTARY FUNCTIONS---------------


#Process MRNs.
process_maternal_mrns <- function(mrn) {
    
    #Remove 'R0' endings
    altered_mrn <- str_split(mrn, "-R")[[1]][1] 

    if(startsWith(altered_mrn, "0")) {
        altered_mrn <- str_remove(altered_mrn, "0")
    }
    return(altered_mrn)
}


#Define an asterisk function.
asterisk <- function(pval) {
    res <- case_when(
                pval < 0.001 ~ "***",
                pval < 0.01 ~ "**",
                pval < 0.05 ~ "*",
                .default=""
            )
    return(res)
 }


#Define a function to pull out gestation age in weeks.
gest_age <- function(ga) {
    split <- str_split_1(ga, "w")
    if (length(split)==2) {
        return(as.numeric(split[1]) + as.numeric(str_remove(split[2], "d"))/7)
    } else {return(NA)}
}

#Define a function to pull out rupture to delivery time in minutes.
rupture_to_delivery <- function(rup) {
    split <- str_split_1(rup, "h ")
    if (length(split)==2) {
        return(as.numeric(split[1])*60 + as.numeric(str_remove(split[2], "m")))
    } else {return(NA)}
}

#Extract the duration of each stage of labour.
labour_stage <- function(time) {
    days <- as.numeric(str_extract(time, "\\d+(?= day?)"))
    hours <- as.numeric(str_extract(time, "\\d+(?= hour?)"))
    mins <- as.numeric(str_extract(time, "\\d+(?= minute)"))
    total_duration_in_mins <- sum(24*60*days + 60*hours + mins, na.rm = TRUE)
    if(total_duration_in_mins==0) {
        total_duration_in_mins <- NA
    }
    return(total_duration_in_mins)
}

#Extract number of PREVIOUS pregnancies; a woman on her first pregnancy has a gravidity of 1
previous_pregnancies <- function(gp) {
    matches <- str_match(gp, "G(\\d+)P(\\d+)")
    g <- as.numeric(matches[2])
    return(g-1)
}

#Extract parity
parity <- function(gp) {
    matches <- str_match(gp, "G(\\d+)P(\\d+)")
    p <- as.numeric(matches[3])
    return(p)
}

#Does she have any previous C-sections indicated indirectly?
identify_previous_c_sections <- function(births, i) {
    mrn <- births$maternal_mrn[i] #Identify this mother
    baby_mrn <- births$baby_mrn[i] #Identify the baby
    mat_age <- births$raw_maternal_age[i] #Identify her age
    num_preg <- births$num_previous_pregnancies[i] #Identify number of previous pregnancies

    previous <- sum(
        births$maternal_mrn == mrn & #Previous births to the same mother
        !(births$baby_mrn == baby_mrn) & #Different baby
        births$raw_maternal_age <= mat_age & 
        births$num_previous_pregnancies < num_preg & (
            (births$delivery_method %in% c("C-Section (unplanned, with FTI)", "C-Section (planned)", "C-Section (unplanned, without FTI)")) | #Either THIS birth was by C-section, or
            (births$previous_c_section_indicated_directly == 1)) #C-section has already taken place by then
        )
    return(ifelse(previous>0, 1, 0))
}

#Define a function that implements Fisher's exact test if the contingency table is small enough, and the chisq elements otherwise.
measure_binary_differences <- function(x, y) {
    tryCatch({
        #Calculate the contingency table from vectors x and y.
        tab <- table(x, y)
        if(all(tab > 10)) { #This is sufficient to turn off the warnings.
            return(chisq.test(x, y)$p.value)

        } else {
            return(fisher.exact(x, y)$p.value)
        }}, error=function(e) {1}) #If these don't work, return 1- no s.f.
}


#Evaluate VIFs for a set of predictors.
evaluate_vifs <- function(df, outcomes, predictors) {

    #Pick an outcome.
    dummy_outcome <- outcomes[1] #It doesn't matter what this outcome is, it just needs to be real.

    #Define formula.
    form <- as.formula(paste0(dummy_outcome, " ~ ", paste(predictors, collapse=" + ")))

    #Fit a linear model. We're not going to use this, as VIFs only dependent on predictors; this is just for VIF calculation.
    lin_model_for_vifs <- lm(form, data=df)

    #Calculate VIFs.
    vifs <- car::vif(lin_model_for_vifs)


    return(max(vifs)) #This is the multicollinearity threshold we need to pass.
}


#Define a function to get odds ratios from a set of predictors.
evaluate_predictors <- function(df, outcomes, predictors) {
    overall_results <- data.frame()
    for (outcome in outcomes) {
        form <- as.formula(paste0(outcome, " ~ ", paste(predictors, collapse=" + ")))
        model <- glm(form, family=binomial(link="logit"), data=df)

        #Collate results.
        coef_summary <- as.data.frame(summary(model)$coefficients)
        ci <- as.data.frame(confint.default(model))

        results_df <- data.frame(
            variable=rownames(coef_summary),
            estimate=exp(coef_summary$Estimate),
            lower_conf=exp(ci[['2.5 %']]),
            upper_conf=exp(ci[['97.5 %']]),
            pval=as.numeric(coef_summary[["Pr(>|z|)"]]),
            outcome=outcome) %>% filter(!(variable=="(Intercept)"))
        overall_results <- rbind(overall_results, results_df)
    }

    #Perform FDR on the odds ratios.
    overall_results <- overall_results %>%
        mutate(pval=p.adjust(pval, method='fdr'), signif=asterisk(pval)) 
}


#Now evaluate a single model by checking VIFs and computing odds ratios.
evaluate_model <- function(df, outcomes, predictors) {

    #Check VIFs of the baseline variables (not the interaction terms.)
    max_vif <- evaluate_vifs(df, outcomes, predictors)
    assert(max_vif < 5) #If the VIFs are too high, this assertion will break the code.

    #Now
    overall_results <- evaluate_predictors(df, outcomes, predictors)

    return(overall_results) 

}

#Define a function which identifies 1:1 matches based on a certain set of variables.
try_exact_matching <- function(df, treatment_var, variable_prefixes, outcome_vars) {

    #Define result.
    odds_ratios <- data.frame(
        treatment_var=character(0),
        outcome_var=character(0),
        estimate=numeric(0),
        lower_conf=numeric(0),
        upper_conf=numeric(0),
        pval=numeric(0),
        num_matches=numeric(0),
        frac_matched=numeric(0)
        )
    index <- 1

    #Define relevant variables.
    variables_to_match <- c()
    for (prefix in variable_prefixes) {
        variables_to_match <- c(variables_to_match, colnames(df)[startsWith(colnames(df), prefix)])
    }


    #Define formula.
    form <- as.formula(paste0(treatment_var, " ~ ", paste(variables_to_match, collapse=" + ")))

    #Now match exactly, with replacement.
    tryCatch({
        match_object <- matchit(form, data=df, method="nearest", replace=TRUE,  exact=variables_to_match, distance=df$arbitrary_timestamp)

        matched <- as.data.table(match_data(match_object))

        #Run a weighted logistic regression model to identify the effect of the treatment variable. (McNemar's is hard using matchit, which doesn't want to give you a subclass if you're using replacement.)
        for (outcome_var in outcome_vars) {
            
            #Fit model.
            lm <- glm(as.formula(paste(outcome_var, "~", treatment_var)), family=binomial(link="logit"), data=matched, weights=matched$weights)    
            summary <- summary(lm)
            coefs <- as.data.frame(coef(summary))
            confs <- confint.default(lm)

            estimate <- exp(coefs[treatment_var, 'Estimate'])
            lower_conf <- exp(confs[treatment_var, '2.5 %'])
            upper_conf <- exp(confs[treatment_var, '97.5 %'])
            pval <- as.numeric(coefs[treatment_var, 'Pr(>|z|)'])

            #Save in table.
            log_result <- c(treatment_var, outcome_var, estimate, lower_conf, 
                    upper_conf, pval, sum(matched[[treatment_var]]==1), sum(matched[[treatment_var]]==1)/sum(df[[treatment_var]]==1))


            odds_ratios[index,] <- log_result
            index <- index + 1

            }
        }, error=function(e) {NULL})
    
    return(odds_ratios)
}



#Now put this together: given a set of outcomes, a set of strata, and a set of prefixes we'd ideally like to match, get ratios.
get_exact_match_odds_ratios <- function(df, outcome_vars, treatment_var, prefixes_to_match=c(), stratify_by=NULL) {

    #Define overall results.
    overall_results <- data.frame()


    #Go first by stratum (if strata are given)
    if (is.null(stratify_by)) {


        for (j in 1:length(prefixes_to_match)) {

            #Get odds ratios, matching for all ratios including this one.
            ors <- try_exact_matching(df, treatment_var, prefixes_to_match[1:j], outcome_vars) %>%
                mutate(subgroup=NULL, matched_up_to=prefixes_to_match[j], num_matched_vars=j) 
            #Add to overall results.
            overall_results <- rbind(overall_results, ors)
        }
    } else {
        #Define strata.
        strata <- c(colnames(df)[startsWith(colnames(df), stratify_by)])
        df[[paste0(stratify_by, "_ref")]] <- ifelse(rowSums(df[strata])==0, 1, 0)
        strata <- c(strata, paste0(stratify_by, "_ref"))

    
        for (stratum in strata) {

            subdf <- df %>% filter(.data[[stratum]]==1)
            #Then, for each prefix:
            for (j in 1:length(prefixes_to_match)) {
                
                #Get odds ratios, matching for all ratios including this one.
                ors <- try_exact_matching(subdf, treatment_var, prefixes_to_match[1:j], outcome_vars) %>%
                    mutate(subgroup=stratum, matched_up_to=prefixes_to_match[j], num_matched_vars=j) 

                #Add to overall results.
                overall_results <- rbind(overall_results, ors)

            } 
        }
    }

    return(overall_results)

}


#Define a function which identifies a matched cohort based on a certain set of variables (shortened version of above)
get_matched_cohort <- function(df, treatment_var, variable_prefixes) {


    #Define relevant variables.
    variables_to_match <- c()
    for (prefix in variable_prefixes) {
        variables_to_match <- c(variables_to_match, colnames(df)[startsWith(colnames(df), prefix)])
    }


    #Define formula.
    form <- as.formula(paste0(treatment_var, " ~ ", paste(variables_to_match, collapse=" + ")))

    #Now match exactly, with replacement.
    #tryCatch({
        match_object <- matchit(form, data=df, method="nearest", replace=TRUE, exact=variables_to_match, distance=df$arbitrary_timestamp)

        matched <- as.data.table(match_data(match_object))

        return(matched)

        # }, error=function(e) {
        #     return(NULL)
        # })
    

}



#Divide pH/acidemia into categories (We're just using the pH definition here)
divide_ph <- function(pH, base_excess) {
    res <- case_when(
        pH <= pH_threshold & base_excess < base_excess_threshold ~ "Acidemia",
        .default = "No acidemia"
    )
    return(res)
}

# # # #######----------PREPROCESS COVARIATES----------

# #Load raw data.
# births <- read.csv("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/combined_covariates-2025-11-03.csv") %>%
#     rename(baby_mrn=Baby.MRN, maternal_mrn=MRN) #Initially: 83953 unique births, 68110 unique maternal MRNs- roughly 1.23 births per mother


# #Define necessary covariates.
# births <- births %>% 
#     distinct(baby_mrn, .keep_all = TRUE) %>% #Get rid of multiple records per baby.- now 68110 unique maternal MRNs (83953 unique births)
#      mutate(

#         apgar_score = as.numeric(Apgar.5), #5-minute Apgar score--in final cohort, 65 have none recorded (0.75% of final cohort)
#         depressed_apgar = ifelse(!is.na(apgar_score) & (apgar_score < 7), 1, 0),

#          gestation_age_in_weeks=purrr::map_vec(GA, gest_age),

#          pregnancy_term=case_when(
#             gestation_age_in_weeks < 28 ~ "extremely preterm",
#             gestation_age_in_weeks < 32 ~ "very preterm",
#             gestation_age_in_weeks < 37 ~ "preterm",
#             gestation_age_in_weeks < 39 ~ "early term",
#             gestation_age_in_weeks < 41 ~ "full term",
#             gestation_age_in_weeks >= 41 ~ "late term", #Not enough (post-term (post 41 weeks) to be useful)
#             .default=NA     
#         ),

#         race=case_when(
#             Race %in% c('Black', 'White', 'Asian', 'American Indian', 'Pacific Island', 'East Indian') ~ Race,
#             Race == 'HLW-Hispanic Latino/White' ~ "White",
#             Race == 'HLB-Hispanic Latino/Black' ~ "Black",
#             Race == "" | is.na(Race) | is.null(Race) | grepl("Patient Declined", Race) | grepl("Unknown", Race) ~ "Unknown",
#             .default = "Other"
#         ),

#         #Record the presence of a previous cesarean
#         previous_c_section_indicated_directly=ifelse(
#             Delivery.Method=="VBAC" | Delivery.Method.1 == "VBAC" | (!is.na(as.numeric(Prior.C.S.Ct)) & (as.numeric(Prior.C.S.Ct)>0)) , 1, 0),

#         delivery_method = case_when(
#             Delivery.Method %in% c('VBAC, Spontaneous', 
#                 'Vaginal, Spontaneous', 
#                 'Vaginal, Forceps',
#                 'Vaginal, Vacuum (Extractor)', 
#                 'Vaginal Delivery', 
#                 'Vaginal, Breech', 
#                 'SVD*',
#                 'FAVD', 
#                 'Vaginal Delivery`', 
#                 'Vacuum Assisted Delivery',
#                 'Forceps Assisted Delivery') |
#             Delivery.Method.1 %in% c('VBAC, Spontaneous', 
#                 'Vaginal, Spontaneous', 
#                 'Vaginal, Forceps',
#                 'Vaginal, Vacuum (Extractor)', 
#                 'Vaginal Delivery', 
#                 'Vaginal, Breech', 
#                 'SVD*',
#                 'FAVD', 
#                 'Vaginal Delivery`', 
#                 'Vacuum Assisted Delivery',
#                 'Forceps Assisted Delivery') ~ "Vaginal",

#             Delivery.Method %in% c(
#                 'C-Section, Classical',
#                 'C-Section, Low Transverse',
#                 'C-Section, Unspecified',
#                 'C-Section, Low Vertical', 
#                 'Caesarean Section') |
#             Delivery.Method.1 %in% c(
#                 'C-Section, Classical',
#                 'C-Section, Low Transverse',
#                 'C-Section, Unspecified',
#                 'C-Section, Low Vertical', 
#                 'Caesarean Section') |
#             !(C.S.Categorization == "") |
#             !(C.S.Priority == "") |
#             !(Indications.for.Cesarean == "") |
#             !(C.Section.Indications== "") |
#             grepl("general", Anesthesia) |
#             grepl("spinal", Anesthesia)
#              ~ "C-Section",
#             .default = "Other/Unknown"
#             ),


#         #Temporal variables
#         Delivery.Date.Time=ymd_hms(Delivery.Date.Time),
#         original_delivery_date=Delivery.Date.Time,
#         arbitrary_timestamp=as.numeric(difftime(Delivery.Date.Time, ymd_hms("2000-01-01 00:00:00"), units="mins")), #Define each birth by an arbitrary timestamp 
#         year_of_delivery=year(Delivery.Date.Time),
#         month_of_delivery=month(Delivery.Date.Time),
#         season_of_delivery=case_when(
#             month_of_delivery > 2 & month_of_delivery <= 5 ~ "spring",
#             month_of_delivery > 5 & month_of_delivery <= 8 ~ "summer",
#             month_of_delivery > 8 & month_of_delivery <= 11 ~ "autumn",
#             month_of_delivery > 11 | month_of_delivery <= 2 ~ "winter",
#             .default = NA
#         ),
#         is_weekend=ifelse(wday(Delivery.Date.Time) %in% c(6, 7), 1, 0),
#         hour_of_delivery=hour(Delivery.Date.Time),
#         time_of_delivery = case_when(
#             hour_of_delivery <= 5 ~ "small hours",
#             hour_of_delivery <= 11 ~ "morning",
#             hour_of_delivery <= 17 ~ "afternoon",
#             hour_of_delivery <= 23 ~ "evening",
#             .default = NA
#         ),

#         #Length of stage 1 (before 10cm dilation) and stage 2 (after)
#         stage_1_length = purrr::map_vec(Stage.1.Length, labour_stage),
#         stage_2_length = purrr::map_vec(Stage.2.Length, labour_stage),

#         #Location of delivery
#         location=case_when(
#             Delivery.Location %in% c("CHESTER COUNTY HOSPITAL", "HUP", "MEDICAL CENTER OF PRINCETON", "PENNSYLVANIA HOSPITAL") ~ Delivery.Location,
#             .default="Other/Unknown"
#         ),


#         #If someone does get a C-section, note the presence of specific explanations:
#         fetal_tracing_indicated=ifelse(grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", Indications.for.Cesarean) | grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", C.Section.Indications), 1, 0),
#         other_unplanned=ifelse(grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", Indications.for.Cesarean) | grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", C.Section.Indications), 1, 0),

#         #Note stillbirths
#         stillbirth_recorded = ifelse(
#             Neonatal.Demise.=="Yes" |
#             Liv.Stat %in% c("ND", "FD"), 1, 0), #'living at delivery' seems to indicate something different.

#         #Now categorise maternal age.
#         maternal_age=case_when(
#             as.numeric(Age..Years.) < 20 ~ "19 and under",
#             as.numeric(Age..Years.) < 30 ~ "20-29",
#             as.numeric(Age..Years.) < 35 ~ "30-35",
#             as.numeric(Age..Years.) < 40 ~ "35-39",
#             as.numeric(Age..Years.) >= 40 ~ "40+",
#             .default=NA
#         ),

#         #Keep raw maternal age for later.
#         raw_maternal_age=as.numeric(Age..Years.),
        
#         #Extract number of previous pregnancies, and previous parity
#         num_previous_pregnancies=purrr::map_vec(GP, previous_pregnancies), 
#         num_parity=purrr::map_vec(GP, parity)

#         ) %>% #Still 68110 mothers.

#         #A check against multiple births: discard births where there's another in the last 30 days.
#         arrange(maternal_mrn, arbitrary_timestamp) %>% group_by(maternal_mrn) %>%
#         mutate(
#             num_births_within_1_month = purrr::map_dbl(row_number(), function(i) {

#                 current_time <- arbitrary_timestamp[i]
#                 current_baby <- baby_mrn[i]

#                 sum(!(baby_mrn==current_baby) &
#                      (abs(arbitrary_timestamp-current_time) < 24*60*30)) #Count any birth of a DIFFERENT baby within 30 days to the same mother.
#              })) %>% ungroup() %>% filter((year_of_delivery == 2017 & month_of_delivery >=4) | year_of_delivery > 2017, location=="HUP") 

# print(paste("Initially:", length(unique(births$baby_mrn)), "births at HUP after April 2017")) 


# births <- births %>%
#     filter(num_births_within_1_month==0) #%>% #Now 66813 records- around 2% of maternal MRNs lost, around 3% of births are multiple, fine- still 1.20 births per mother



# print(paste("After filtering out births to the same mother < 1 month apart:", length(unique(births$baby_mrn)), "births")) 


# births <- births %>%
#         select(baby_mrn, maternal_mrn, gestation_age_in_weeks, race, delivery_method, year_of_delivery,
#             season_of_delivery, month_of_delivery, is_weekend, hour_of_delivery, time_of_delivery, location, 
#              fetal_tracing_indicated, other_unplanned, stillbirth_recorded,
#              maternal_age,  num_previous_pregnancies, num_parity, apgar_score, depressed_apgar,
#              stage_1_length, stage_2_length, arbitrary_timestamp, original_delivery_date, raw_maternal_age, pregnancy_term,
#                 previous_c_section_indicated_directly, Induction, Augmentation, L.D.Complications, Rupture.Type, Forceps.Attempted, Vacuum.Attempted) #%>% 

# #print(paste("After filtering out births not at HUP:", length(unique(births$baby_mrn)), "births")) 

# births <- births %>%
#              mutate(delivery_method=case_when(
#                 delivery_method=="C-Section" & fetal_tracing_indicated ~ "C-Section (Nonreassuring EFM)",
#                 delivery_method=="C-Section" & other_unplanned ~ "C-Section (Labor Arrest)",
#                 delivery_method=="C-Section" ~ "C-Section (Planned)",
#                 delivery_method=="Vaginal" ~ "Vaginal",
#                 .default = "Other/Unknown") #Classify delivery methods.
#                 )


# #Save raw covariates.
# write.csv(births, "preprocessed-births-without-outcomes.csv")
# print("Covariates preprocessed")

# # #########-----------LINK ACIDEMIA OUTCOMES----------------------

# #Load preprocessed visits.
# births <- read.csv("preprocessed-births-without-outcomes.csv") 


# #Refine maternal MRNs. (We have 14270 unique maternal MRNs beforehand, and afterwards; we're not filing anything down.)
# births <- births %>% mutate(maternal_mrn=
#     purrr::map_vec(as.character(maternal_mrn), ~tryCatch(process_maternal_mrns(.x), error=function(e) {NA})))

# #Load offsets from individual files.
# offsets <- data.frame()
# for (k in 1:9) {
#     df <- read.csv(paste0("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/id_map/CTGData_offsets/0", k, "_offsetmap.csv"))
#     offsets <- rbind(offsets, df)
# }
# offsets <- offsets %>%
#         rename(filename=file)

# #Link trace-end timestamps to maternal MRNs, joining the ID-map to the trace file on filenames.
# id_map <- read.csv("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/id_map/Project Data trace_start_end_ctgdata_MASTER identified dataset 124776 tracings JM merge 12.5.22.csv") %>%
#     mutate(filename=purrr::map_vec(filename, function(x){str_split(x, "/")[[1]][2]})) %>%
#     rename(maternal_mrn=Mother.MRN) %>%
#     select(maternal_mrn, filename, trace_end) %>% #Keep the trace_end in order to number births consistently. Assume that each trace_end uniquely identifies a birth.
#     left_join(offsets, by="filename")

# #Now load lab results and link them to maternal MRNs.
# lab_results <- read.csv("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/CTGData_labs.csv") %>%
#     inner_join(id_map, by="filename") %>%
#     mutate(parsed_trace_end_date=parse_date_time(trace_end, "mdy HM"), parsed_order_date=ymd_hms(LabOrderDtime)) %>% #1.129 filenames per mother; 1.101 valid delivery dates per mother; 1.283 unique order timestamps per mother.
#     filter(LabTestName %in% c("pH Cord", "Base Excess Cord")) %>% #1.073 filenames per mother; 1.04 valid delivery dates per mother; 1.106 unique order timestamps per mother.

#     #If there are multiple labs of the same type ordered at the same timestamp, discard them as a proxy for multiple births.
#     group_by(maternal_mrn, parsed_order_date, LabTestName) %>% #1.073967 filenames/1.043493 delivery dates/1.106282 order timestamps
#     filter(n()==1) %>% 
#     ungroup() %>% 
    
#     #If there are multiple labs of the same type ordered at DIFFERENT timestamps, take the first one.
#     group_by(maternal_mrn, parsed_trace_end_date, LabTestName) %>%
#     filter(parsed_order_date==first(parsed_order_date)) %>% 
#     ungroup() %>% #1.072218 filenames, 1.043493 order dates, 1.068346 order timestamps

#     tidyr::pivot_wider(id_cols=c(maternal_mrn, parsed_trace_end_date, parsed_order_date, filename, offset_hours), names_from=LabTestName, values_from=LabResultValueFloat) %>%
#     clean_names() 


# #Now link results to births.
# births <- inner_join(births, lab_results, by="maternal_mrn", relationship="many-to-many") %>% #1.08944357774311 births per mother; 1.07106950944704 delivery dates per mother; 1.07115617958052 order dates per mother."
#     ungroup() %>% arrange(maternal_mrn) %>%
#     mutate(parsed_original_delivery_date=ymd_hms(original_delivery_date), offset_hours=ifelse(is.na(offset_hours), 0, offset_hours)) %>%
#     mutate(parsed_trace_end_date_offset = parsed_trace_end_date + hours(offset_hours),
#         parsed_original_delivery_date_offset = parsed_original_delivery_date + hours(offset_hours),
#         parsed_order_date_offset = parsed_order_date + hours(offset_hours)) %>%
#    mutate(mins_from_end_of_tracing_to_delivery=abs(as.numeric(difftime(parsed_original_delivery_date, parsed_trace_end_date, units="mins")))) 
   
# print(paste("After filtering out births without unique cord pH and cord base excess lab orders (where unique means one timestamp per lab type):", length(unique(births$baby_mrn)), "births")) 

# births <- births %>%
#    filter(mins_from_end_of_tracing_to_delivery <= 24*60*30) #%>% #Match within 1 month-- as the tracing often ends a few minutes or hours before the birth. #1.016 births per mother.


# births <- births %>%
#    group_by(baby_mrn) %>% filter(mins_from_end_of_tracing_to_delivery==min(mins_from_end_of_tracing_to_delivery)) %>% #Take the closest pair of tracing-end/delivery-date if there are multiple in a month
#    ungroup() %>%
#    distinct(maternal_mrn, baby_mrn, .keep_all = TRUE)

# print(paste("After filtering out births with lab results more than 1 month from delivery:", length(unique(births$baby_mrn)), "births")) 




# #Keep the first birth associated with every patient.
# births <- births %>% mutate(recorded=1)
# raw_births <- read.csv("preprocessed-births-without-outcomes.csv") %>%
#     left_join(select(births, c(baby_mrn, recorded)), by="baby_mrn") %>%
#     arrange(arbitrary_timestamp) %>%
#     group_by(maternal_mrn) %>% mutate(birth_number=row_number()) %>%
#     filter(recorded==1) %>% ungroup()

# births <- births %>%
#     inner_join(select(raw_births, baby_mrn, birth_number), by="baby_mrn") %>%
#     filter(birth_number==1)

# print(paste("After taking the first birth to each mother:", length(unique(births$baby_mrn)), "births")) 



# #Identify whether patients have records of previous C-sections.
# births$previous_c_section_indicated_indirectly <- purrr::map_vec(1:nrow(births), ~identify_previous_c_sections(births, .x), .progress=TRUE)
# births$previous_c_section_indicated <- ifelse(births$previous_c_section_indicated_directly + births$previous_c_section_indicated_indirectly > 0, 1, 0)

# #Categorise acidemia outcomes.
# births <- births %>%
#     mutate(p_h_cord = as.numeric(p_h_cord), base_excess_cord = as.numeric(base_excess_cord),
#         acidemia = ifelse((p_h_cord <= pH_threshold) & (base_excess_cord < base_excess_threshold), 1, 0),
#         acidemia_none = ifelse((p_h_cord <= pH_threshold) & (base_excess_cord < base_excess_threshold), 0, 1)) %>% #Main acidemia def is just pH, with base excess as sens. analysis
#     mutate(gestation_age_in_weeks=as.numeric(gestation_age_in_weeks)) %>%
#     filter(!is.na(gestation_age_in_weeks)) #%>% #Discards 30 births.


# print(paste("After filtering out births with no recorded gestational age:", length(unique(births$baby_mrn)), "births")) 


# births <- births %>%
#     filter(!is.na(p_h_cord)) %>% #Discards 16 births
#     filter(!is.na(base_excess_cord)) #%>% #Discards 3 births 

# print(paste("After filtering out births with either a non-existent pH or base excess results:", length(unique(births$baby_mrn)), "births")) 


# births <- births %>%
#     mutate(race=case_when(
#         race %in% c("Asian", "East Indian") ~ "Asian", #Groups together 75 East Indian births
#         race %in% c("American Indian", "Pacific Island") ~ "Other", #Groups together 18 births
#         .default = race
#     ), pregnancy_term = case_when(
#         pregnancy_term=="extremely preterm" ~ "very preterm", #only 70 births in this category
#         .default=pregnancy_term 
#     )) %>%
#     filter(!(delivery_method=="Other/Unknown")) #Removes 1 birth.

# print(paste("After filtering out births with unknown delivery method:", length(unique(births$baby_mrn)), "births")) 

# births <- births %>%
#     mutate(categorised_num_pregnancies = ifelse(is.na(num_previous_pregnancies), NA, ifelse(as.numeric(num_previous_pregnancies) > 1, "2+", as.character(num_previous_pregnancies))),
#             categorised_num_births = ifelse(is.na(num_parity), NA, ifelse(as.numeric(num_parity) > 0, 1, 0))) %>% #Categorises continuous variables.
#         filter(!is.na(categorised_num_births))

# print(paste("After filtering out births with unknown parity:", length(unique(births$baby_mrn)), "births")) 

# print(paste("After attaching lab results, we have ", nrow(births), "records.")) #8703 records

# #Save preprocessed births with outcomes.
# write.csv(births, "preprocessed-births-with-outcomes.csv")


# ######------DESCRIBE COHORT--------

#Load births back out of memory.
births <- read.csv("preprocessed-births-with-outcomes.csv")


#Assign complications to births
all_complications <- purrr::map_df(1:nrow(births), function(i) {
    comps <- str_remove(str_split_1(births$L.D.Complications[i], "\n"), "[.]")
    comps <- ifelse(comps=="Chorioamnionitis", "Chorioamionitis", ifelse(
        comps %in% c("None", ""), "None Recorded", comps)) #Correct obvious spelling error

    return(data.frame(baby_mrn=births$baby_mrn[i], complication=comps))
})

#Save complications.
write.csv(all_complications, "all_complications.csv")

#Pull out the top 5 complications
top_complications <- all_complications %>% 
    group_by(complication) %>% 
    summarise(NumBirths=n()) %>% 
    arrange(desc(NumBirths)) %>%
    rename(Characteristic=complication) %>%
    mutate(Attribute="Complications")

#Now summarise cohorts
summary_table <- top_complications[1:6,]
outcomes <- c(colnames(births)[endsWith(colnames(births), "acidemia")], "depressed_apgar")

#Define variables we are interested in summarising, and their names.
variables_to_summarise <- c(
    "delivery_method", "maternal_age", "pregnancy_term", "race", "categorised_num_pregnancies", "categorised_num_births",
    "previous_c_section_indicated",
      outcomes
)
variable_names <- c(
    "Delivery Method", "Maternal Age", "Pregnancy Term", "Race", "Previous Pregnancies", "Previous Births", 
    "Previous C-section?",
     paste("Outcome:", outcomes)
)

#Attach raw births.
for (i in 1:length(variables_to_summarise)) {
    var <- variables_to_summarise[i]
    sum <- births %>%
        group_by(!!sym(var)) %>%
        summarise(NumBirths=n()) %>%
        rename(Characteristic=!!sym(var)) %>%
        mutate(Attribute=variable_names[i])
    
    if (var=="pregnancy_term") {
        sum <- sum %>%
            mutate(Characteristic=factor(Characteristic, 
            levels=c("very preterm", "preterm", "early term", "full term", "late term"),
            labels=c("Very preterm (<32w)", "Preterm (32-37w)", "Early term (37-39w)", "Full term (39-41w)", "Late term (>41w)"))) %>%
            arrange(Characteristic)
    }

    if (startsWith(var, "pH")) {
        sum <- sum %>%
            filter(Characteristic==1)
        }

    summary_table <- rbind(summary_table, sum)
}

#Now attach cohort percentages.
total_births <- nrow(births)

summary_table <- summary_table %>% 
    mutate(NumBirths = paste0(NumBirths, " (", signif(NumBirths/total_births*100, 2), "%)"))

#Attach size of cohort to the top of the table
summary_table <- rbind(data.frame(Attribute="Overall", Characteristic="", NumBirths=total_births), summary_table)

#Print a summary
print(paste("In this cohort, we have", nrow(births), "births to", length(unique(births$maternal_mrn)), "unique parents."))

#Save the summary table.
write.csv(summary_table, "concise_summary_table.csv")

#Print interquartile range of age and gestational age.
print(paste("Maternal age, median:", median(births$raw_maternal_age)))
print(quantile(births$raw_maternal_age, c(0.25, 0.75)))
print(paste("Gestational age:", median(as.numeric(births$gestation_age_in_weeks))))
print(quantile(as.numeric(births$gestation_age_in_weeks), c(0.25, 0.75)))

#Plot Apgar scores.
p <- ggplot(births, aes(x=apgar_score)) +
    geom_histogram() + labs(
        x="5-Minute Apgar Score",
        y="Number of Births"
    ) +
theme_minimal(base_size = 18) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
    scale_fill_brewer(palette="Dark2") 
    ggsave("unadjusted_apgar.pdf", width=8, height=5, plot=p)



######-----CALCULATE UNADJUSTED DIFFERENCES-----------------

#Load preprocessed data.
births <- read.csv("preprocessed-births-with-outcomes.csv")


###PART A: UNADJUSTED DIFFERENCES IN OUTCOMES BY DELIVERY METHOD

#First, define outcomes of interest and their names.
outcomes <- c("acidemia", "depressed_apgar")
outcome_names <- c("Acidemia", "Apgar Score < 7")
delivery_method_exposures <- c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)")


#Define overall results.
unadjusted_differences_in_outcome_by_method <- data.frame()

#All differences in outcomes by method measured relative to vaginal delivery
for (dm in delivery_method_exposures) {
    #Select reference and exposure cohort.
    df <- births %>%
        filter(delivery_method %in% c(dm, "Vaginal")) 
    
    #Find significant differences.
    diffs <- df %>%
        summarise(across(all_of(outcomes), ~measure_binary_differences(delivery_method, .x))) %>%
        mutate(delivery_method=dm) %>%
        tidyr::pivot_longer(!delivery_method, names_to="outcome", values_to="pval")

    unadjusted_differences_in_outcome_by_method <- rbind(unadjusted_differences_in_outcome_by_method, diffs)
}

unadjusted_differences_in_outcome_by_method <- unadjusted_differences_in_outcome_by_method %>%
    mutate(
        signif=asterisk(pval)) 

#Get overall outcome rates by outcome.
rates <- births %>%
    group_by(delivery_method) %>%
    summarise(across(all_of(outcomes), mean)) %>%
    tidyr::pivot_longer(!delivery_method, names_to="outcome", values_to="rate") 

#Join this to the original dataframe.
unadjusted_differences_in_outcome_by_method <- unadjusted_differences_in_outcome_by_method %>%
    right_join(rates, by=c("delivery_method", "outcome")) %>%
        mutate(outcome=factor(outcome, levels=outcomes, labels=outcome_names),
        delivery_method=factor(delivery_method, levels=c(delivery_method_exposures, "Vaginal"), labels=c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal")))

#Save work.
write.csv(unadjusted_differences_in_outcome_by_method, "unadjusted_differences_in_outcome_by_method.csv")

#Plot this
p <- ggplot(unadjusted_differences_in_outcome_by_method, 
    aes(x=outcome, y=rate*100, fill=delivery_method)) +
     geom_col(position=position_dodge(width=0.85)) +
     labs(
        x = "Outcome",
        y = "Rate (% of Births)",
        fill = "Delivery Method",
        title = "Acidemia Rates by Delivery Method",
        subtitle = "Differences measured relative to vaginal delivery",
        caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    geom_text(aes(x=outcome, y=rate*100+1.5, label=signif(rate*100, 2)), position=position_dodge(width=0.85), size=4) +
    geom_text(aes(x=outcome, y=rate*100+3, label=signif), position=position_dodge(width=0.85), size=4) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
    ggsave("unadjusted_outcome_rates_by_method.pdf", width=15, height=6, plot=p)


###PART B: UNADJUSTED DIFFERENCES IN DELIVERY METHOD BY RACE

#Define overall results.
unadjusted_differences_in_method_by_race <- data.frame()
race_exposures <- c("Asian", "Black", "Other", "Unknown")

#All differences by race measured relative to white.
for (rc in race_exposures) {
    #Select reference and exposure cohort.
    df <- births %>%
        filter(race %in% c(rc, "White")) 
    
    #Find significant differences in each method of delivery.
    for (dm in c(delivery_method_exposures, "Vaginal")) {
        diffs <- df %>%
            summarise(pval= measure_binary_differences(race, delivery_method==dm)) %>%
            mutate(race=rc, delivery_method=dm) 
    
        unadjusted_differences_in_method_by_race <- rbind(unadjusted_differences_in_method_by_race, diffs)
    }
}

unadjusted_differences_in_method_by_race <- unadjusted_differences_in_method_by_race %>%
    mutate(
        signif=asterisk(pval)) 


#Now get rates of each outcome by race.
rates_by_race <- births %>%
    group_by(race) %>%
    mutate(TotalBirths=n()) %>%
    group_by(race, delivery_method) %>%
    summarise(rate=n()/first(TotalBirths)) 

#Join this to the original dataframe.
unadjusted_differences_in_method_by_race <- unadjusted_differences_in_method_by_race %>%
    right_join(rates_by_race, by=c("delivery_method", "race")) %>%
    mutate(delivery_method=factor(delivery_method, levels=c(delivery_method_exposures, "Vaginal"), labels=c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal")))

#Save work.
write.csv(unadjusted_differences_in_method_by_race, "unadjusted_differences_in_method_by_race.csv")

#Plot this
p <- ggplot(unadjusted_differences_in_method_by_race, 
    aes(fill=race, y=rate*100, x=delivery_method)) +
     geom_col(position=position_dodge(width=0.85)) +
     labs(
        x = "Delivery Method",
        y = "Rate (% of Births)",
        fill = "Race",
        title = "Delivery Method by Race",
        subtitle = "Differences in method measured relative to white race",
        caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    geom_text(aes(x=delivery_method, y=rate*100+3, label=signif(rate*100, 2)), position=position_dodge(width=0.85), size=5) +
    geom_text(aes(x=delivery_method, y=rate*100+6, label=signif), position=position_dodge(width=0.85), size=5) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
    scale_fill_brewer(palette="Dark2") 
    ggsave("unadjusted_method_by_race.pdf", width=15, height=6, plot=p)


# ###PART C: UNADJUSTED DIFFERENCES IN OUTCOME BY RACE

unadjusted_differences_in_outcome_by_race <- data.frame()

#All differences by race measured relative to white.
for (rc in race_exposures) {
    #Select reference and exposure cohort.
    df <- births %>%
        filter(race %in% c(rc, "White")) 
    
    #Find significant differences in each outcome.
    diffs <- df %>%
        summarise(across(all_of(outcomes), ~measure_binary_differences(race, .x))) %>%
        mutate(race=rc) %>%
        tidyr::pivot_longer(!race, names_to="outcome", values_to="pval")

    unadjusted_differences_in_outcome_by_race <- rbind(unadjusted_differences_in_outcome_by_race, diffs)
}

#No longer performing FDR correction as we're only interested in Black vs white comparisons.
unadjusted_differences_in_outcome_by_race <- unadjusted_differences_in_outcome_by_race %>% mutate(
        signif=asterisk(pval))


#Now get rates of each outcome by race.
rates_by_race <- births %>%
    group_by(race) %>%
    summarise(across(all_of(outcomes), mean)) %>%
    tidyr::pivot_longer(!race, names_to="outcome", values_to="rate")

#Join this to the original dataframe.
unadjusted_differences_in_outcome_by_race <- unadjusted_differences_in_outcome_by_race %>%
    right_join(rates_by_race, by=c("outcome", "race")) %>%
    mutate(outcome=factor(outcome, levels=outcomes, labels=outcome_names))

#Save work.
write.csv(unadjusted_differences_in_outcome_by_race, "unadjusted_differences_in_outcome_by_race.csv")

#Plot this
p <- ggplot(unadjusted_differences_in_outcome_by_race, 
    aes(fill=race, y=rate*100, x=outcome)) +
     geom_col(position=position_dodge(width=0.9)) +
     labs(
        x = "Outcome",
        y = "Rate (% of Births)",
        fill = "Race",
        title = "Acidemia Outcomes by Race",
        subtitle = "Differences in outcomes measured relative to white race"
    ) +
    geom_text(aes(x=outcome, y=rate*100+2, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=5) +
    geom_text(aes(x=outcome, y=rate*100, label=signif), position=position_dodge(width=0.9), size=10) +
    theme_minimal(base_size = 18) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
    scale_fill_brewer(palette="Dark2") 
    ggsave("unadjusted_outcome_by_race.pdf", width=15, height=6, plot=p)

##PART D: SUBGROUP ANALYSIS: ACIDEMIA OUTCOMES BY RACE, SUBGROUPED BY DELIVERY METHOD
#For sample size reasons, focus only on Black and white patients.

outcomes <- c("acidemia", "depressed_apgar")
outcome_names <- c("Acidemia", "Apgar Score < 7")
delivery_method_exposures <- c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)")

#Now measure racial differences, subgrouped by delivery_method.
#Select reference and exposure cohort.
df <- births %>%
    filter(race %in% c("Black", "White")) 
    
#Find significant differences in each outcome.
unadjusted_differences_in_outcome_by_race_and_method <- df %>%
    group_by(delivery_method) %>%
    summarise(across(all_of(outcomes), ~measure_binary_differences(race, .x))) %>%
    mutate(race="Black") %>%
    tidyr::pivot_longer(!c(race, delivery_method), names_to="outcome", values_to="pval")


#No longer performing FDR as we're interested in Black vs white comparisons, single outcome.
unadjusted_differences_in_outcome_by_race_and_method <- unadjusted_differences_in_outcome_by_race_and_method %>%
    mutate(
        signif=asterisk(pval)) 


#Now get rates of each outcome by race.
rates_by_race <- births %>%
    filter(race %in% c("Black", "White")) %>%
    group_by(race, delivery_method) %>%
    summarise(across(all_of(outcomes), mean)) %>%
    tidyr::pivot_longer(!c(race, delivery_method), names_to="outcome", values_to="rate")

#Join this to the original dataframe.
unadjusted_differences_in_outcome_by_race_and_method <- unadjusted_differences_in_outcome_by_race_and_method %>%
    right_join(rates_by_race, by=c("outcome", "race", "delivery_method")) %>%
    mutate(outcome=factor(outcome, levels=outcomes, labels=outcome_names),
    delivery_method=factor(delivery_method, 
        levels=c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal"), 
            labels=c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal"))) %>%
    group_by(outcome, delivery_method) %>%
    mutate(max_rate=max(rate)) %>% ungroup()

#Save work.
write.csv(unadjusted_differences_in_outcome_by_race_and_method, "unadjusted_differences_in_outcome_by_race_and_method.csv")

#Plot this
p <- ggplot(unadjusted_differences_in_outcome_by_race_and_method, 
    aes(fill=race, y=rate*100, x=delivery_method)) +
     geom_col(position=position_dodge(width=0.9)) +
     facet_wrap(~outcome, nrow=1) +
     labs(
        x = "Delivery Method",
        y = "Rate (% of Births)",
        fill = "Race",
        title = "Acidemia Outcomes by Race and Delivery Method",
        subtitle = "Differences in outcomes measured relative to white race"
    ) +
    geom_text(aes(x=delivery_method, y=rate*100+2, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=5) +
    geom_text(aes(x=delivery_method, y=max_rate*100+10, label=signif), size=5) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
    scale_fill_brewer(palette="Dark2") 
    ggsave("unadjusted_outcome_by_race_and_delivery_method.pdf", width=15, height=6, plot=p)


#####-----LOGISTIC REGRESSION MODELS OF DELIVERY METHOD AND OUTCOME-----------------

#Define potential predictors of delivery method or outcome (categorical and binary). [RIGHT NOW, leave out previous pregnancies as a predictor]
ctg_predictors <- c("race", "year_of_delivery", "season_of_delivery",
    "time_of_delivery", "maternal_age", "pregnancy_term") 
bin_predictors <- c("is_weekend", "previous_c_section_indicated", "categorised_num_births")
delivery_method_exposures <- c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)")

#Load preprocessed births.
births <- read.csv("preprocessed-births-with-outcomes.csv")

#Binarise categorical predictors with respect to a reference (except for delivery method, where we'll binarise all outcomes.)
ctg_refs <- c("White", 2017, "spring", "afternoon", "20-29", "full term")
vars_to_keep <- c("baby_mrn", bin_predictors, outcomes)


#Categorical predictors:
for (i in 1:length(ctg_refs)) {
    var <- ctg_predictors[i]
    ref <- ctg_refs[i]
    vals <- unique(births[[var]]) 
    for (val in vals) {
        if (!(val==ref)) {
            new_var <- paste0(var, "_", val)
            births[[new_var]] <- ifelse(births[[var]]==val, 1, 0)
            vars_to_keep <- c(vars_to_keep, new_var)
        }
    }
}

#Delivery method:
for (dm in c(delivery_method_exposures, "Vaginal")) {
    new_var <- paste0("delivery_method_", dm)
    births[[new_var]] <- ifelse(births$delivery_method==dm, 1, 0)
    vars_to_keep <- c(vars_to_keep, new_var)
}

#Keep only important columns and save.
df <- births %>%
    select(all_of(vars_to_keep)) %>%
    clean_names()
write.csv(df, "binarised-births.csv")

#Load binarised data.
births <- read.csv("binarised-births.csv")

#Now define the feasible predictors of delivery method.
method_predictors <- c(bin_predictors)
for (pred in ctg_predictors) {
    method_predictors <- c(method_predictors, colnames(births)[startsWith(colnames(births), pred)])
}

#The feasible predictors of outcome are just the above, plus method itself.
delivery_methods <- colnames(births)[startsWith(colnames(births), "delivery_method_")]


# ###PART A: PREDICTING DELIVERY METHOD

method_predictors_in_order <- c("maternal_age_19_and_under", "maternal_age_30_35",
    "maternal_age_35_39", "maternal_age_40", "pregnancy_term_very_preterm", "pregnancy_term_preterm",
    "pregnancy_term_early_term", "pregnancy_term_late_term",
    "race_asian", "race_black", 
    "race_other", "race_unknown", 
    "categorised_num_births", "previous_c_section_indicated",
    "year_of_delivery_2018", "year_of_delivery_2019", "year_of_delivery_2020",
    "season_of_delivery_summer", "season_of_delivery_autumn", "season_of_delivery_winter",
    "is_weekend", "time_of_delivery_small_hours", "time_of_delivery_morning", "time_of_delivery_evening"
    )
method_predictors_in_order_names <- c("Maternal Age: <20", "Maternal Age: 30-35", "Maternal Age: 35-39", "Maternal Age: 40+",
    "Gest. Age: <32w", "Gest. Age: 32-37w", "Gest. Age: 37-39w", "Gest. Age: 41w+",
    "Race: Asian", "Race: Black", "Race: Other", "Race: Unknown",
    "Previous Births", "Previous CD",
    "Y: 2018", "Y: 2019", "Y: 2020", "S: Jun-Aug", "S: Sep-Nov", "S: Dec-Feb",
    "Weekend", "H: 12a-6a", "H: 6a-12p", "H: 6p-12a") 


#Plot the model.
plot_lr_model <- function(overall_results, predictors, predictor_names, 
    outcomes, outcome_names, savename, title, nrow=1, width=10, height=10, base_size=15, 
    errorbar_height=0.2, asterisk_offset=0.5, x_axis_size=10, subtitle="") {
    
    #Label the predictors.
    overall_results <- overall_results %>%
        mutate(variable=factor(variable, levels=rev(predictors), labels=rev(predictor_names)),
            outcome=factor(outcome, levels=outcomes, labels=outcome_names)) %>%
        group_by(outcome) %>% mutate(annotation_position=max(upper_conf)*(1+asterisk_offset)) %>% ungroup() %>%
        filter(lower_conf > 0 & is.finite(upper_conf))
    


    #Now make facet plots.
    p <- ggplot(overall_results, aes(x=estimate, y=variable, color=ifelse(pval < 0.05, "red", "black"))) +
        geom_vline(xintercept=1, linetype="dashed") +
        geom_point(show.legend = FALSE) + geom_errorbarh(aes(xmin=lower_conf, xmax=upper_conf, 
            color=ifelse(pval < 0.05, "red", "black")), height=errorbar_height, show.legend = FALSE) +
        labs(
            x="Odds Ratio", 
            y="Risk Factor", 
            #title=title,
            #subtitle=subtitle,
            #caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
            ) +
        scale_x_log10() +
        coord_cartesian(clip = "off") +
        facet_wrap(~outcome, nrow=nrow, scales = "free_x") + 
        geom_text(aes(x=annotation_position, y=variable, label=signif, color=ifelse(pval<0.05, "red", "black")), show.legend = FALSE) +
        theme_minimal(base_size = base_size) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            plot.tag = element_text(angle=45),
            axis.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, size=x_axis_size),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()) +
        scale_color_identity() + guides(color="none")
        ggsave(paste0(savename, ".pdf"), width=width, height=height, plot=p)

}

#Assign names to delivery method for this graph.
delivery_method_names <- c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal")

#Calculate odds ratios relevant to DELIVERY METHOD, save, and plot.
dm_results <- evaluate_model(births, delivery_methods, method_predictors)
write.csv(dm_results, "delivery-method-lr.csv")
plot_lr_model(dm_results, method_predictors_in_order, method_predictors_in_order_names,
    delivery_methods, delivery_method_names,  "lr-delivery-method", "Risk Factors for Delivery Method", width=12)


###PART B: PREDICTING ACIDEMIA OUTCOME
#Right now, we're not testing for interactions, but we could!
outcomes <- c("acidemia", "depressed_apgar")
outcome_names <- c("Acidemia", "Apgar Score < 7")


outcome_predictors <- c(method_predictors_in_order, delivery_methods[1:3])
outcome_predictor_names <- c(method_predictors_in_order_names, delivery_method_names[1:3]) #Vaginal birth is the default

#Calculate odds ratios.
outcome_results <- evaluate_model(births, outcomes, outcome_predictors) #LR models use 
write.csv(outcome_results, "acidemia-outcome-lr.csv")

#Plot these.
plot_lr_model(outcome_results, outcome_predictors, outcome_predictor_names,
    outcomes, outcome_names, "lr-acidemia-outcome", "Risk Factors for Poor Fetal Health",
    nrow=1, width=9, base_size=15, asterisk_offset = 1)



####---------------COARSENED EXACT MATCHING FOR ACIDEMIA OUTCOMES------

#Match Black and white patients according to: maternal age, pregnancy term, mode of delivery, previous birth, whether or not they have previous C-sections; if possible, year of delivery and time of day.
#We use nearest-neighbour matching where each treated (Black) patient is matched to a white patient.

#First, identify exposure and reference cohort.
raw_births <- read.csv("preprocessed-births-with-outcomes.csv")
births <- read.csv("binarised-births.csv")

#Identify delivery methods.
delivery_methods <- colnames(births)[startsWith(colnames(births), "delivery_method_")]
delivery_method_names <- c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "C-Section (Planned)", "Vaginal")


#Arrange births in order they occurred, so nearest-neighbour matching will match births which occurred close together in time.
timestamps <- read.csv("preprocessed-births-without-outcomes.csv") %>%
    select(baby_mrn, arbitrary_timestamp)
births <- births %>% inner_join(timestamps, by="baby_mrn") %>%
    arrange(arbitrary_timestamp)

#Define outcomes.
outcomes <- c("acidemia", "depressed_apgar")
outcome_names <- c("Acidemia", "Apgar Score < 7")

#Define exposed and treated cohort.
exposed_births <- raw_births %>%
    filter(race=="Black") %>%
    select(baby_mrn) %>% inner_join(births, by="baby_mrn")
ref_births <- raw_births %>%
    filter(race=="White") %>%
    select(baby_mrn) %>% inner_join(births, by="baby_mrn")

#Define the maximum number of matches (i.e. number in the treatment group.)
max_matches <- nrow(exposed_births)

#Combine reference and exposed births.
df <- rbind(exposed_births, ref_births) %>%
    arrange(arbitrary_timestamp)

#Save this cohort.
write.csv(df, 'narrow-cohort.csv')

#Define the ideal variables to match, in order of importance.
prefixes_to_match <- c("maternal_age", "pregnancy_term", "previous_c_section_indicated",
    "categorised_num_births")


###PART A: UNSTRATIFIED DIFFERENCES IN DELIVERY METHOD BY RACE, WITH EXACT MATCHING
#Here we have no strata- we're just interested in matching Black to White patients as closely as possible.
dm_results_exact_matched <- get_exact_match_odds_ratios(df, delivery_methods, "race_black", prefixes_to_match=prefixes_to_match) %>%
    mutate(signif=asterisk(as.numeric(pval)))
write.csv(dm_results_exact_matched, "exact-matching-delivery-method-by-race.csv")

###PART B: DIFFERENCES IN DELIVERY METHOD BY RACE, WITH EXACT MATCHING- STRATIFIED BY INDIVIDUAL VARIABLES
stratified_dm_results <- data.frame()
for (var in prefixes_to_match) {
    dm_results_exact_matched <- get_exact_match_odds_ratios(df, delivery_methods, "race_black", stratify_by=var, prefixes_to_match=prefixes_to_match[!(prefixes_to_match==var)]) %>%
        mutate(signif=asterisk(as.numeric(pval)))
    stratified_dm_results <- rbind(stratified_dm_results, dm_results_exact_matched)
}
write.csv(stratified_dm_results, "exact-matching-delivery-method-by-race-stratified.csv")

#Plot overall results.
plot_matched_results <- function(results, outcome_name, title, filename,  with_stratification=FALSE, width=12, height=8, cex=1.2, xpos=c(0.2, 1.2)) {


    #Open PDF
    pdf(file = paste0("forest-plots/", filename, ".pdf"), 
            width = width, 
            height = height, 
            family = "Helvetica")

        if (with_stratification) {

            #Get only the full-matched results.
            res <- results %>% group_by(outcome_var, subgroup) %>% 
            mutate(estimate=as.numeric(estimate), lower_conf=as.numeric(lower_conf), upper_conf=as.numeric(upper_conf), 
                frac_matched=as.numeric(frac_matched), pval=as.numeric(pval)) %>%
            filter(num_matched_vars==max(num_matched_vars), is.finite(log(lower_conf)), is.finite(log(upper_conf))) %>%
            ungroup() %>%
                arrange(outcome_var, subgroup) %>% group_by(outcome_var) %>%
                mutate(outcome_label=ifelse(row_number()==ceiling(n()/2), as.character(outcome_var), "")) %>%
                ungroup() %>%
            mutate(signif=asterisk(pval)) #pval=p.adjust(pval, method='fdr'),
            
            

            #Work out where to draw lines
            boundaries <- res %>%
                group_by(outcome_var) %>% summarise(size=n()) %>%
                ungroup() %>% mutate(lines=sum(size)-cumsum(size)+0.5) 


            forest(x=log(res$estimate), 
                ci.lb=log(res$lower_conf),
                ci.ub=log(res$upper_conf),
                slab = res$outcome_label,
                header = paste("Outcome", outcome_name),
                psize=1,
                xlab = "Odds Ratio: Black relative to white",
                ilab = data.frame(stratum=res$subgroup, matches=paste0(res$num_matches, " (", signif(res$frac_matched*100, 2), ")"), pval=res$signif),
                ilab.lab = c("Subgroup", "Matches (% of Total)", ""),
                ilab.pos = c(4, 4, 2),
                ilab.xpos= log(xpos),
                col=ifelse(res$pval < 0.05, "#FF0000", "#000000"),
                col.lab = "#000000",  
                annotate=TRUE,
                refline = 0,
                atransf = function(x) exp(x),
                cex=cex,
                #main=title,
                cex.main=cex*1.1)
        abline(h=boundaries$lines, lty=2,col='gray70')
        }

     else {

            #Get only the full-matched results.
            res <- results %>% group_by(outcome_var) %>% 
            filter(num_matched_vars==max(num_matched_vars)) #%>%
            #mutate(pval=p.adjust(pval, method='fdr'))

            forest(x=log(res$estimate), 
            ci.lb=log(res$lower_conf),
            ci.ub=log(res$upper_conf),
            slab = res$outcome_var,
            header = paste("Outcome:", outcome_name),
            psize=1,
            xlab = "Odds Ratio: Black relative to white",
            ilab = data.frame(pval=res$signif),
            ilab.lab = c(""),
            ilab.pos = c(2, 2, 2),
            ilab.xpos= log(xpos),
            col=ifelse(res$pval < 0.05, "#FF0000", "#000000"),
            annotate=TRUE,
            refline = 0,
            efac=0.5,
            atransf = function(x) exp(x),
            cex=cex)
            #main=title)
    }


    dev.off()

}

#Plot overall, full-matched results stratified by individual variables.
unstratified_results <- read.csv("exact-matching-delivery-method-by-race.csv") %>%
    mutate(outcome_var=factor(outcome_var, levels=delivery_methods, labels=delivery_method_names))
plot_matched_results(unstratified_results, "Delivery Method", "Overall effect of Black race (relative to white) on mode of delivery",
     "delivery-method-race-unstratified", width=10, height=6, cex=1.4,  xpos=c(4)) 


#Define subgroups.

subgroups_in_order <- c("maternal_age_19_and_under", "maternal_age_ref", "maternal_age_30_35",
    "maternal_age_35_39", "maternal_age_40", "pregnancy_term_very_preterm", "pregnancy_term_preterm",
    "pregnancy_term_early_term", "pregnancy_term_ref", "pregnancy_term_late_term",
    "categorised_num_births", "categorised_num_births_ref", "previous_c_section_indicated", "previous_c_section_indicated_ref")

subgroups_in_order_names <- c("Maternal Age: <20", "Maternal Age: 20-29", "Maternal Age: 30-35", "Maternal Age: 35-39", "Maternal Age: 40+",
    "Gest. Age: <32w", "Gest. Age: 32-37w", "Gest. Age: 37-39w", "Gest. Age: 39-41w", "Gest. Age: 41w+",
    "Previous Births: Yes", "Previous Births: No", "Previous C-section: Yes", "Previous C-section: No") 

# Plot overall, full-matched results stratified by individual variables.
stratified_results <- read.csv("exact-matching-delivery-method-by-race-stratified.csv") %>%
    mutate(outcome_var=factor(outcome_var, levels=delivery_methods, labels=delivery_method_names)) %>%
    filter(!(startsWith(subgroup, "year") | startsWith(subgroup, "time"))) %>% #We're not interested in these subgroups
    mutate(subgroup=factor(subgroup, levels=subgroups_in_order, labels=subgroups_in_order_names))

plot_matched_results(stratified_results, "", "Subgroup-specific effect of Black race (relative to white) on mode of delivery",
     "delivery-method-race-stratified", with_stratification = TRUE, width=12, height=12, xpos=c(0.0007, 0.02, 54)) 


###PART C: DIFFERENCES IN OUTCOMES, WITH EXACT MATCHING (INC ON DELIVERY METHOD)
dm_results_exact_matched <- get_exact_match_odds_ratios(df, outcomes, "race_black", prefixes_to_match=c("delivery_method", prefixes_to_match)) %>%
    mutate(signif=asterisk(as.numeric(pval)))
write.csv(dm_results_exact_matched, "exact-matching-outcomes-by-race.csv") 

unstratified_results <- read.csv("exact-matching-outcomes-by-race.csv") %>%
    mutate(outcome_var=factor(outcome_var, levels=outcomes, labels=outcome_names))
plot_matched_results(unstratified_results, "Acidemia", "Overall effect of Black race (relative to white) on acidemia",
     "acidemia-race-unstratified", width=10, height=6, cex=1.4, xpos=c(3)) 

# ###PART D: DIFFERENCES IN OUTCOMES, WITH EXACT MATCHING (INC ON DELIVERY METHOD)- STRATIFIED BY VARIABLE

stratified_dm_results <- data.frame()
prefixes_with_delivery <- c(prefixes_to_match, "delivery_method")
for (var in prefixes_with_delivery) {
    dm_results_exact_matched <- get_exact_match_odds_ratios(df, outcomes, "race_black", stratify_by=var, prefixes_to_match=prefixes_with_delivery[!(prefixes_with_delivery==var)]) %>%
        mutate(signif=asterisk(as.numeric(pval)))
    stratified_dm_results <- rbind(stratified_dm_results, dm_results_exact_matched)
}
write.csv(stratified_dm_results, "exact-matching-outcome-by-race-stratified.csv")


#Plot overall, full-matched results stratified by DELIVERY METHOD ONLY.
stratified_results <- read.csv("exact-matching-outcome-by-race-stratified.csv") %>%
    filter(startsWith(outcome_var, "acidemia")) %>%
    mutate(outcome_var=factor(outcome_var, levels=outcomes, labels=outcome_names)) %>%
    filter(startsWith(subgroup, "delivery_method")) %>% #We're not interested in these subgroups
    mutate(subgroup=factor(subgroup, levels=c(delivery_methods, subgroups_in_order), labels=c(delivery_method_names, subgroups_in_order_names)))

plot_matched_results(stratified_results, "", "Delivery-mode-specific effect of Black race (relative to white) on acidemia",
     "acidemia-race-stratified-dm-only", with_stratification = TRUE, width=12, height=6, xpos=c(0.0015, 0.07, 20), cex=1.3) 

#Plot overall, full-matched results stratified by DELIVERY METHOD ONLY (with apgar score).
stratified_results <- read.csv("exact-matching-outcome-by-race-stratified.csv") %>%
    mutate(outcome_var=factor(outcome_var, levels=outcomes, labels=outcome_names)) %>%
    filter(startsWith(subgroup, "delivery_method")) %>% #We're not interested in these subgroups
    mutate(subgroup=factor(subgroup, levels=c(delivery_methods, subgroups_in_order), labels=c(delivery_method_names, subgroups_in_order_names)))

plot_matched_results(stratified_results, "", "Delivery-mode-specific effect of Black race (relative to white) on birth outcomes",
     "acidemia-race-stratified-dm-only-with-apgar", with_stratification = TRUE, width=10, height=6, xpos=c(0.004, 0.06, 20)) 



####PART E: DELIVERY METHOD, STRATIFIED BY ACIDEMIA
stratified_results <- get_exact_match_odds_ratios(df, delivery_methods, "race_black", stratify_by="acidemia", prefixes_to_match=prefixes_to_match) %>%
        mutate(signif=asterisk(as.numeric(pval))) %>%
        mutate(outcome_var=factor(outcome_var, levels=delivery_methods, labels=delivery_method_names)) %>%
        mutate(subgroup=factor(subgroup, levels=c("acidemia_ref", "acidemia"), labels=c("No acidemia", "Acidemia")))

write.csv(stratified_results, "delivery-method-stratified-outcome.csv")

plot_matched_results(stratified_results, "", "Acidemia-specific effect of Black race (relative to white) on delivery method",
     "dm-race-stratified-outcome-only", with_stratification = TRUE, width=12, height=12, xpos=c(0.05, 0.1, 50)) 


######-----------LOOK FOR DIFFERENCES IN PH------------------------

#Load narrow cohort (Black and White patients only.)
births <- read.csv("preprocessed-births-with-outcomes.csv")
df <- read.csv("narrow-cohort.csv") %>%
    inner_join(select(births, c(baby_mrn, p_h_cord, base_excess_cord)), by="baby_mrn")

print(colnames(df))

#Narrow this down to race, delivery method, and continuous results.
df <- df %>%
    mutate(Race=ifelse(race_black==1, "Black", "White"),
        DeliveryMethod=case_when(
            delivery_method_c_section_nonreassuring_efm == 1 ~ "C-Section (Nonreassuring EFM)", 
            delivery_method_c_section_labor_arrest == 1 ~ "C-Section (Labor Arrest)",
            delivery_method_c_section_planned == 1 ~ "C-Section (Planned)",
            delivery_method_vaginal == 1 ~ "Vaginal"
        )) %>% rename(pH=p_h_cord, BaseExcess=base_excess_cord)

#Match Black to white patients.
matched_df <- get_matched_cohort(df, "race_black", c("maternal_age", "pregnancy_term", "previous_c_section_indicated",
    "categorised_num_births", "delivery_method"))

#Get overall mean pH by race; test for differences in overall mean.
overall_mean_pH <-  df %>% group_by(Race) %>% summarise(mean_pH=mean(pH))
overall_mean_pH_diff = glm(pH ~ race_black, family="gaussian", data=df)
print("Overall differences in pH BEFORE MATCHING")
print(overall_mean_pH)
print(summary(overall_mean_pH_diff))

overall_mean_pH <- matched_df %>% group_by(Race) %>% summarise(mean_pH=weighted.mean(pH, w=weights))
overall_mean_pH_diff = glm(pH ~ race_black, family="gaussian", data=matched_df, weight=weights)
print("Overall differences in pH AFTER MATCHING")
print(overall_mean_pH)
print(summary(overall_mean_pH_diff))



#Categorise pH in the overall matched database (this also includes base excess now.)
matched_df$categorised_pH <- purrr::map2_chr(
  matched_df$pH,
  matched_df$BaseExcess,
  divide_ph,
  .progress = TRUE
)


#######--------PART A: PH DIFFERENCES BY RACE. Get the rates of pH-category in each race.---------
matched_rates <- matched_df %>%
    group_by(Race) %>% mutate(OverallNumUnits=sum(weights)) %>%
    group_by(Race, categorised_pH) %>% summarise(rate=sum(weights)/first(OverallNumUnits)) %>%
    ungroup()

#Test for SIGNIFICANT differences in pH-category in each race.
diffs_in_ph <- data.frame(categorised_pH=character(0), estimate=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0), pval=numeric(0))
index <- 1
for (categ in unique(matched_df$categorised_pH)) {
    copy_matched_df <- data.frame(matched_df)
    copy_matched_df$target <- ifelse(copy_matched_df$categorised_pH==categ, 1, 0)
    lm <- glm(target ~ race_black, data=copy_matched_df, family=binomial(link="logit"), weights=weights)
    summary <- summary(lm)
    coefs <- as.data.frame(coef(summary))
    confs <- confint.default(lm)

    estimate <- exp(coefs["race_black", 'Estimate'])
    lower_conf <- exp(confs["race_black", '2.5 %'])
    upper_conf <- exp(confs["race_black", '97.5 %'])
    pval <- as.numeric(coefs["race_black", 'Pr(>|z|)'])
    diffs_in_ph[index,] <- c(categ, estimate, lower_conf, upper_conf, pval)
    index <- index + 1
}

#Link rates and significance together.
matched_rates <- matched_rates %>%
    inner_join(diffs_in_ph, by="categorised_pH") %>%
    mutate(signif=asterisk(pval)) %>%
    mutate(categorised_pH=factor(categorised_pH, levels = c("Acidemia", "No acidemia"))) %>%
    group_by(categorised_pH) %>% mutate(max_rate=max(rate)) %>% ungroup()

#Plot categorised pH by race in the overall matched cohort.
p <- ggplot(matched_rates, aes(x=categorised_pH, y=rate*100, fill=Race)) +
    geom_col(position=position_dodge(width=0.9)) +
     labs(
        x = "pH",
        y = "Rate (% of Births)",
        fill = "Race",
        title = "pH by Race in Overall Matched Cohort",
        #subtitle = "Differences measured relative to white race",
        #caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    geom_text(aes(x=categorised_pH, y=rate*100+3, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=4) +
    geom_text(aes(x=categorised_pH, y=max_rate*100+5, label=signif), size=4) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
    ggsave("pH_diffs_overall.pdf", width=6, height=4, plot=p)


#######--------PART B: PH DIFFERENCES BY DELIVERY METHOD---------=--
stratified_matched_rates <- data.frame()

for (dm in unique(df$DeliveryMethod)) {
    subdf <- df %>% filter(DeliveryMethod==dm)

    #Match cohort.
    matched_subdf <- get_matched_cohort(subdf, "race_black", c("maternal_age", "pregnancy_term", "previous_c_section_indicated",
    "categorised_num_births"))

    #Look for differences in mean pH.
    lm <- glm(pH ~ race_black, data=matched_subdf, weight=weights)
    print(paste("LOOKING FOR PH DIFFERENCES IN BIRTHS BY", dm))
    print(paste("Mean pH for Black patients:", weighted.mean(matched_subdf$pH[matched_subdf$race_black==1], w=matched_subdf$weights[matched_subdf$race_black==1])))
    print(paste("Mean pH for white patients:", weighted.mean(matched_subdf$pH[matched_subdf$race_black==0], w=matched_subdf$weights[matched_subdf$race_black==0])))
    print(summary(lm))
    print(confint.default(lm))
    


    #Get rates for this delivery method (i.e. what fraction of each delivery method has each pH category?)
    matched_subdf$categorised_pH <- purrr::map2_chr(
        matched_subdf$pH,
        matched_subdf$BaseExcess,
        divide_ph,
        .progress = TRUE
    )

    #Get the rates of pH in each race.
    matched_rates <- matched_subdf %>%
        group_by(Race) %>% mutate(OverallNumUnits=sum(weights)) %>%
        group_by(Race, categorised_pH) %>% summarise(rate=sum(weights)/first(OverallNumUnits)) %>%
        ungroup()


    #Measure and link differences.
    diffs_in_ph <- data.frame(categorised_pH=character(0), estimate=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0), pval=numeric(0))
    index <- 1
    for (categ in unique(matched_subdf$categorised_pH)) {
        copy_matched_df <- data.frame(matched_subdf)
        copy_matched_df$target <- ifelse(copy_matched_df$categorised_pH==categ, 1, 0)
        lm <- glm(target ~ race_black, data=copy_matched_df, family=binomial(link="logit"), weights=weights)
        summary <- summary(lm)
        coefs <- as.data.frame(coef(summary))
        confs <- confint.default(lm)

        estimate <- exp(coefs["race_black", 'Estimate'])
        lower_conf <- exp(confs["race_black", '2.5 %'])
        upper_conf <- exp(confs["race_black", '97.5 %'])
        pval <- as.numeric(coefs["race_black", 'Pr(>|z|)'])
        diffs_in_ph[index,] <- c(categ, estimate, lower_conf, upper_conf, pval)
        index <- index + 1
    }

    #Link this.
    matched_rates <- matched_rates %>%
        inner_join(diffs_in_ph, by="categorised_pH") %>%
        mutate(categorised_pH=factor(categorised_pH, levels = c("Acidemia", "No acidemia"))) %>%
        group_by(categorised_pH) %>% mutate(max_rate=max(rate)) %>% ungroup() %>% mutate(DeliveryMethod=dm)
    
    #Add to overall results.
    stratified_matched_rates <- rbind(stratified_matched_rates, matched_rates)
}

#Now plot likelihood of pH category by delivery method.
stratified_matched_rates <- stratified_matched_rates %>%
    mutate(signif=asterisk(pval)) 

p <- ggplot(stratified_matched_rates, aes(x=categorised_pH, y=rate*100, fill=Race)) +
    geom_col(position=position_dodge(width=0.9)) +
     labs(
        x = "pH",
        y = "Rate (% of Births by this Method)",
        fill = "Race",
        title = "pH by Delivery Method",
        subtitle = "Differences measured after exact matching of Black patients to white patients on risk factors"
        #caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    facet_wrap(~DeliveryMethod) +
    geom_text(aes(x=categorised_pH, y=rate*100+4, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=4) +
    geom_text(aes(x=categorised_pH, y=max_rate*100+7, label=signif), size=4) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
    ggsave("pH_diffs_overall_by_dm.pdf", width=12, height=6, plot=p)

##########-----PART C: DIFFERENCES IN DELIVERY METHOD BY ACIDEMIA CATEGORY.--------

stratified_matched_rates <- data.frame()
df$categorised_pH <- purrr::map2_chr(
        df$pH,
        df$BaseExcess,
        divide_ph,
        .progress = TRUE
    )


for (categ in unique(df$categorised_pH)) {
    subdf <- df %>% filter(categorised_pH==categ)

    #Match cohort.
    matched_subdf <- get_matched_cohort(subdf, "race_black", c("maternal_age", "pregnancy_term", "previous_c_section_indicated",
        "categorised_num_births"))

    #Get the rates of pH in each race (likelihood of being born by each method).
    matched_rates <- matched_subdf %>%
        group_by(Race) %>% mutate(OverallNumUnits=sum(weights)) %>%
        group_by(Race, DeliveryMethod) %>% summarise(rate=sum(weights)/first(OverallNumUnits)) %>%
        ungroup()


    #Measure and link differences.
    diffs_in_ph <- data.frame(DeliveryMethod=character(0), estimate=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0), pval=numeric(0))
    index <- 1
    for (dm in unique(matched_subdf$DeliveryMethod)) {
        copy_matched_df <- data.frame(matched_subdf)
        copy_matched_df$target <- ifelse(copy_matched_df$DeliveryMethod==dm, 1, 0)
        lm <- glm(target ~ race_black, data=copy_matched_df, family=binomial(link="logit"), weights=weights)
        summary <- summary(lm)
        coefs <- as.data.frame(coef(summary))
        confs <- confint.default(lm)

        estimate <- exp(coefs["race_black", 'Estimate'])
        lower_conf <- exp(confs["race_black", '2.5 %'])
        upper_conf <- exp(confs["race_black", '97.5 %'])
        pval <- as.numeric(coefs["race_black", 'Pr(>|z|)'])
        diffs_in_ph[index,] <- c(dm, estimate, lower_conf, upper_conf, pval)
        index <- index + 1
    }

    #Link this.
    matched_rates <- matched_rates %>%
        inner_join(diffs_in_ph, by="DeliveryMethod") %>%
        group_by(DeliveryMethod) %>% mutate(max_rate=max(rate)) %>% ungroup() %>% mutate(categorised_pH=categ)
    
    #Add to overall results.
    stratified_matched_rates <- rbind(stratified_matched_rates, matched_rates)
}


#Now plot this.
stratified_matched_rates <- stratified_matched_rates %>% 
    mutate(signif=asterisk(pval)) %>%
    #mutate(signif=asterisk(p.adjust(pval, method="fdr"))) %>%
    mutate(categorised_pH=factor(categorised_pH, levels = c("Acidemia", "No acidemia"))) 

p <- ggplot(stratified_matched_rates, aes(x=DeliveryMethod, y=rate*100, fill=Race)) +
    geom_col(position=position_dodge(width=0.9)) +
     labs(
        x = "pH",
        y = "Rate (% of Births by this Method)",
        fill = "Race",
        title = "Delivery Method by pH",
        #caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    facet_wrap(~categorised_pH, nrow=1) +
    geom_text(aes(x=DeliveryMethod, y=rate*100+4, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=4) +
    geom_text(aes(x=DeliveryMethod, y=max_rate*100+7, label=signif), size=4) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
    ggsave("dm_diffs_overall_by_pH.pdf", width=12, height=6, plot=p)



#TEST FOR DIFFERENCES IN RUPTURE TO DELIVERY TIME
#Stage 1 and 2 time basically incomprehensible (i.e. NA for almost all patients.)
#Link rupture to delivery time
unprocessed_data <- read.csv("/Volumes/chip-lacava/Groups/mfm/penn-maternal-fetal-monitoring/combined_covariates-2025-11-03.csv") %>%
     rename(baby_mrn=Baby.MRN, maternal_mrn=MRN) %>% select(baby_mrn, Rupture.to.Delivery) %>%
     mutate(baby_mrn=as.integer(baby_mrn), rupture_to_delivery=purrr::map_vec(Rupture.to.Delivery, rupture_to_delivery)) %>%
     distinct(baby_mrn, rupture_to_delivery)
write.csv(unprocessed_data, "rupture_to_delivery_time.csv")

#Now link these in the matched cohort
df <- df %>% inner_join(unprocessed_data, by="baby_mrn") %>%
    mutate(log_t = log(rupture_to_delivery))

#Look for differences in delivery time by mode of delivery
diffs_in_time <- data.frame(DeliveryMethod=character(0), estimate=numeric(0), lower_conf=numeric(0), upper_conf=numeric(0), pval=numeric(0))
stratified_matched_cohorts <- data.frame()
index <- 1
for (dm in c("C-Section (Nonreassuring EFM)", "C-Section (Labor Arrest)", "Vaginal")) {
    subdf <- df %>% filter(DeliveryMethod==dm) #%>% filter(is.finite(log_t) & !is.na(log_t)) #Filter out those with no recorded time after matching
    print(dm) #rupture to delivery time is 90% complete for everyone
    print(sum(is.finite(subdf$log_t) & !is.na(subdf$log_t))/nrow(subdf))
    subdf <- subdf %>% filter(is.finite(log_t) & !is.na(log_t)) 



    #Match cohort.
    matched_subdf <- get_matched_cohort(subdf, "race_black", c("maternal_age", "pregnancy_term", "previous_c_section_indicated",
        "categorised_num_births")) 


    #Measure and link differences.
    stratified_matched_cohorts <- rbind(stratified_matched_cohorts, matched_subdf)
    
    #Use accelerated failure time.
    lm <- glm(log_t ~ race_black, data=matched_subdf, family="gaussian", weights=weights)
    summary <- summary(lm)
    coefs <- as.data.frame(coef(summary))
    confs <- confint.default(lm)

    estimate <- exp(coefs["race_black", 'Estimate'])
    lower_conf <- exp(confs["race_black", '2.5 %'])
    upper_conf <- exp(confs["race_black", '97.5 %'])
    pval <- as.numeric(coefs["race_black", 'Pr(>|t|)'])
    diffs_in_time[index,] <- c(dm, estimate, lower_conf, upper_conf, pval)
    index <- index + 1

}

#Now plot difference in time by delivery method
diffs_in_time <- diffs_in_time %>%
    mutate(signif=asterisk(pval)) 

print(diffs_in_time) #No significant differences in time to delivery!
    #mutate(signif=asterisk(p.adjust(pval, method="fdr"))) 
p <- ggplot(stratified_matched_cohorts, aes(x=log_t, fill=Race)) +
    geom_density() +
     labs(
        x = "Rupture to Delivery Time (log)",
        y = "Density",
        fill = "Race",
        title = "Differences in Time to Delivery",
        #caption = "FTI: Fetal Tracing Indication; AD: Arrest of Descent/Dilatation"
    ) +
    facet_wrap(~DeliveryMethod) +
    #geom_text(aes(x=categorised_pH, y=rate*100+4, label=signif(rate*100, 2)), position=position_dodge(width=0.9), size=4) +
    #geom_text(aes(x=, y=max_rate*100+7, label=signif), size=4) +
    theme_minimal(base_size = 15) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())
    ggsave("time_diffs_overall_by_dm.pdf", width=12, height=6, plot=p)
