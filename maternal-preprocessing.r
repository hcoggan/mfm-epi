
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
setwd("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/ch265493/mfm-penn/epi-work/mar-res-main/")
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

#Does she have any previous Cesarean Deliverys indicated indirectly?
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
            (births$delivery_method %in% c("Cesarean Delivery (unplanned, with FTI)", "Cesarean Delivery (planned)", "Cesarean Delivery (unplanned, without FTI)")) | #Either THIS birth was by Cesarean Delivery, or
            (births$previous_c_section_indicated_directly == 1)) #Cesarean Delivery has already taken place by then
        )
    return(ifelse(previous>0, 1, 0))
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

# # #######----------PREPROCESS COVARIATES----------

#Load raw data.
births <- read.csv("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/mfm/penn-maternal-fetal-monitoring/combined_covariates-2025-11-03.csv") %>%
    rename(baby_mrn=Baby.MRN, maternal_mrn=MRN) #Initially: 83953 unique births, 68110 unique maternal MRNs- roughly 1.23 births per mother


#Define necessary covariates.
births <- births %>% 
    distinct(baby_mrn, .keep_all = TRUE) %>% #Get rid of multiple records per baby.- now 68110 unique maternal MRNs (83953 unique births)
     mutate(

        apgar_score = as.numeric(Apgar.5), #5-minute Apgar score--in final cohort, 65 have none recorded (0.75% of final cohort)
        depressed_apgar = ifelse(!is.na(apgar_score) & (apgar_score < 7), 1, 0),

         gestation_age_in_weeks=purrr::map_vec(GA, gest_age),

         pregnancy_term=case_when(
            gestation_age_in_weeks < 28 ~ "extremely preterm",
            gestation_age_in_weeks < 32 ~ "very preterm",
            gestation_age_in_weeks < 37 ~ "preterm",
            gestation_age_in_weeks < 39 ~ "early term",
            gestation_age_in_weeks < 41 ~ "full term",
            gestation_age_in_weeks >= 41 ~ "late term", #Not enough (post-term (post 41 weeks) to be useful)
            .default=NA     
        ),

        race=case_when(
            Race %in% c('Black', 'White', 'Asian', 'American Indian', 'Pacific Island', 'East Indian') ~ Race,
            Race == 'HLW-Hispanic Latino/White' ~ "White",
            Race == 'HLB-Hispanic Latino/Black' ~ "Black",
            Race == "" | is.na(Race) | is.null(Race) | grepl("Patient Declined", Race) | grepl("Unknown", Race) ~ "Unknown",
            .default = "Other"
        ),

        #Record the presence of a previous cesarean
        previous_c_section_indicated_directly=ifelse(
            Delivery.Method=="VBAC" | Delivery.Method.1 == "VBAC" | (!is.na(as.numeric(Prior.C.S.Ct)) & (as.numeric(Prior.C.S.Ct)>0)) , 1, 0),

        delivery_method = case_when(
            Delivery.Method %in% c('VBAC, Spontaneous', 
                'Vaginal, Spontaneous', 
                'Vaginal, Forceps',
                'Vaginal, Vacuum (Extractor)', 
                'Vaginal Delivery', 
                'Vaginal, Breech', 
                'SVD*',
                'FAVD', 
                'Vaginal Delivery`', 
                'Vacuum Assisted Delivery',
                'Forceps Assisted Delivery') |
            Delivery.Method.1 %in% c('VBAC, Spontaneous', 
                'Vaginal, Spontaneous', 
                'Vaginal, Forceps',
                'Vaginal, Vacuum (Extractor)', 
                'Vaginal Delivery', 
                'Vaginal, Breech', 
                'SVD*',
                'FAVD', 
                'Vaginal Delivery`', 
                'Vacuum Assisted Delivery',
                'Forceps Assisted Delivery') ~ "Vaginal",

            Delivery.Method %in% c(
                'C-Section, Classical',
                'C-Section, Low Transverse',
                'C-Section, Unspecified',
                'C-Section, Low Vertical', 
                'Caesarean Section') |
            Delivery.Method.1 %in% c(
                'C-Section, Classical',
                'C-Section, Low Transverse',
                'C-Section, Unspecified',
                'C-Section, Low Vertical', 
                'Caesarean Section') |
            !(C.S.Categorization == "") |
            !(C.S.Priority == "") |
            !(Indications.for.Cesarean == "") |
            !(C.Section.Indications== "") |
            grepl("general", Anesthesia) |
            grepl("spinal", Anesthesia)
             ~ "C-Section",
            .default = "Other/Unknown"
            ),


        #Temporal variables
        Delivery.Date.Time=ymd_hms(Delivery.Date.Time),
        original_delivery_date=Delivery.Date.Time,
        arbitrary_timestamp=as.numeric(difftime(Delivery.Date.Time, ymd_hms("2000-01-01 00:00:00"), units="mins")), #Define each birth by an arbitrary timestamp 
        year_of_delivery=year(Delivery.Date.Time),
        month_of_delivery=month(Delivery.Date.Time),
        season_of_delivery=case_when(
            month_of_delivery > 2 & month_of_delivery <= 5 ~ "spring",
            month_of_delivery > 5 & month_of_delivery <= 8 ~ "summer",
            month_of_delivery > 8 & month_of_delivery <= 11 ~ "autumn",
            month_of_delivery > 11 | month_of_delivery <= 2 ~ "winter",
            .default = NA
        ),
        is_weekend=ifelse(wday(Delivery.Date.Time) %in% c(6, 7), 1, 0),
        hour_of_delivery=hour(Delivery.Date.Time),
        time_of_delivery = case_when(
            hour_of_delivery <= 5 ~ "small hours",
            hour_of_delivery <= 11 ~ "morning",
            hour_of_delivery <= 17 ~ "afternoon",
            hour_of_delivery <= 23 ~ "evening",
            .default = NA
        ),

        #Length of stage 1 (before 10cm dilation) and stage 2 (after)
        stage_1_length = purrr::map_vec(Stage.1.Length, labour_stage),
        stage_2_length = purrr::map_vec(Stage.2.Length, labour_stage),

        #Location of delivery
        location=case_when(
            Delivery.Location %in% c("CHESTER COUNTY HOSPITAL", "HUP", "MEDICAL CENTER OF PRINCETON", "PENNSYLVANIA HOSPITAL") ~ Delivery.Location,
            .default="Other/Unknown"
        ),


        #If someone does get a C-section, note the presence of specific explanations:
        fetal_tracing_indicated=ifelse(grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", Indications.for.Cesarean) | grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", C.Section.Indications), 1, 0),
        other_unplanned=ifelse(grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", Indications.for.Cesarean) | grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", C.Section.Indications), 1, 0),

        #Note stillbirths
        stillbirth_recorded = ifelse(
            Neonatal.Demise.=="Yes" |
            Liv.Stat %in% c("ND", "FD"), 1, 0), #'living at delivery' seems to indicate something different.

        #Now categorise maternal age.
        maternal_age=case_when(
            as.numeric(Age..Years.) < 20 ~ "19 and under",
            as.numeric(Age..Years.) < 30 ~ "20-29",
            as.numeric(Age..Years.) < 35 ~ "30-35",
            as.numeric(Age..Years.) < 40 ~ "35-39",
            as.numeric(Age..Years.) >= 40 ~ "40+",
            .default=NA
        ),

        #Keep raw maternal age for later.
        raw_maternal_age=as.numeric(Age..Years.),
        
        #Extract number of previous pregnancies, and previous parity
        num_previous_pregnancies=purrr::map_vec(GP, previous_pregnancies), 
        num_parity=purrr::map_vec(GP, parity)

        ) %>% #Still 68110 mothers.

        #A check against multiple births: discard births where there's another in the last 30 days.
        arrange(maternal_mrn, arbitrary_timestamp) %>% group_by(maternal_mrn) %>%
        mutate(
            num_births_within_1_month = purrr::map_dbl(row_number(), function(i) {

                current_time <- arbitrary_timestamp[i]
                current_baby <- baby_mrn[i]

                sum(!(baby_mrn==current_baby) &
                     (abs(arbitrary_timestamp-current_time) < 24*60*30)) #Count any birth of a DIFFERENT baby within 30 days to the same mother.
             })) %>% ungroup() %>% filter((year_of_delivery == 2017 & month_of_delivery >=4) | year_of_delivery > 2017, location=="HUP") 

#Identify term births.
births <- births %>%
    filter(gestation_age_in_weeks > 20)

print(paste("Initially:", length(unique(births$baby_mrn)), "births at HUP after April 2017 with a GA of more than 20 weeks")) 

births <- births %>%
    filter(num_births_within_1_month==0) #%>% #Now 66813 records- around 2% of maternal MRNs lost, around 3% of births are multiple, fine- still 1.20 births per mother

print(paste("After filtering out births to the same mother < 1 month apart:", length(unique(births$baby_mrn)), "births")) 

births <- births %>%
        select(baby_mrn, maternal_mrn, gestation_age_in_weeks, race, delivery_method, year_of_delivery,
            season_of_delivery, month_of_delivery, is_weekend, hour_of_delivery, time_of_delivery, location, 
             fetal_tracing_indicated, other_unplanned, stillbirth_recorded,
             maternal_age,  num_previous_pregnancies, num_parity, apgar_score, depressed_apgar,
             stage_1_length, stage_2_length, arbitrary_timestamp, original_delivery_date, raw_maternal_age, pregnancy_term,
                previous_c_section_indicated_directly, Induction, Augmentation, L.D.Complications, Rupture.Type, Forceps.Attempted, Vacuum.Attempted) #%>% 

#print(paste("After filtering out births not at HUP:", length(unique(births$baby_mrn)), "births")) 



#----EXTRACT INDUCTIONS/AUGMENTATIONS/COMPLICATIONS------

#A function to extract multi-valued attributes-- labor complications, induction, and so on.
extract_vars <- function(birth, col, col_name) {

    #Written for the example of induction.

    #Each induction becomes a vector of 1s and 0s (i.e. induction_Oxytocin...)
    inductions <- unique(unlist(str_split(births[[col]], '\n')))
    print(inductions)
    induction_cols <- c()
    for (induction in inductions) {
        if (!((induction=='None') | (induction=='') | (induction=='None.'))) { #i.e. filter out blank strings
            births[[paste0(col_name, '_', induction)]] <- ifelse(grepl(induction, births[[col]]), 1, 0)
            induction_cols <- c(induction_cols, paste0(col_name, '_', induction))
        }
    }
    #Also add an any-induction column.
    births[[paste0('any_', col_name)]] <- ifelse(rowSums(births[induction_cols]) > 0, 1, 0)
    return(births)
}

#Extract everything we need.
births <- extract_vars(births, 'Induction', 'induction')
births <- extract_vars(births, 'Augmentation', 'augmentation')
births <- extract_vars(births, 'L.D.Complications', 'complication')



#----CATEGORISE DELIVERY METHOD/LABOR ONSET-------

births <- births %>%
             mutate(delivery_method=case_when(
                delivery_method=="C-Section" & fetal_tracing_indicated ~ "Cesarean Delivery (FHT)",
                delivery_method=="C-Section" & other_unplanned ~ "Cesarean Delivery (Labor Arrest)",
                delivery_method=="C-Section" ~ "Cesarean Delivery (Planned)",
                delivery_method=="Vaginal" ~ "Vaginal",
                .default = "Other/Unknown") #Classify delivery methods.
                , labor_onset = case_when(
                    delivery_method=='Cesarean Delivery (Planned)' ~ 'None',
                    any_induction==1 ~ 'Induced',                               #Excludes planned C-sections
                    any_augmentation==1 ~ 'Spontaeneous (With Augmentation)',   #Excludes induction
                    any_augmentation==0 ~ 'Spontaeneous (Without Augmentation)',
                    .default = 'Other/Unknown'
                )
            )


#Save raw covariates.
write.csv(births, "preprocessed-births-without-outcomes.csv")

# #########-----------LINK ACIDEMIA OUTCOMES----------------------

#Load preprocessed visits.
births <- read.csv("preprocessed-births-without-outcomes.csv") 


#Refine maternal MRNs. (We have 14270 unique maternal MRNs beforehand, and afterwards; we're not filing anything down.)
births <- births %>% mutate(maternal_mrn=
    purrr::map_vec(as.character(maternal_mrn), ~tryCatch(process_maternal_mrns(.x), error=function(e) {NA})))

#Load offsets from individual files.
offsets <- data.frame()
for (k in 1:9) {
    df <- read.csv(paste0("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/mfm/penn-maternal-fetal-monitoring/id_map/CTGData_offsets/0", k, "_offsetmap.csv"))
    offsets <- rbind(offsets, df)
}
offsets <- offsets %>%
        rename(filename=file)

#Link trace-end timestamps to maternal MRNs, joining the ID-map to the trace file on filenames.
id_map <- read.csv("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/mfm/penn-maternal-fetal-monitoring/id_map/Project Data trace_start_end_ctgdata_MASTER identified dataset 124776 tracings JM merge 12.5.22.csv") %>%
    mutate(filename=purrr::map_vec(filename, function(x){str_split(x, "/")[[1]][2]})) %>%
    rename(maternal_mrn=Mother.MRN) %>%
    select(maternal_mrn, filename, trace_end) %>% #Keep the trace_end in order to number births consistently. Assume that each trace_end uniquely identifies a birth.
    left_join(offsets, by="filename")

#Now load lab results and link them to maternal MRNs.
lab_results <- read.csv("/Users/helenacoggan/mnt/chip-lacava-e2/Groups/mfm/penn-maternal-fetal-monitoring/CTGData_labs.csv") %>%
    inner_join(id_map, by="filename") %>%
    mutate(parsed_trace_end_date=parse_date_time(trace_end, "mdy HM"), parsed_order_date=ymd_hms(LabOrderDtime)) %>% #1.129 filenames per mother; 1.101 valid delivery dates per mother; 1.283 unique order timestamps per mother.
    filter(LabTestName %in% c("pH Cord", "Base Excess Cord")) %>% #1.073 filenames per mother; 1.04 valid delivery dates per mother; 1.106 unique order timestamps per mother.

    #If there are multiple labs of the same type ordered at the same timestamp, discard them as a proxy for multiple births.
    group_by(maternal_mrn, parsed_order_date, LabTestName) %>% #1.073967 filenames/1.043493 delivery dates/1.106282 order timestamps
    filter(n()==1) %>% 
    ungroup() %>% 
    
    #If there are multiple labs of the same type ordered at DIFFERENT timestamps, take the first one.
    group_by(maternal_mrn, parsed_trace_end_date, LabTestName) %>%
    filter(parsed_order_date==first(parsed_order_date)) %>% 
    ungroup() %>% #1.072218 filenames, 1.043493 order dates, 1.068346 order timestamps

    tidyr::pivot_wider(id_cols=c(maternal_mrn, parsed_trace_end_date, parsed_order_date, filename, offset_hours), names_from=LabTestName, values_from=LabResultValueFloat) %>%
    clean_names() 


#Now link results to births.
births <- inner_join(births, lab_results, by="maternal_mrn", relationship="many-to-many") %>% #1.08944357774311 births per mother; 1.07106950944704 delivery dates per mother; 1.07115617958052 order dates per mother."
    ungroup() %>% arrange(maternal_mrn) %>%
    mutate(parsed_original_delivery_date=ymd_hms(original_delivery_date), offset_hours=ifelse(is.na(offset_hours), 0, offset_hours)) %>%
    mutate(parsed_trace_end_date_offset = parsed_trace_end_date + hours(offset_hours),
        parsed_original_delivery_date_offset = parsed_original_delivery_date + hours(offset_hours),
        parsed_order_date_offset = parsed_order_date + hours(offset_hours)) %>%
   mutate(mins_from_end_of_tracing_to_delivery=abs(as.numeric(difftime(parsed_original_delivery_date, parsed_trace_end_date, units="mins")))) 
   
print(paste("After filtering out births without unique cord pH and cord base excess lab orders (where unique means one timestamp per lab type):", length(unique(births$baby_mrn)), "births")) 

births <- births %>%
   filter(mins_from_end_of_tracing_to_delivery <= 24*60*30) #%>% #Match within 1 month-- as the tracing often ends a few minutes or hours before the birth. #1.016 births per mother.


births <- births %>%
   group_by(baby_mrn) %>% filter(mins_from_end_of_tracing_to_delivery==min(mins_from_end_of_tracing_to_delivery)) %>% #Take the closest pair of tracing-end/delivery-date if there are multiple in a month
   ungroup() %>%
   distinct(maternal_mrn, baby_mrn, .keep_all = TRUE)

print(paste("After filtering out births with lab results more than 1 month from delivery:", length(unique(births$baby_mrn)), "births")) 




#Keep the first birth associated with every patient.
births <- births %>% mutate(recorded=1)
raw_births <- read.csv("preprocessed-births-without-outcomes.csv") %>%
    left_join(select(births, c(baby_mrn, recorded)), by="baby_mrn") %>%
    arrange(arbitrary_timestamp) %>%
    group_by(maternal_mrn) %>% mutate(birth_number=row_number()) %>%
    filter(recorded==1) %>% ungroup()

births <- births %>%
    inner_join(select(raw_births, baby_mrn, birth_number), by="baby_mrn") %>%
    filter(birth_number==1)

print(paste("After taking the first birth to each mother:", length(unique(births$baby_mrn)), "births")) 


#Identify whether patients have records of previous C-sections.
births$previous_c_section_indicated_indirectly <- purrr::map_vec(1:nrow(births), ~identify_previous_c_sections(births, .x), .progress=TRUE)
births$previous_c_section_indicated <- ifelse(births$previous_c_section_indicated_directly + births$previous_c_section_indicated_indirectly > 0, 1, 0)

#Categorise acidemia outcomes.
births <- births %>%
    mutate(p_h_cord = as.numeric(p_h_cord), base_excess_cord = as.numeric(base_excess_cord),
        acidemia = ifelse((p_h_cord <= pH_threshold) & (base_excess_cord < base_excess_threshold), 1, 0),
        acidemia_none = ifelse((p_h_cord <= pH_threshold) & (base_excess_cord < base_excess_threshold), 0, 1)) %>% #Main acidemia def is just pH, with base excess as sens. analysis
    mutate(gestation_age_in_weeks=as.numeric(gestation_age_in_weeks)) %>%
    filter(!is.na(gestation_age_in_weeks)) #%>% #Discards 30 births.


print(paste("After filtering out births with no recorded gestational age:", length(unique(births$baby_mrn)), "births")) 


births <- births %>%
    filter(!is.na(p_h_cord)) %>% #Discards 16 births
    filter(!is.na(base_excess_cord)) #%>% #Discards 3 births 

print(paste("After filtering out births with either a non-existent pH or base excess results:", length(unique(births$baby_mrn)), "births")) 


births <- births %>%
    mutate(race=case_when(
        race %in% c("Asian", "East Indian") ~ "Asian", #Groups together 75 East Indian births
        race %in% c("American Indian", "Pacific Island") ~ "Other", #Groups together 18 births
        .default = race
    ), pregnancy_term = case_when(
        pregnancy_term=="extremely preterm" ~ "very preterm", #only 70 births in this category
        .default=pregnancy_term 
    )) %>%
    filter(!(delivery_method=="Other/Unknown")) #Removes 1 birth.

print(paste("After filtering out births with unknown delivery method:", length(unique(births$baby_mrn)), "births")) 

births <- births %>%
    mutate(categorised_num_pregnancies = ifelse(is.na(num_previous_pregnancies), NA, ifelse(as.numeric(num_previous_pregnancies) > 1, "2+", as.character(num_previous_pregnancies))),
            categorised_num_births = ifelse(is.na(num_parity), NA, ifelse(as.numeric(num_parity) > 0, 1, 0))) %>% #Categorises continuous variables.
        filter(!is.na(categorised_num_births)) %>%
    clean_names() 

print(paste("After filtering out births with unknown parity:", length(unique(births$baby_mrn)), "births")) 

print(paste("After attaching lab results, we have ", nrow(births), "records.")) #8703 records

#Save preprocessed births with outcomes.
write.csv(births, "preprocessed-births-with-outcomes.csv")
