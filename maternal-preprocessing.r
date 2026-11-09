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

#Define thresholds.
arterial_threshold <- 7.10
cantu_venous_threshold <- 7.15

#-----LOAD DATASETS-----------------------------------
deliveries  <- fread(paste0(load_filepath, "Deliveries Dataset 2026.11.01.csv"))
# labs        <- read.csv(paste0(load_filepath, "Labs Dataset 2026.11.01.csv"))
pregnancies <- fread(paste0(load_filepath, "Pregnancies Dataset 2026.11.01.csv"))
newborns    <- fread(paste0(load_filepath, "Newborns Dataset 2026.11.01.csv"))

#------DEFINE COLUMNS--------------------
id_cols <- c("PregnancyId", "NewbornId", "DeliveryMrn")

covariate_cols <- c("maternal_age", "pregnancy_term", "year_of_delivery",
                "gestation_age_in_weeks",
                "race", "bmi_at_last_prenatal_visit", "FetalAnomalies",
                "parity", "prior_cd", "insurance", "prenatal_care_utilisation",
                "hypertensive_disorders", "diabetes", "substance_use_during_pregnancy", 
                "DeliveryHospital")

outcome_cols <- c(
    "delivery_method", "depressed_apgar", "neonatal_demise",
    "nicu_adm", "dx_hie", 
    "arrival_to_delivery", #Current labor duration measure
    "arterial_ph_recorded", "any_ph_recorded",
    "arterial_acidemia", "univ_acidemia", "cantu_acidemia"
)


#######------FILTER VISITS-----------------------------------

# Check that records are not duplicated between pregnancies.
assert(nrow(deliveries)==length(unique((deliveries$PregnancyId))))

# Keep singleton births only (using both the direct indication in the record and the indirect indication from pregnancies/cesarean indication.)
# Don't filter out stillbirths; this is a form of conditioning on the outcome.
deliveries <- deliveries %>%
    group_by(PregnancyId) %>% 
    mutate(DistinctFetuses=n_distinct(NewbornId)) %>%
    ungroup() %>%
    filter(NumberOfFetuses==1, DistinctFetuses==1) %>%
    select(-DistinctFetuses) %>%
    filter(!grepl("Multiple Gestation", CesareanIndication))

print(paste0("At the start of the process, we have ", nrow(deliveries), " distinct singleton deliveries."))

# Keep only those after May 2015 (when arterial pH recording starts).
deliveries <- deliveries %>% 
    filter((
        month(BirthInstant) >= 5 & year(BirthInstant)== 2015) | year(BirthInstant) >= 2016
    )

print(paste0("After filtering out births before May 2015, we have ", nrow(deliveries), " distinct singleton deliveries."))

#Filter out abortions and those with unspecified delivery method.
deliveries <- deliveries %>% 
    filter(!is.na(DeliveryMethod), !(DeliveryMethod %in% c("SAB", "TAB")))

print(paste0("After filtering out abortions, we have ", nrow(deliveries), " distinct singleton deliveries."))

######-----PROCESS COVARIATES--------------
deliveries <- deliveries %>%
    mutate(
        AgeAtDelivery = as.numeric(AgeAtDelivery),
        maternal_age = case_when(
            AgeAtDelivery < 20 ~ "19 and under",
            AgeAtDelivery < 30 ~ "20-29",
            AgeAtDelivery < 35 ~ "30-35",
            AgeAtDelivery < 40 ~ "35-39",
            AgeAtDelivery >= 40 ~ "40+",
            .default=NA
        ),
        gestation_age_in_weeks = floor(as.numeric(GestationalAgeDays)/7), #Gestational age in completed weeks
        pregnancy_term=case_when(
            gestation_age_in_weeks < 28 ~ "extremely preterm",
            gestation_age_in_weeks < 32 ~ "very preterm",
            gestation_age_in_weeks < 37 ~ "preterm",
            gestation_age_in_weeks < 39 ~ "early term",
            gestation_age_in_weeks < 41 ~ "full term",
            gestation_age_in_weeks >= 41 ~ "late term", #Not enough (post-term (post 41 weeks) to be useful)
            .default=NA  
        ),
        year_of_delivery = year(as_datetime(BirthInstant)),
        ) %>% left_join(
            select(pregnancies, 
                c(PregnancyId, DeliveryMrn, OrganizationPrenatalEncounterCount, BmiAtLastPnc, #use BMI at delivery
                DxPreeclampsia, DxSeverePreeclampsia, DxSuperimposedPreeclampsia, DxEclampsia, DxHELLP, 
                DxHtnDisorderOfPregnancy, DxGestationalHtn, DxChronicHtnDuringPregnancy, DxGestationalDiabetes, DxPregestationalDiabetes,
                PregnancyParaCount, PriorCesarean, PatientRace, Ethnicity, InsuranceClass, SmokingTobaccoUse, SmokelessTobaccoUse, AlcoholUse, 
                CocaineCrackHeroinUse, IllicitDrugUse
            )), by=c("PregnancyId", "DeliveryMrn")) %>%
        
        #Transform race using Carlson priority system.
        #141 American Indian and 40 Native Hawaiian patients; we group them with Other.
        mutate(
            race = case_when(
                Ethnicity == 'Hispanic Latino' ~ "Hispanic",
                grepl("Some Other Race|American Indian or Alaskan Native|Native Hawaiian or Other Pacific Islander", PatientRace) ~ "Other",
                grepl("Asian", PatientRace) ~ "Asian",
                grepl("Black or African American", PatientRace) ~ "Black",
                grepl("White", PatientRace) ~ "White",
                .default = "Unknown"
            ),

            #Categorise BMI according to Carlson.
            BmiAtLastPnc = as.numeric(BmiAtLastPnc),
            bmi_at_last_prenatal_visit = case_when( 
                (BmiAtLastPnc < 15) |  (BmiAtLastPnc > 100) ~ "Unknown",
                BmiAtLastPnc < 25 ~ "Normal", #As in Carlson, <20 people would be underweight, so we group them with normal weight.
                BmiAtLastPnc < 30 ~ "Overweight",
                BmiAtLastPnc < 35 ~ "Obese I",
                BmiAtLastPnc < 40 ~ "Obese II",
                BmiAtLastPnc >= 40 ~ "Obese III",
                .default = "Unknown"
            ),

            #Birth history (Prior C-section often null for nulliparous laborers).
            parity = ifelse(PregnancyParaCount > 0, 1, 0),
            prior_cd = ifelse(!is.na(parity) & !is.na(PriorCesarean) & PriorCesarean == 1, 1, 0), #

            #Insurance
            insurance = case_when(
                InsuranceClass == "Commercial" ~ "private",
                InsuranceClass %in% c("Medicaid", "Medicare", "Managed Medicaid") ~ "public",
                .default = "Other/Unknown" #includes TRICARE, VA, self-pay, or unknown
            ),

            #Prenatal care.
            #10th percentile of visits is 0, 25th is 5, median is 8, 75th is 11, 90th is 14
            #Use of the Kotelchuck visits is difficult because we don't know when the prenatal care began
            #ACOG guidelines recommend visits every 4 weeks, every 2 weeks until 36 weeks, and weekly thereafter
            #Use Carter 2016 guidelines.
            prenatal_care_utilisation = case_when(
                OrganizationPrenatalEncounterCount <= 5 ~ "low",
                OrganizationPrenatalEncounterCount <= 10 ~ "medium",
                OrganizationPrenatalEncounterCount > 10 ~ "high",
                .default = "unknown"
            ),

            #Pre-existing diagnoses
            hypertensive_disorders = ifelse(
                DxPreeclampsia + DxSeverePreeclampsia + DxSuperimposedPreeclampsia + DxEclampsia +
                DxHELLP + DxHtnDisorderOfPregnancy + DxGestationalHtn + DxChronicHtnDuringPregnancy > 0,
                1, 0
            ),
            diabetes = ifelse(
                DxGestationalDiabetes + DxPregestationalDiabetes > 0, 1, 0
            ),

            #Confirmed substance use during pregnancy
            substance_use_during_pregnancy = ifelse(
                SmokingTobaccoUse %in% c("Every Day", "Heavy Smoker", "Light Smoker", "Some Days") | #Around 600
                #excludes 'Smoker, Current Status Unknown', 'Unknown', 'Former', 'Never', 'Never Assessed', 'Passive Smoke Exposure - Never Smoker', null 
                SmokelessTobaccoUse %in% c('Current') | #Around 39
                #excludes 'Former', 'Never', 'Unknown', null
                AlcoholUse %in% c('Yes') | #2300
                #excludes 'Never', 'No', "Not Asked", 'Not Currently', and null
                CocaineCrackHeroinUse %in% c('Yes') | # 45
                #excludes null
                IllicitDrugUse %in% c('Yes'), #785
                #excludes 'Never', 'No', "Not Asked", "Not Currently", and "No"
                1, 0
            )) 

###---------FILTER FOR KNOWN COVARIATES---------
deliveries <- deliveries %>%
    filter(!is.na(maternal_age), !is.na(pregnancy_term), !is.na(DeliveryMethod))

print(paste0("After filtering out those with unknown maternal age, gestational age or delivery method, we have ", nrow(deliveries), " distinct singleton deliveries."))


#######------ATTACH OUTCOMES-----------------------------------

# NOTE: the 'labs' dataset does not include pH that I can see-- it's mostly hemoglobin.

#First, load outcomes.
outcomes <- newborns %>% rename(DeliveryMrn=MomDeliveryMrn) 
delivery_times <- deliveries %>% select(BirthInstant, PregnancyId, DeliveryMrn, NewbornId)


#Extract JSON measurements.
get_json_measurements <- function(outcomes, delivery_times, measurement='CordpH') {

    #Extract from JSON format.
    cord_gas <- outcomes %>% select(DeliveryMrn, PregnancyId, NewbornId, !!sym(measurement)) %>% 
        mutate(!!sym(measurement) := map(!!sym(measurement), ~ {
            if (is.na(.x)) NULL 
            else fromJSON(gsub('""', '"', .x))
        })) %>%
        unnest(!!sym(measurement)) %>% 
        mutate(!!sym(measurement) := as.numeric(!!sym(measurement))) %>%
        filter(!is.na(!!sym(measurement)))

    return(cord_gas)
}

# Get first cord gas and base excess measurements, stored in JSON format.
get_first_gas_measurements <- function(outcomes, delivery_times, measurement='CordpH') {

    #Apply JSON extraction function.
    cord_gas <- get_json_measurements(outcomes, delivery_times, measurement=measurement)

    #Get first measurements.
    cord_gas <- cord_gas %>% 
        inner_join(delivery_times, by=c("PregnancyId", "NewbornId", "DeliveryMrn")) %>%
        mutate(
            CollectionInstant=as_datetime(CollectionInstant),
            BirthInstant=as_datetime(BirthInstant),
            cord_gas_interval = as.numeric(CollectionInstant-BirthInstant, units='mins')) %>%
        #Keep only measurements within the 30-min timeframe.
        filter(cord_gas_interval <= 30 & cord_gas_interval >= 0) %>% 
        #Drop non-arterial measurements where arterial ones are available.
        group_by(PregnancyId, DeliveryMrn, NewbornId) %>%
        mutate(
            arterial_available = any(Source=='Arterial')) %>%
        ungroup() %>%
        filter(
            (arterial_available & (Source == 'Arterial')) | (!arterial_available)
        ) %>%
        #Take the first of the remaining measurements.
        group_by(PregnancyId, DeliveryMrn, NewbornId) %>%
        filter(cord_gas_interval == min(cord_gas_interval)) %>% 
        ungroup() %>%
        #If multiple measurements at the same time are available, prefer the venous one; if not, break ties at random.
        mutate(is_venous = (Source == 'Venous')) %>%
        group_by(PregnancyId, DeliveryMrn, NewbornId) %>%
        arrange(PregnancyId, DeliveryMrn, NewbornId, desc(as.integer(is_venous))) %>%
        dplyr::slice(1) %>%
        ungroup() %>%
        select(DeliveryMrn, PregnancyId, NewbornId, Source, !!sym(measurement)) %>%
        rename(!!sym(paste0(measurement, "_Source")) := Source) 

    return(cord_gas)
}

#---Extract cord gas and base excess measurements.
cord_gas <- get_first_gas_measurements(outcomes, delivery_times)
base_excess <- get_first_gas_measurements(outcomes, delivery_times, measurement='CordBaseExcess')


#---Extract other outcomes from the newborn dataset.
newborn_outcomes <- outcomes %>%
    rename(nicu_los = NicuLengthOfStay) %>%
    mutate(
        nicu_adm = ifelse(HadNicuAdmission + HadNicuReadmission > 0, 1, 0),
        dx_hie = !is.na(DxHie)) %>%
    select(PregnancyId, DeliveryMrn, NewbornId, nicu_los, nicu_adm, dx_hie, ExpiredInNicu)

#----Extract outcomes from the pregnancy dataset.
pregnancy_outcomes <- pregnancies %>% 
    mutate(
        fetal_or_neonatal_demise = ifelse(HadFetalDemise + HadNeonatalDemise > 0, 1, 0)
        ) %>%
    select(PregnancyId, DeliveryMrn, fetal_or_neonatal_demise, PlannedDeliveryMethod)

# #---Join these outcomes to the delivery dataset.
deliveries <- deliveries %>% left_join(
    newborn_outcomes, by=c("PregnancyId", "DeliveryMrn", "NewbornId")
) %>% left_join(
    pregnancy_outcomes, by=c("PregnancyId", "DeliveryMrn") 
) %>% left_join(
    cord_gas, by=c("PregnancyId", "DeliveryMrn", "NewbornId")
) %>% left_join(
    base_excess, by=c("PregnancyId", "DeliveryMrn", "NewbornId")
)

assert(nrow(deliveries)==length(unique(deliveries$NewbornId)))

#----ASSESS AVAILABILITY OF CORD GAS DATA------

hosp_cord_gas_dist <- deliveries %>%
    group_by(DeliveryHospital) %>%
    mutate(NumBirths = n()) %>%
    group_by(DeliveryHospital, CordpH_Source) %>%
    summarise(RecordingRate=round(n()/first(NumBirths)*100, 1), RecordingN=n()) %>%
    filter(!is.na(CordpH_Source))
write.csv(hosp_cord_gas_dist, "pH_recording_rates_by_hospital.csv")

hosp_base_excess_dist <- deliveries %>%
    group_by(DeliveryHospital) %>%
    mutate(NumBirths = n()) %>%
    group_by(DeliveryHospital, CordBaseExcess_Source) %>%
    summarise(RecordingRate=round(n()/first(NumBirths)*100, 1), RecordingN=n()) %>%
    filter(!is.na(CordBaseExcess_Source))
write.csv(hosp_base_excess_dist, "base_excess_recording_rates_by_hospital.csv")

#------FILTER FOR HOSPITALS-------
deliveries <- deliveries %>%
    filter(DeliveryHospital %in% c("HUP", "PAH"))

print(paste0("After filtering out those outside HUP/PAH, we have ", nrow(deliveries), " distinct singleton deliveries."))


#---Extract delivery method (from DeliveryMethod; if unspecified, use CesareanDelivery column to confirm C-sections).
# #Births are assumed to be C-sections under any of the following conditions:
    # #C-section actively documented in DeliveryMethod;
    # #CesareanDeliveryCesareanDelivery column set to 1 (indicating positive documentation of a cesarean);
    # #C-section indications documented;
    # #General or spinal anesthesia given.
    # All other births are assumed to be vaginal (i.e. CDCD ==0 is assumed to indicate vaginal birth.)
deliveries <- deliveries %>%
    mutate(
        delivery_type = case_when(
            DeliveryMethod %in% c(
                "C-Section, Classical",
                "C-Section, Low Transverse",
                "C-Section, Low Vertical",
                "C-Section, Unspecified",
                "Caesarean Section") | 
                CesareanDeliveryCesareanDelivery == 1 |
                (!is.na(CesareanIndication) & CesareanIndication != 'null' & CesareanIndication != '') |
                grepl('Spinal', BirthAnesthesia, ignore.case = TRUE) |
                grepl('General', BirthAnesthesia, ignore.case = TRUE)
                ~ "C-Section",
            DeliveryMethod %in% c(
                "FAVD",
                "Forceps Assisted Delivery",
                "SVD*",
                "Vacuum Assisted Delivery",
                "Vaginal Delivery",
                "Vaginal Delivery`",
                "Vaginal, Breech",
                "Vaginal, Forceps",
                "Vaginal, Spontaneous",
                "Vaginal, Vacuum (Extractor)",
                "VBAC, Spontaneous") |
                CesareanDeliveryCesareanDelivery == 0
                ~ "Vaginal",
            .default = 'Unknown'
            ),
        #Use the 2025 rules: fetal tracing carries, then labor indication; all other C-sections assumed to be planned.
        delivery_method = case_when(
            delivery_type=='C-Section' & grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", CesareanIndication) ~ 'CD (NRFHT)',
            delivery_type=='C-Section' & grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", CesareanIndication) ~ "CD (Labor Arrest)",
            delivery_type=='C-Section' ~ 'CD (Planned)',
            delivery_type=='Vaginal' ~ 'Vaginal',
            .default = 'Unknown'
        )
        ) 


#---MEASURE DELIVERY METHOD INCONSISTENCIES------

# #Check-- cross-compare indications and indications.
write.csv(table(deliveries$CesareanIndication, deliveries$delivery_method), 'checks/crosstable-methods-and-indication.csv')
write.csv(table(deliveries$LaborAttempted, deliveries$delivery_method), 'checks/crosstable-methods-and-labor-attempt.csv')
write.csv(table(deliveries$CesareanExpected, deliveries$delivery_method), 'checks/crosstable-methods-and-cd-expected.csv')

# #Measure inconsistencies between delivery methods and other checks, using the previous set of rules.
overlaps <- deliveries %>% 
    group_by(delivery_method) %>%
    summarise(
        nrfht_indicated = round(mean(grepl("Fetal Tracing Indication|Cord Prolapse|Prolapsed Cord|Suspected Uterine Rupture", CesareanIndication))*100),
        labor_arrest_indicated = round(mean(grepl("Failed Induction|Arrest of Descent|Arrest of Dilatation", CesareanIndication))*100),
        repeat_elective_indicated = round(mean(grepl("Repeat Elective", CesareanIndication))*100),
        labor_attempt_documented = round(mean(!is.na(LaborAttempted) & LaborAttempted==1)*100),
        absence_of_labor_attempt_documented = round(mean(!is.na(LaborAttempted) & LaborAttempted==0)*100),
        expected_c_section_documented = round(mean(!is.na(CesareanExpected) & CesareanExpected==1)*100)
    )
write.csv(overlaps, 'checks/methods-crosschecks.csv')

#---Define final outcomes from the delivery method.

deliveries <- deliveries %>% 
    mutate(
        depressed_apgar = ifelse(TotalApgarFiveMinute < 7, 1, 0),
        neonatal_demise = ifelse(
            fetal_or_neonatal_demise == 1 |
            BornAlive == 0 | LivingStatus %in% c("Fetal Demise", "Neonatal Demise") |
            ExpiredInNicu == 1, 
            1, 0
        ),

        BirthInstant = as_datetime(BirthInstant),
        PushingStartInstant = as_datetime(PushingStartInstant),
        MotherArrivalInstant = as_datetime(MotherArrivalInstant),
        LaborStartInstant = as_datetime(LaborStartInstant),

        #Labor duration measures.
        second_stage_duration = as.numeric(SecondStageLengthMinutes),
        first_stage_duration = as.numeric(FirstStageLengthMinutes),
        first_and_second_stages_duration = first_stage_duration + second_stage_duration,
        implied_labor_duration = as.numeric(BirthInstant - LaborStartInstant, units='mins'),
        pushing_to_delivery = as.numeric(BirthInstant - PushingStartInstant, units='mins'),
        arrival_to_delivery = as.numeric(BirthInstant - MotherArrivalInstant, units='mins'),
    )


#-----ASSESS CORD GAS DISTRIBUTIONS-----
cord_gas_dist <- ggplot(
    deliveries, aes(x=CordpH, fill=CordpH_Source)
) + geom_density(alpha = 0.3) + 
    labs(x='Cord gas pH', y='Probability density', fill='Source')
    theme_minimal() +
    theme(legend.position = "bottom")
ggsave('cord-gas-dist.pdf', width=6, height=4, cord_gas_dist, dpi=300)

any_ph_recorded <- mean(!is.na(deliveries$CordpH_Source))
arterial_ph_recorded <- mean(!is.na(deliveries$CordpH_Source) & (deliveries$CordpH_Source=='Arterial'))

print(paste0("Of these, we have ", round(any_ph_recorded*100, 1), "% recording rate in pH; for arterial, this is ", round(arterial_ph_recorded*100, 1), "%."))


#---DEFINE CORD GAS OUTCOMES-------
deliveries <- deliveries %>%
    mutate(

        #Note whether any reading exists.
        arterial_ph_recorded = ifelse(!is.na(CordpH_Source) & CordpH_Source == 'Acidemia', 1, 0),
        any_ph_recorded = ifelse(!is.na(CordpH_Source), 1, 0),

        #ACIDEMIA.From most to least conservative:

        #Just use arterial cord gas pH readings.
        arterial_acidemia = ifelse(
            arterial_ph_recorded, ifelse(CordpH <= arterial_threshold, 1, 0), NA_real_
        ),

        #Apply arterial thresholds to all-source pH readings.
        univ_acidemia = ifelse(
            any_ph_recorded, 
            ifelse(CordpH <= arterial_threshold, 1, 0), NA_real_
        ),
        
        #Apply Cantu thresholds to venous and unspecified cord pH readings.
        cantu_acidemia = ifelse(
            arterial_ph_recorded, arterial_acidemia, 
            ifelse(
                any_ph_recorded, 
                ifelse(CordpH <= cantu_venous_threshold, 1, 0), NA_real_)
        )


    )

#----FILTER FOR COVARIATES WE CARE ABOUT AND SAVE-------
deliveries <- deliveries %>% 
    select(all_of(c(id_cols, covariate_cols, outcome_cols)))
for (col in colnames(deliveries)) {
    print(col)
    assert(sum(is.na(deliveries$col))==0)
}
write.csv(deliveries, "preprocessed-deliveries-with-outcomes.csv")





