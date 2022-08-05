library(tidyverse)
library(dplyr)

# TODO: Have this script be like a data formatter that can be added to egfr plotter github repo
# TODO: Some data loss selecting DIST_ID over Patient_ID, could change that.
# Basically, you want this script to produce 4 files like Scott's listed below,
# based off of all the data thats currently in the database.

# Scott's OG files for reference
# master.file = read_tsv("Y:/IRB_00096551/egfr-table-files/slope_cohort_summary_right_censored_vitals_labs_meds_insulin.tsv")
master.file <- master.file.temp
bmi.file = read_tsv("Y:/IRB_00096551/egfr-table-files/bmi_study_pop.tsv")
bp.file = read_tsv("Y:/IRB_00096551/egfr-table-files/bp_study_pop.tsv")
outpatient.egfr.3om = read_tsv("Y:/IRB_00096551/egfr-table-files/slope_cohort_egfr.tsv")

ids <- read_csv("Y:/IRB_00096551/UDDB_Table_Files_For_Plotting_2021/uddb_patientID_2021.csv", col_types = "cnccc") %>%
  separate(UPDB_ID, c(NA, "UPDB_ID"), sep = "_") %>%
  separate(DIST_ID, c(NA, "DIST_ID"), sep = "_")


################################################################################
#get bmi file
################################################################################
# Reading in BMI file pulled straight from the database
currentBMI = read_csv("Y:/IRB_00096551/UDDB_Table_Files_For_Plotting_2021/uddb_bmi_height_weight.csv")
currentBMI$BMI = as.numeric(currentBMI$BMI)
currentBMI$HEIGHT_CM = as.numeric(currentBMI$HEIGHT_CM)
currentBMI$WEIGHT_KG = as.numeric(currentBMI$WEIGHT_KG)
currentBMI = left_join(currentBMI, ids, by = "Patient_ID")

currentBMISorted = arrange(currentBMI, OBSERVATION_DATE)
currentBMIIndexed = currentBMISorted %>%
  group_by(DIST_ID) %>%
  mutate(BMI_N = row_number())

# Calc BMI by height and weight incase BMI field is NA
# BMI = kg/m^2
currentBMIIndexed$CALC_BMI = currentBMIIndexed$WEIGHT_KG / ((currentBMIIndexed$HEIGHT_CM / 100) ^ 2)
# Grab values present in BMI. If missing, grab values present in CALC_BMI
currentBMIIndexed$BMI_CALC_BMI = ifelse(!is.na(currentBMIIndexed$BMI), currentBMIIndexed$BMI, currentBMIIndexed$CALC_BMI)
# Remove missing data, keep anything with a BMI or that was able to have a BMI calculated
currentBMIMissing = currentBMIIndexed %>%
  filter(is.na(BMI_CALC_BMI))
currentBMINotMissing = currentBMIIndexed %>%
  filter(!is.na(BMI_CALC_BMI))

currentBMINotMissing = select(currentBMINotMissing, DIST_ID, BMI, HEIGHT_CM, WEIGHT_KG, CALC_BMI, BMI_CALC_BMI, BMI_N)
# Remove intermediate dataframes
rm(currentBMI, currentBMIIndexed, currentBMIMissing, currentBMISorted)


################################################################################
#get bp file
################################################################################
# Reading in BP file pulled straight from the database
currentBP = read_csv("Y:/IRB_00096551/UDDB_Table_Files_For_Plotting_2021/uddb_bp.csv")
currentBP = left_join(currentBP, ids, by = "Patient_ID")

currentBPSorted = arrange(currentBP, OBSERVATION_DATE)
currentBPIndexed = currentBPSorted %>%
  group_by(DIST_ID) %>%
  mutate(BP_N = row_number())

currentBPMissing <- currentBPIndexed %>%
  filter(is.na(BP_SYS) | is.na(BP_DIAS))
currentBPNotMissing <- currentBPIndexed %>%
  filter(!is.na(BP_SYS) &!is.na(BP_DIAS))

currentBPNotMissing = select(currentBPNotMissing, DIST_ID, BP_SYS, BP_DIAS, BP_N)
colnames(currentBPNotMissing)[2:3] = c("BP_SYST", "BP_DIAST")
rm(currentBP, currentBPIndexed, currentBPMissing, currentBPSorted)


################################################################################
#get master file
################################################################################
# Need the following variables below and uncorrected_slope
#
#eGFR_outpatient_data <- data.frame(id = as.integer(master.file$dist_id),
#                                   sex = as.character(master.file$sex),
#                                   ethnicity = as.character(master.file$ethn),
#                                   first_age = master.file$age_first_egfr,
#                                   first_egfr = master.file$first_egfr,
#                                   first_ckd_stage = as.character(master.file$first_ckd_stage),
#                                   last_age = master.file$age_last_egfr,
#                                   last_egfr = master.file$last_egfr,
#                                   last_ckd_stage = as.character(master.file$last_ckd_stage),
#                                   number_egfrs = master.file$n_egfr,
#                                   followup_years = master.file$followup_yrs,
#                                   slope = master.file$corrected_slope,
#                                   r2 = master.file$corrected_r2,
#                                   obsv_esrd_followup = master.file$obsv_esrd_fu_yrs,
#                                   esrd_followup = master.file$icd_esrd_fu_yrs,
#                                   dialysis_followup = master.file$icd_dialysis_fu_yrs,
#                                   kidney_followup = master.file$icd_kidney_tx_fu_yrs)

alleGFRs = read_csv("Y:/IRB_00096551/UDDB_Table_Files_For_Plotting_2021/uddb_egfrs_2021.csv")
alleGFRs = left_join(alleGFRs, ids, by = "Patient_ID")

patients = read_csv("Y:/IRB_00096551/UDDB_Table_Files_For_Plotting_2021/Patients_IHC_UUHSC_with_Birth_DT.csv")

patientEGFRs = left_join(alleGFRs, patients, by = "Patient_ID", all.x = TRUE)
patientEGFRs = group_by(patientEGFRs, Patient_ID) %>%
  mutate(first_egfr = min(RESULT_DTM))

alleGFRs = select(alleGFRs, DIST_ID, SEX_CD, ETHNICITY_DSC, )


################################################################################
#get outpatient.egfr.3om
################################################################################
