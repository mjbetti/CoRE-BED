library("data.table")
library("dplyr")

gwas_array <- c("PASS_ADHD_Demontis2018", "PASS_AgeFirstBirth", "PASS_Alzheimers_Jansen2019", "PASS_Anorexia", "PASS_AtrialFibrillation_Nielsen2018", "PASS_Autism", "PASS_BDSCZ_Ruderfer2018", "PASS_BMI1", "PASS_CD_deLange2017", "PASS_Celiac", "PASS_CigarettesPerDay_Liu2019", "PASS_Coronary_Artery_Disease", "PASS_DrinksPerWeek_Liu2019", "PASS_Ever_Smoked", "PASS_GeneralRiskTolerance_KarlssonLinner2019", "PASS_HDL", "PASS_Height1", "PASS_IBD_deLange2017", "PASS_Insomnia_Jansen2019", "PASS_Intelligence_SavageJansen2018", "PASS_IschemicStroke_Malik2018", "PASS_LDL", "PASS_Lupus", "PASS_MDD_Wray2018", "PASS_MedicationUse_Wu2019", "PASS_NumberChildrenEverBorn", "PASS_ProstateCancer", "PASS_ReactionTime_Davies2018", "PASS_Rheumatoid_Arthritis", "PASS_SCZvsBD_Ruderfer2018", "PASS_SleepDuration_Dashti2019", "PASS_Type_2_Diabetes", "PASS_UC_deLange2017", "PASS_Years_of_Education1", "UKB_460K.biochemistry_AlkalinePhosphatase", "UKB_460K.biochemistry_AspartateAminotransferase", "UKB_460K.biochemistry_Cholesterol", "UKB_460K.biochemistry_Creatinine", "UKB_460K.biochemistry_IGF1", "UKB_460K.biochemistry_Phosphate", "UKB_460K.biochemistry_Testosterone_Male", "UKB_460K.biochemistry_TotalBilirubin", "UKB_460K.biochemistry_TotalProtein", "UKB_460K.biochemistry_VitaminD", "UKB_460K.blood_EOSINOPHIL_COUNT", "UKB_460K.blood_PLATELET_COUNT", "UKB_460K.blood_RBC_DISTRIB_WIDTH", "UKB_460K.blood_RED_COUNT", "UKB_460K.blood_WHITE_COUNT", "UKB_460K.bmd_HEEL_TSCOREz", "UKB_460K.body_BALDING1", "UKB_460K.body_BMIz", "UKB_460K.body_HEIGHTz", "UKB_460K.body_WHRadjBMIz", "UKB_460K.bp_DIASTOLICadjMEDz", "UKB_460K.cancer_BREAST", "UKB_460K.cancer_PROSTATE", "UKB_460K.cov_EDU_YEARS", "UKB_460K.disease_AID_SURE", "UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED", "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP", "UKB_460K.lung_FEV1FVCzSMOKE", "UKB_460K.lung_FVCzSMOKE", "UKB_460K.mental_NEUROTICISM", "UKB_460K.other_MORNINGPERSON", "UKB_460K.pigment_SUNBURN", "UKB_460K.repro_MENARCHE_AGE", "UKB_460K.repro_MENOPAUSE_AGE", "UKB_460K.repro_NumberChildrenEverBorn_Pooled")

#Declare the path of the working directory
in_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G"

cat_df <- data.frame()

for (gwas in gwas_array) {
	print(gwas)
	total_file <- fread(paste0(gwas, "_h2_total.txt"), header = FALSE, sep = " ", quote = "")
	total_df <- as.data.frame(total_file)
	total_df <- total_df[,1]
	
	results_file <- fread(paste0(gwas, "_h2.results"), header = TRUE, sep = "\t", quote = "")
	results_df <- as.data.frame(results_file)
	results_df <- data.frame(results_df[,3])
	results_df <- transpose(results_df)
	results_df <- results_df * total_df
	cat_df <- rbind(cat_df, results_df)
}

cat_df <- data.frame(gwas_array, cat_df)
names(cat_df) <- c("GWAS", "Class1", "Class2", "Class3", "Class4", "Class5", "Class6", "Class7", "Class8", "any_corebed")

write.table(cat_df, file = "core_bed_h2_absolute_values.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)