import os, sys, requests, argparse, pybedtools, pandas as pd
from itertools import chain
import subprocess

ref_dir = "/home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/playground/ref_files_lifted_hg38/test"

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-t", "--tissue", type = str, required = True, help = "the input tissue type")
args = parser.parse_args()

#tissues = ["adipose_adipocyte", "adipose_34m", "adipose_omental_fat_pad_51f", "adipose_omental_fat_pad_53f", "adipose_subcutaneous_25f", "adipose_subcutaneous_41f", "adipose_subcutaneous_49f", "adipose_subcutaneous_59f", "adipose_subcutaneous_81f", "blood_cd4_alpha_beta_memory_t_cell", "blood_cd4_alpha_beta_memory_t_cell_blood_origin", "blood_cd4_alpha_beta_t_cell", "blood_cd4_alpha_beta_t_cell_33f", "blood_cd4_alpha_beta_t_cell_21m", "blood_cd4_alpha_beta_t_cell_37m", "blood_cd4_alpha_beta_t_cell_phorbol", "blood_cd4_cd25_alpha_beta_t_cell", "blood_cd8_alpha_beta_memory_t_cell", "blood_cd8_alpha_beta_t_cell", "blood_cd8_alpha_beta_t_cell_33f", "blood_cd8_alpha_beta_t_cell_34f", "blood_cd8_alpha_beta_t_cell_21m", "blood_cd8_alpha_beta_t_cell_28m", "blood_cd8_alpha_beta_t_cell_37m", "blood_effector_memory_cd4_alpha_beta_t_cell", "blood_mononuclear_cell_male", "blood_naive_thymus_cd4_alpha_beta_t_cell", "blood_naive_thymus_cd4_alpha_beta_t_cell_35f", "blood_naive_thymus_cd4_alpha_beta_t_cell_26m", "blood_peripheral_mononuclear_cell_28f", "blood_peripheral_mononuclear_cell_27m", "blood_peripheral_mononuclear_cell_28m", "blood_peripheral_mononuclear_cell_32m", "blood_peripheral_mononuclear_cell_39m", "blood_regulatory_t_cell_35f", "blood_regulatory_t_cell_28m", "blood_regulatory_t_cell_blood_origin", "blood_t_cell", "blood_t_cell_21m", "blood_t_cell_36m", "blood_t_cell_37m", "blood_t1_helper_cell", "blood_t1_helper_cell_26f", "blood_t1_helper_cell_33m", "blood_t17_helper_cell", "blood_t17_helper_cell_blood_origin", "blood_t17_helper_cell_phorbol", "blood_t2_helper_cell_26f", "blood_t2_helper_cell_33m", "bone_arm", "bone_femur", "bone_marrow_stroma", "bone_leg", "bone_osteoblast", "brain_ammons_horn_84m", "brain_angular_gyrus_75f", "brain_angular_gyrus_81m", "brain_astrocyte", "brain_astrocyte_cerebellum", "brain_astrocyte_hippocampus", "brain_astrocyte_spinal_cord", "brain_embryo_112_days", "brain_embryo_56_58_days", "brain_embryo_80_days", "brain_embryo_105_days_f", "brain_embryo_109_days_f", "brain_embryo_117_days_f", "brain_embryo_120_days_f", "brain_embryo_142_days_f", "brain_embryo_17_weeks_f", "brain_embryo_85_days_f", "brain_embryo_96_days_f", "brain_embryo_101_days_m", "brain_embryo_104_days_m", "brain_embryo_105_days_m", "brain_embryo_122_days_m", "brain_embryo_72_76_days_m", "brain_caudate_nucleus_75f", "brain_caudate_nucleus_78m", "brain_caudate_nucleus_81m", "brain_cerebellar_cortex_78_81m", "brain_cerebellum_27_35m", "brain_cerebellum_53m", "brain_cingulate_gyrus_75f", "brain_cingulate_gyrus_81m", "brain_frontal_cortex_67_80f", "brain_frontal_cortex_27_35m", "brain_germinal_matrix_20_weeks_m", "brain_globus_pallidus_78_84m", "brain_inferior_parietal_cortex_84m", "brain_hippocampus_75f", "brain_hippocampus_73m", "brain_medulla_oblongata_78_84m", "brain_midbrain_78_84m", "brain_middle_frontal_area_75f", "brain_middle_frontal_area_81m", "brain_middle_frontal_gyrus_78m", "brain_occipital_lobe_84m", "brain_pons_78m", "brain_putamen_78m", "brain_substantia_nigra_75f", "brain_substantia_nigra_81m", "brain_superior_temporal_gyrus_84m", "brain_temporal_lobe_75f", "brain_temporal_lobe_81m", "cancer_prostate_epithelial_22rv1", "cancer_pancreas_adenocarcinoma_8988t", "cancer_glioblastoma_a172", "cancer_lung_epithelial_carcinoma_a549", "cancer_lung_epithelial_carcinoma_a549_treated_ethanol_1hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_1hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_12hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_15min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_2hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_20min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_25min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_3hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_30min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_4hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_6hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_7hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_8hr", "cancer_muscle_ewing_sarcoma_a673", "cancer_adenoid_cystic_carcinoma_acc112", "cancer_renal_cell_carcinoma_achn", "cancer_neuroblastoma_be2c", "cancer_prostate_c4-2b", "cancer_colorectal_adenocarcinoma_caco-2", "cancer_kidney_clear_cell_carcinoma_caki2", "cancer_myelogenous_leukemia_cmk", "cancer_melanoma_colo829", "cancer_desmoplastic_medulloblastoma_daoy", "cancer_acute_lymphoblastic_leukemia_dnd-41", "cancer_b_cell_lymphoma_dohh2", "cancer_kidney_rhaboid_tumor_g401", "cancer_neuroglioma_h4", "cancer_glioblastoma_h54", "cancer_haploid_myelogenous_leukemia_hap-1", "cancer_colorectal_adenocarcinoma_hct116", "cancer_cervix_adenocarcinoma_hela-s3", "cancer_cervix_adenocarcinoma_hela-s3_g1b_phase", "cancer_cervix_adenocarcinoma_hela-s3_treated_interferon_alpha_4hr", "cancer_hepatocellular_carcinoma_hepg2", "cancer_acute_promyelocytic_leukemia_hl-60", "cancer_colorectal_adenocarcinoma_ht-29", "cancer_fibrosarcoma_ht1080", "cancer_hepatocellular_carcinoma_huh-7.5", "cancer_endometrial_adenocarcinoma_ishikawa_treated_dmso_1hr", "cancer_endometrial_adenocarcinoma_ishikawa_treated_17b-estradiol_30min", "cancer_endometrial_adenocarcinoma_ishikawa_treated_afimoxifene_30min", "cancer_myelogenous_leukemia_k562", "cancer_myelogenous_leukemia_k562_g1_phase", "cancer_myelogenous_leukemia_k562_g2_phase", "cancer_myelogenous_leukemia_k562_treated_dmso_72hr", "cancer_myelogenous_leukemia_k562_treated_vorinostat_72hr", "cancer_myelogenous_leukemia_k562_treated_sodium_butyrate_72hr", "cancer_b_cell_lymphoma_karpas-422", "cancer_myelogenous_leukemia_kbm-7", "cancer_myeloma_kms-11", "cancer_acute_lymphoblastic_leukemia_kopt-k1", "cancer_prostate_adenocarcinoma_lncap_clone_fgc", "cancer_prostate_adenocarcinoma_lncap_clone_fgc_treated_17b_12hr", "cancer_acute_lymphoblastic_leukemia_loucy", "cancer_colorectal_adenocarcinoma_lovo", "cancer_glioblastoma_m059j", "cancer_mammary_gland_adenocarcinoma_mcf-7", "cancer_mammary_gland_adenocarcinoma_mcf-7_originated_from_mcf-7", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_lactate_24hr", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_estradiol_1hr", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_ctcf_shrna_knockown", "cancer_medulloblastoma", "cancer_osteosarcoma_mg63", "cancer_burkitt_lymphoma_namalwa", "cancer_burkitt_lymphoma_namalwa_treated_sendai_virus_2hr", "cancer_acute_promyelocytic_leukemia_nb4", "cancer_squamous_cell_carcinoma_nci-h226", "cancer_large_cell_lung_nci-h460", "cancer_myeloma_nci-h929", "cancer_testicular_embryonal_carcinoma_nt2_d1", "cancer_b_cell_lymphoma_oci-ly1", "cancer_b_cell_lymphoma_oci-ly3", "cancer_b_cell_lymphoma_oci-ly7", "cancer_pancreas_duct_epithelial_carcinoma_panc1", "cancer_parathyroid_adenoma_62m", "cancer_prostate_adenocarcinoma_pc-3", "cancer_lung_adenocarcinoma_pc-9", "cancer_renal_cell_adenocarcinoma_rcc_7860", "cancer_renal_cell_carcinoma", "cancer_colon_carcinoma_rko", "cancer_melanoma_rpmi-7951", "cancer_plasma_cell_myeloma_rpmi8226", "cancer_rhabdomyosarcoma_sjcrh30", "cancer_osteosarcoma_sjsa1", "cancer_melanoma_sk-mel-5", "cancer_neuroblastoma_sk-n-dz", "cancer_neuroblastoma_sk-n-dz_treated_dmso_72hr", "cancer_neuroepithelioma_sk-n-mc", "cancer_neuroblastoma_sk-n-sh", "cancer_neuroblastoma_sk-n-sh_treated_retinoic_acid_48hr", "cancer_b_cell_lymphoma_su-dhl-6", "cancer_colorectal_adenocarcinoma_sw480", "cancer_mammary_gland_ductal_carcinoma_t47d", "cancer_mammary_gland_ductal_carcinoma_t47d_treated_17b-estradiol_30min", "cancer_prostate_epithelial_carcinoma_vcap", "cancer_eye_retinoblastoma_weri-rb-1", "digestive_colon_mucosa_56f", "digestive_colon_mucosa_73f", "digestive_duodenum_mucosa_59m", "digestive_duodenum_mucosa_76m", "digestive_esophagus_30f", "digestive_esophagus_34m", "digestive_esophagus_muscularis_mucosa_53f", "digestive_esophagus_muscularis_mucosa_37m", "digestive_esophagus_squamous_epithelium_53f", "digestive_esophagus_squamous_epithelium_37m", "digestive_gastroesophageal_sphincter_53f", "digestive_large_intestine_embryo_103_days_f", "digestive_large_intestine_embryo_105_days_f", "digestive_large_intestine_embryo_107_days_f", "digestive_large_intestine_embryo_108_days_f", "digestive_large_intestine_embryo_110_days_f", "digestive_large_intestine_embryo_120_days_f", "digestive_large_intestine_embryo_91_days_f", "digestive_large_intestine_embryo_98_days_f", "digestive_large_intestine_embryo_105_days_m", "digestive_large_intestine_embryo_108_days_m", "digestive_large_intestine_embryo_113_days_m", "digestive_large_intestine_embryo_115_days_m", "digestive_large_intestine_embryo_91_days_m", "digestive_rectum_mucosa_50f", "digestive_rectum_mucosa_61f", "digestive_rectum_mucosa_59m", "digestive_stomach_mucosa_59m", "digestive_peyers_patch_53f", "digestive_peyers_patch_37m", "digestive_peyers_patch_54m", "digestive_sigmoid_colon_51f", "digestive_sigmoid_colon_53f", "digestive_sigmoid_colon_34m", "digestive_sigmoid_colon_54m", "digestive_sigmoid_colon_3m", "digestive_small_intestine_30f", "digestive_small_intestine_105_days_f", "digestive_small_intestine_107_days_f", "digestive_small_intestine_108_days_f", "digestive_small_intestine_110_days_f", "digestive_small_intestine_120_days_f", "digestive_small_intestine_91_days_f", "digestive_small_intestine_98_days_f", "digestive_small_intestine_34m", "digestive_small_intestine_3m", "digestive_small_intestine_105_days_m", "digestive_small_intestine_108_days_m", "digestive_small_intestine_115_days_m", "digestive_small_intestine_87_days_m", "digestive_small_intestine_91_days_m", "digestive_stomach_101_days", "digestive_stomach_30f", "digestive_stomach_51f", "digestive_stomach_53f", "digestive_stomach_f", "digestive_stomach_105_days_f", "digestive_stomach_107_days_f", "digestive_stomach_108_days_f", "digestive_stomach_121_days_f", "digestive_stomach_147_days_f", "digestive_stomach_96_days_f", "digestive_stomach_98_days_f", "digestive_stomach_34m", "digestive_stomach_54m", "digestive_stomach_3m", "digestive_stomach_108_days_m", "digestive_stomach_127_days_m", "digestive_stomach_58_76_days_m", "digestive_stomach_91_days_m", "digestive_transverse_colon_51f", "digestive_transverse_colon_53f", "digestive_transverse_colon_37m", "digestive_transverse_colon_54m", "endocrine_adrenal_gland_96_days", "endocrine_adrenal_gland_30f", "endocrine_adrenal_gland_51f", "endocrine_adrenal_gland_53f", "endocrine_adrenal_gland_108_days_f", "endocrine_adrenal_gland_113_days_f", "endocrine_adrenal_gland_85_days_f", "endocrine_adrenal_gland_34m", "endocrine_adrenal_gland_37m", "endocrine_adrenal_gland_54m", "endocrine_adrenal_gland_101_days_m", "endocrine_adrenal_gland_108_days_m", "endocrine_adrenal_gland_85_days_m", "endocrine_adrenal_gland_97_days_m", "endocrine_pancreas_59", "endocrine_pancreas_m", "endocrine_pancreas_45m", "endocrine_pancreas_46m", "endocrine_ovary_30f", "endocrine_ovary_51f", "endocrine_ovary_53f", "endocrine_ovary_embryo_f", "endocrine_testis_37m", "endocrine_testis_54m", "endocrine_testis_embryo_m", "endocrine_thyroid_gland_51f", "endocrine_thyroid_gland_53f", "endocrine_thyroid_gland_37m", "endocrine_thyroid_gland_54m", "endothelial_brain_microvascular", "endothelial_dermis_blood_vessel_adult_f", "endothelial_dermis_blood_vessel_newborn_m", "endothelial_dermis_microvascular_lymphatoc_vessel_f", "endothelial_dermis_microvascular_lymphatoc_vessel_m", "endothelial_umbilical_vein_newborn_m", "endothelial_umbilical_vein_newborn", "endothelial_glomerulus", "endothelial_kidney_capillary_113_days_f", "endothelial_lung_microvascular_f", "endothelial_pulmonary_artery_f", "epithelial_amnion", "epithelial_bronchial", "epithelial_bronchial_f_treated_retinoic_acid", "epithelial_choroid_plexus", "epithelial_colon", "epithelial_esophagus", "epithelial_prostate", "epithelial_prostate_m", "epithelial_proximal_tubule", "epithelial_foreskin_keratinocyte_newborn_m", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_2.5_days", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_5.5_days", "epithelial_glomerulus_visceral_3", "epithelial_proximal_tubule_hk-2", "epithelial_pancreatic_duct_dpde6-e6e7", "epithelial_bone_marrow_hs-27a", "epithelial_iris_pigment", "epithelial_keratinocyte_f", "epithelial_keratinocyte_m", "epithelial_kidney", "epithelial_glomerulus_43_62_m", "epithelial_tubule_80f_62m", "epithelial_tubule_80f_treated_cisplatin", "epithelial_skin_leg_53f", "epithelial_skin_leg_37m", "epithelial_mammary_luminal_33f", "epithelial_mammary_f", "epithelial_mammary_18f", "epithelial_mammary_50f", "epithelial_breast_mcf_10a", "epithelial_breast_mcf_10a_treated_tamoxifen_24hr", "epithelial_breast_mcf_10a_treated_tamoxifen_6hr", "epithelial_mammary_myoepithelial_33f", "epithelial_mammary_myoepithelial_36f", "epithelial_non-pigmented_ciliary", "epithelial_renal_cortical", "epithelial_retinal", "epithelial_prostate_rwpe1", "epithelial_prostate_rwpe2", "epithelial_skin_of_body_82_days_f", "esc_deriv_bipolar_neuron_gm23338_origin_treated_doxycycline_4_days", "esc_deriv_cardiac_mesoderm_h7-hesc_origin", "esc_deriv_cardiac_muscle_rues2_origin", "esc_deriv_neural_progenitor_h9_origin_1", "esc_deriv_ectodermal", "esc_deriv_endodermal", "esc_deriv_endodermal_hues64_origin", "esc_deriv_hepatocyte_h9_origin", "esc_deriv_mesenchymal_stem_h1-hesc_origin", "esc_deriv_mesendoderm_h1-hesc_origin", "esc_deriv_mesodermal_hues64_origin", "esc_deriv_neural_crest_h1-hesc_origin", "esc_deriv_neural_progenitor_h9_origin_2", "esc_deriv_neural_progenitor_h1-hesc_origin", "esc_deriv_neural_progenitor_h9_origin_3", "esc_deriv_neuron_h9_origin", "esc_deriv_smooth_muscle_h9_origin", "esc_deriv_trophoblast_h1-hesc_origin", "esc_elf-1", "esc_es-i3", "esc_h1-hesc", "esc_h7-hesc", "esc_h9", "esc_hues48", "esc_hues6", "esc_hues64", "esc_ucsf-4", "eye_56_days_76_days_m", "eye_76_days_f", "eye_125_days_103_days_m", "eye_74_days_85_days", "eye_89_days_f", "heart_aorta_30f", "heart_aorta_34m", "heart_ascending_aorta_51f", "heart_ascending_aorta_53f", "heart_coronary_artery_51f", "heart_coronary_artery_53f", "heart_101_days", "heart_59_days_76_days_f", "heart_80_days", "heart_96_days", "heart_103_days_f", "heart_105_days_f", "heart_110_days_f", "heart_116_days_98_days_f", "heart_116_days_117_days_f", "heart_147_days_f", "heart_91_days_f", "heart_left_ventricle_53f", "heart_left_ventricle_101_days_103_days_f", "heart_left_ventricle_136_days_f", "heart_left_ventricle_34m", "heart_left_ventricle_3m", "heart_27m_35m", "heart_3m", "heart_105_days_m", "heart_110_days_m", "heart_120_days_m", "heart_72_days_76_days_m", "heart_91_days_m", "heart_96_days_m", "heart_right_ventricle_101_days_103_days_f", "heart_right_ventricle_34m", "heart_right_ventricle_3m", "heart_left_atrium_101_days_f", "heart_right_atrium_51f", "heart_right_atrium_53f", "heart_right_atrium_34m", "heart_thoracic_aorta_37m", "heart_thoracic_aorta_54m", "heart_tibial_artery_53f", "heart_tibial_artery_37m", "hsc_and_b_cell_b_cell", "hsc_and_b_cell_b_cell_27f", "hsc_and_b_cell_b_cell_27f_43f", "hsc_and_b_cell_b_cell_34f", "hsc_and_b_cell_b_cell_43f", "hsc_and_b_cell_b_cell_21m", "hsc_and_b_cell_b_cell_37m", "hsc_and_b_cell_cd14_monocyte_f", "hsc_and_b_cell_cd14_monocyte_34f", "hsc_and_b_cell_cd14_monocyte_21m", "hsc_and_b_cell_cd14_monocyte_37m", "hsc_and_b_cell_cd1c_myeloid_dendritic", "hsc_and_b_cell_cmp_cd34", "hsc_and_b_cell_cmp_cd34_f", "hsc_and_b_cell_cmp_cd34_27f", "hsc_and_b_cell_cmp_cd34_33f", "hsc_and_b_cell_cmp_cd34_50f", "hsc_and_b_cell_cmp_cd34_m", "hsc_and_b_cell_cmp_cd34_adult_m", "hsc_and_b_cell_cmp_cd34_23m", "hsc_and_b_cell_cmp_cd34_36m", "hsc_and_b_cell_cmp_cd34_37m", "hsc_and_b_cell_cmp_cd34_42m", "hsc_and_b_cell_cmp_cd34_43m", "hsc_and_b_cell_cmp_cd34_49m", "hsc_and_b_cell_germinal_center", "hsc_and_b_cell_mpp", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_20_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_11_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_13_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_15_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_17_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_18_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_4_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_6_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_8_days", "hsc_and_b_cell_lymphocyte_jurkat_clone_e61", "hsc_and_b_cell_lymphocyte_naive_b_cell", "hsc_and_b_cell_nk_34f", "hsc_and_b_cell_nk_21m", "hsc_and_b_cell_nk_37m", "hsc_and_b_cell_neutrophil", "hsc_and_b_cell_neutrophil_m", "ipsc_cwru1_m", "ipsc_gm23338_53m_gm23248_origin", "ipsc_ips_df_19.11_newborn_m", "ipsc_ips_df_19.7_newborn_m", "ipsc_ips_df_14.7_newborn_m", "ipsc_ips_df_6.9_newborn_m", "ipsc_ips-11a_36m", "ipsc_ips-15b_48f", "ipsc_ips-18a_48f", "ipsc_ips-18c_48f", "ipsc_ips-20b_55m", "ipsc_ips-nihi11_71m_ag20443_origin", "ipsc_ips-nihi7_85f_ag08395_origin", "ipsc_l1-s8", "ipsc_l1-s8r", "kidney_hek293", "kidney_hek293t", "kidney_59_days_59_days_f", "kidney_80_days", "kidney_105_days_f", "kidney_108_days_f", "kidney_113_days_f", "kidney_120_days_f", "kidney_121_days_f", "kidney_76_days_f_76_days_m", "kidney_85_days_f", "kidney_50m", "kidney_67m", "kidney_105_days_m", "kidney_85_days_m", "kidney_87_days_m", "kidney_left_107_days_f", "kidney_left_110_days_f", "kidney_left_147_days_f", "kidney_left_59_days_f_91_days_m", "kidney_left_87_days_f", "kidney_left_89_days_f", "kidney_left_115_days_m", "kidney_left_87_days_m", "kidney_left_96_days_m", "kidney_left_renal_cortex_interstitium_105_days_m", "kidney_left_renal_cortex_interstitium_120_days_m", "kidney_left_renal_pelvis_105_days_m", "kidney_left_renal_pelvis_120_days_m", "kidney_renal_cortex_interstitium_103_days_f", "kidney_renal_cortex_interstitium_120_days_f", "kidney_renal_cortex_interstitium_89_days_f", "kidney_renal_cortex_interstitium_96_days_f", "kidney_renal_cortex_interstitium_108_days_m", "kidney_renal_cortex_interstitium_113_days_m", "kidney_renal_cortex_interstitium_127_days_m", "kidney_renal_cortex_interstitium_91_days_m", "kidney_renal_cortex_interstitium_97_days_m", "kidney_renal_pelvis_103_days_f", "kidney_renal_pelvis_105_days_f", "kidney_renal_pelvis_89_days_f", "kidney_renal_pelvis_96_days_f", "kidney_renal_pelvis_108_days_m", "kidney_renal_pelvis_113_days_m", "kidney_renal_pelvis_127_days_m", "kidney_renal_pelvis_91_days_m", "kidney_renal_pelvis_97_days_m", "kidney_right_107_days_f", "kidney_right_117_days_f", "kidney_right_147_days_f", "kidney_right_87_days_f", "kidney_right_98_days_f", "kidney_right_108_days_m", "kidney_right_115_days_m", "kidney_right_87_days_m", "kidney_right_91_days_m", "kidney_right_96_days_m", "kidney_right_renal_cortex_interstitium_105_days_m", "kidney_right_renal_cortex_interstitium_120_days_m", "kidney_right_renal_pelvis_105_days_m", "kidney_right_renal_pelvis_120_days_m", "liver_hepatic_stellate_cell_59f", "liver_hepatocyte", "liver_59_days_80_days", "liver_25f", "liver_101_days_113_days_f", "liver_31m", "liver_32m", "liver_78m", "liver_right_lobe_53f", "lung_left_105_days_f", "lung_left_107_days_f", "lung_left_108_days_f", "lung_left_110_days_f", "lung_left_117_days_f", "lung_left_91_days_f", "lung_left_98_days_f", "lung_left_105_days_m", "lung_left_113_days_m", "lung_left_115_days_m", "lung_left_87_days_m", "lung_left_91_days_m", "lung_left_96_days_m", "lung_101_days", "lung_112_days", "lung_67_days", "lung_80_days_76_days_m", "lung_30f", "lung_108_days_f", "lung_120_days_f", "lung_76_days_f", "lung_82_days_f", "lung_85_days_f", "lung_96_days_f", "lung_3m", "lung_103_days_m", "lung_108_days_m", "lung_54_days_58_days_m", "lung_82_days_m", "lung_right_105_days_f", "lung_right_107_days_f", "lung_right_108_days_f", "lung_right_110_days_f", "lung_right_117_days_f", "lung_right_91_days_f", "lung_right_98_days_f", "lung_right_105_days_m", "lung_right_115_days_m", "lung_right_87_days_m", "lung_right_96_days_m", "lung_left_upper_lobe_51f", "lung_left_upper_lobe_53f", "lung_left_upper_lobe_37m", "lymphoblastoid_gm06990", "lymphoblastoid_gm08714", "lymphoblastoid_gm10248", "lymphoblastoid_gm10266", "lymphoblastoid_gm12864", "lymphoblastoid_gm12865", "lymphoblastoid_gm12875", "lymphoblastoid_gm12878", "lymphoblastoid_gm12891", "lymphoblastoid_gm12892", "lymphoblastoid_gm18507", "lymphoblastoid_gm19238", "lymphoblastoid_gm19239", "lymphoblastoid_gm19240", "mesench_adipocyte_msc_origin", "mesench_msc_dedifferentiated_amniotic_fluid_origin", "mesench_embryonic_facial_prominence_53_days_58_days", "mesench_msc_adipose_origin", "muscle_cardiac_myocyte", "muscle_forelimb_108_days_f", "muscle_gastrocnemius_medialis_51f", "muscle_gastrocnemius_medialis_53f", "muscle_gastrocnemius_medialis_37m", "muscle_gastrocnemius_medialis_54m", "muscle_leg_hindlimb_120_days_m", "muscle_arm_101_days", "muscle_arm_105_days_f", "muscle_arm_115_days_f", "muscle_arm_120_days_f", "muscle_arm_85_days_f", "muscle_arm_98_days_f", "muscle_arm_101_days_m", "muscle_arm_104_days_m", "muscle_arm_105_days_m", "muscle_arm_113_days_m", "muscle_arm_115_days_m", "muscle_arm_120_days_m", "muscle_arm_127_days_m", "muscle_arm_96_days_m", "muscle_arm_97_days_m", "muscle_back_105_days_f", "muscle_back_113_days_f", "muscle_back_115_days_f", "muscle_back_85_days_f", "muscle_back_98_days_f", "muscle_back_101_days_m", "muscle_back_104_days_m", "muscle_back_105_days_m", "muscle_back_108_days_m", "muscle_back_127_days_m", "muscle_back_91_days_m", "muscle_back_96_days_m", "muscle_back_97_days_m", "muscle_leg_105_days_f", "muscle_leg_110_days_f", "muscle_leg_113_days_f", "muscle_leg_115_days_f", "muscle_leg_85_days_f", "muscle_leg_101_days_m", "muscle_leg_104_days_m", "muscle_leg_105_days_m", "muscle_leg_115_days_m", "muscle_leg_127_days_m", "muscle_leg_96_days_m", "muscle_leg_97_days_m", "muscle_trunk_113_days_f", "muscle_trunk_115_days_f", "muscle_trunk_120_days_f", "muscle_trunk_121_days_f", "muscle_psoas_30f", "muscle_psoas_27m_35m", "muscle_psoas_34m", "muscle_psoas_3m", "muscle_skeletal_cell", "muscle_skeletal_tissue", "muscle_skeletal_tissue_72f", "muscle_skeletal_tissue_54m", "muscle_tongue_59_days_f_76_days_f", "muscle_tongue_72_days_m", "myosat_skeletal_muscle_myoblast_lhcn-m2", "myosat_myocyte_lhcn-m2_origin", "myosat_myotube_skeletal_muscle_myoblast_origin", "myosat_skeletal_muscle_myoblast", "myosat_skeletal_muscle_myoblast_22m", "myosat_skeletal_muscle_satellite_cell_mesoderm_origin_f", "neurosph_15_weeks_ganglionic_eminence_origin", "neurosph_17_weeks_f_ganglionic_cortex_origin", "neurosph_17_weeks_f_ganglionic_eminence_origin", "neurosph_olfactory_cell_line", "other_breast_epithelium_51f", "other_breast_epithelium_53f", "other_epidermal_melanocyte", "other_foreskin_melanocyte_newborn_m", "other_limb_embryo_53_days_56_days", "other_limb_embryo_58_days_59_days", "other_mammary_stem_cell", "pancreas_body_51f", "pancreas_body_53f", "pancreas_body_37m", "pancreas_body_54m", "pancreas_islet_precursor_cell", "pancreas_30f", "pancreas_34m", "placenta_and_eem_amnion_16_weeks_m", "placenta_and_eem_amnion_stem_cell", "placenta_and_eem_chorion", "placenta_and_eem_chorion_40_weeks_f", "placenta_and_eem_chorion_16_weeks_m", "placenta_and_eem_chorionic_villus_16_weeks", "placenta_and_eem_chorionic_villus_40_weeks_f", "placenta_and_eem_chorionic_villus_16_weeks_m", "placenta_and_eem_chorionic_villus_38_weeks_m", "placenta_and_eem_trophoblast_htr-8_svneo", "placenta_and_eem_placenta_102_days", "placenta_and_eem_placenta_16_weeks", "placenta_and_eem_placenta_53_days", "placenta_and_eem_placenta_56_days_59_days", "placenta_and_eem_placenta_101_days_f_105_days_m", "placenta_and_eem_placenta_105_days_f", "placenta_and_eem_placenta_108_days_f", "placenta_and_eem_placenta_113_days_f", "placenta_and_eem_placenta_85_days_f", "placenta_and_eem_placenta_16_weeks_m", "placenta_and_eem_placenta_85_days_m", "placenta_and_eem_placenta_91_days_m", "placenta_and_eem_placental_basal_plate_40_weeks_f", "placenta_and_eem_placental_basal_plate_38_weeks_m", "placenta_and_eem_trophoblast_cell_17_weeks_18_weeks", "placenta_and_eem_trophoblast_cell_21_weeks", "placenta_and_eem_trophoblast_cell_23_weeks", "placenta_and_eem_trophoblast_cell_39_weeks_40_weeks", "placenta_and_eem_trophoblast_20_weeks_f", "placenta_and_eem_trophoblast_40_weeks_f", "placenta_and_eem_umbilical_cord_59_days_76_days_m", "pns_spinal_cord_108_days_f", "pns_spinal_cord_113_days_f", "pns_spinal_cord_59_days_f_72_days_m", "pns_spinal_cord_87_days_f", "pns_spinal_cord_89_days_f", "pns_spinal_cord_105_days_m", "pns_spinal_cord_96_days_m", "pns_tibial_nerve_51f", "pns_tibial_nerve_53f", "pns_tibial_nerve_37m", "reproductive_prostate_gland_37m", "reproductive_prostate_gland_54m", "reproductive_uterus_53f", "reproductive_vagina_51f", "reproductive_vagina_53f", "sm_muscle_colon_56f", "sm_muscle_colon_77f", "sm_muscle_duodenum_59m", "sm_muscle_duodenum_73m", "sm_muscle_rectum_50f", "sm_muscle_brain_vasculature_smooth_cell", "sm_muscle_stomach_84f", "sm_muscle_stomach_59m", "spleen_112_days", "spleen_30f", "spleen_53f", "spleen_34m", "spleen_54m", "spleen_3m", "stromal_skin_fibroblast_ag04449", "stromal_lung_fibroblast_ag04450", "stromal_skin_fibroblast_ag08395", "stromal_lung_fibroblast_ag08396", "stromal_skin_fibroblast_ag08396", "stromal_gingival_fibroblast_ag09319", "stromal_skin_fibroblast_ag10803", "stromal_skin_fibroblast_ag20443", "stromal_skin_fibroblast_bj", "stromal_brain_pericyte", "stromal_cardiac_fibroblast", "stromal_cardiac_fibroblast_f", "stromal_cardiac_fibroblast_94_days_f_98_days_f", "stromal_skin_fibroblast_eh", "stromal_skin_fibroblast_el", "stromal_skin_fibroblast_elr", "stromal_breast_fibroblast_17f", "stromal_breast_fibroblast_26f", "stromal_dermis_fibroblast", "stromal_dermis_fibroblast_f", "stromal_dermis_fibroblast_none_f", "stromal_gingival_fibroblast", "stromal_lung_fibroblast", "stromal_lung_fibroblast_11f_45m", "stromal_lung_fibroblast_45m", "stromal_mammary_fibroblast_f", "stromal_peridontal_ligament_fibroblast_m", "stromal_pulmonary_artery_fibroblast", "stromal_skin_fibroblast_abdomen_97_days_m", "stromal_aorta_fibroblast_f", "stromal_conjunctiva_fibroblast", "stromal_villous_mesenchyme_fibroblast", "stromal_foreskin_fibroblast_newborn_m", "stromal_skin_fibroblast_gm03348", "stromal_skin_fibroblast_gm04503", "stromal_skin_fibroblast_gm04504", "stromal_skin_fibroblast_gm23248", "stromal_foreskin_fibroblast_hff-myc_foreskin_fibroblast_origin", "stromal_lung_fibroblast_imr-90", "stromal_lung_fibroblast_97_days_m", "stromal_bone_marrow_m", "stromal_lung_fibroblast_wi38", "thymus_embryo_f", "thymus_105_days_f", "thymus_110_days_f", "thymus_113_days_f", "thymus_147_days_f", "thymus_98_days_f", "thymus_3m", "thymus_104_days_m", "thymus_108_days_m", "thymus_113_days_m", "thymus_127_days_m", "urinary_bladder_34m", "urinary_bladder_76_days_m", "urinary_urothelium_cell_line"]

##TSS Overlap
#Import the TSS file, testing various TSS thresholds
bed_tss_path = "/home/bettimj/reference_genomes/refTSS_v3.3_human_coordinate.hg38.bed.gz"

bed_tss_orig_5k1k = pd.read_csv(bed_tss_path, sep = "\t", header = None)
bed_tss_orig_5k1k.iloc[:,1] = bed_tss_orig_5k1k.iloc[:,1] - 5000
bed_tss_orig_5k1k.iloc[:,2] = bed_tss_orig_5k1k.iloc[:,2] + 1000
bed_tss_orig_5k1k = bed_tss_orig_5k1k[bed_tss_orig_5k1k.iloc[:,1] >= 0]
bed_tss_5k1k = pybedtools.BedTool().from_dataframe(bed_tss_orig_5k1k)

bed_tss_orig_25k25k = pd.read_csv(bed_tss_path, sep = "\t", header = None)
bed_tss_orig_25k25k.iloc[:,1] = bed_tss_orig_25k25k.iloc[:,1] - 2500
bed_tss_orig_25k25k.iloc[:,2] = bed_tss_orig_25k25k.iloc[:,2] + 2500
bed_tss_orig_25k25k = bed_tss_orig_25k25k[bed_tss_orig_25k25k.iloc[:,1] >= 0]
bed_tss_25k25k = pybedtools.BedTool().from_dataframe(bed_tss_orig_25k25k)


#Import the binned genome bed
bed_genome_path = "/home/bettimj/reference_genomes/hg38.genome.1kb.bed"
bed_genome = pd.read_csv(bed_genome_path, sep = "\t", header = None)
bed_genome = pybedtools.BedTool().from_dataframe(bed_genome)

#Find overlaps
overlaps_tss_5k1k = bed_genome.intersect(bed_tss_5k1k, u = True)
overlaps_tss_5k1k = overlaps_tss_5k1k.sort()

overlaps_tss_25k25k = bed_genome.intersect(bed_tss_25k25k, u = True)
overlaps_tss_25k25k = overlaps_tss_25k25k.sort()

bed_genome = pybedtools.BedTool.to_dataframe(bed_genome)
r, c = bed_genome.shape
ones = "1"

overlaps_tss_5k1k = pybedtools.BedTool.to_dataframe(overlaps_tss_5k1k)
overlaps_tss_5k1k["tss_5k1k"] = ones
merged_df_5k1k = bed_genome.merge(overlaps_tss_5k1k, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
merged_df_5k1k["tss_5k1k"] = merged_df_5k1k["tss_5k1k"].fillna("0")

overlaps_tss_25k25k = pybedtools.BedTool.to_dataframe(overlaps_tss_25k25k)
overlaps_tss_25k25k["tss_25k25k"] = ones
merged_df_25k25k = merged_df_5k1k.merge(overlaps_tss_25k25k, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
merged_df_25k25k["tss_25k25k"] = merged_df_25k25k["tss_25k25k"].fillna("0")

##Tissue-specific marks overlap
bed_27ac_path = ref_dir + "/" + args.tissue + "_27ac_hg38.bed.gz"
bed_27me3_path = ref_dir + "/" + args.tissue + "_27me3_hg38.bed.gz"
bed_4me1_path = ref_dir + "/" + args.tissue + "_4me1_hg38.bed.gz"
bed_4me3_path = ref_dir + "/" + args.tissue + "_4me3_hg38.bed.gz"
bed_atac_path = ref_dir + "/" + args.tissue + "_ATAC-seq_hg38.bed.gz"
bed_ctcf_path = ref_dir + "/" + args.tissue + "_CTCF_hg38.bed.gz"
bed_dnase_path = ref_dir + "/" + args.tissue + "_DNase-seq_hg38.bed.gz"
bed_ep300_path = ref_dir + "/" + args.tissue + "_EP300_hg38.bed.gz"
bed_h2afz_path = ref_dir + "/" + args.tissue + "_H2AFZ_hg38.bed.gz"
bed_h3k36me3_path = ref_dir + "/" + args.tissue + "_H3K36me3_hg38.bed.gz"
bed_h3k4me2_path = ref_dir + "/" + args.tissue + "_H3K4me2_hg38.bed.gz"
bed_h3k79me2_path = ref_dir + "/" + args.tissue + "_H3K79me2_hg38.bed.gz"
bed_h3k9ac_path = ref_dir + "/" + args.tissue + "_H3K9ac_hg38.bed.gz"
bed_h3k9me3_path = ref_dir + "/" + args.tissue + "_H3K9me3_hg38.bed.gz"
bed_h4k20me1_path = ref_dir + "/" + args.tissue + "_H4K20me1_hg38.bed.gz"
bed_polr2a_path = ref_dir + "/" + args.tissue + "_POLR2A_hg38.bed.gz"
bed_rad21_path = ref_dir + "/" + args.tissue + "_RAD21_hg38.bed.gz"
bed_smc3_path = ref_dir + "/" + args.tissue + "_SMC3_hg38.bed.gz"

bed_27ac = pybedtools.BedTool(bed_27ac_path)
bed_27me3 = pybedtools.BedTool(bed_27me3_path)
bed_4me1 = pybedtools.BedTool(bed_4me1_path)
bed_4me3 = pybedtools.BedTool(bed_4me3_path)
bed_atac = pybedtools.BedTool(bed_atac_path)
bed_ctcf = pybedtools.BedTool(bed_ctcf_path)
bed_dnase = pybedtools.BedTool(bed_dnase_path)
bed_ep300 = pybedtools.BedTool(bed_ep300_path)
bed_h2afz = pybedtools.BedTool(bed_h2afz_path)
bed_h3k36me3 = pybedtools.BedTool(bed_h3k36me3_path)
bed_h3k4me2 = pybedtools.BedTool(bed_h3k4me2_path)
bed_h3k79me2 = pybedtools.BedTool(bed_h3k79me2_path)
bed_h3k9ac = pybedtools.BedTool(bed_h3k9ac_path)
bed_h3k9me3 = pybedtools.BedTool(bed_h3k9me3_path)
bed_h4k20me1 = pybedtools.BedTool(bed_h4k20me1_path)
bed_polr2a = pybedtools.BedTool(bed_polr2a_path)
bed_rad21 = pybedtools.BedTool(bed_rad21_path)
bed_smc3 = pybedtools.BedTool(bed_smc3_path)

bed_genome = pybedtools.BedTool().from_dataframe(bed_genome)

overlaps_27ac = bed_genome.intersect(bed_27ac, u = True)
overlaps_27ac = overlaps_27ac.sort()
overlaps_27ac = pybedtools.BedTool.to_dataframe(overlaps_27ac)
if overlaps_27ac.shape[0] == 0:
	merged_df_25k25k["27ac"] = "0"
	merged_df_27ac = merged_df_25k25k
else:
	overlaps_27ac["27ac"] = ones
	merged_df_27ac = merged_df_25k25k.merge(overlaps_27ac, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_27ac["27ac"] = merged_df_27ac["27ac"].fillna("0")

overlaps_27me3 = bed_genome.intersect(bed_27me3, u = True)
overlaps_27me3 = overlaps_27me3.sort()
overlaps_27me3 = pybedtools.BedTool.to_dataframe(overlaps_27me3)
if overlaps_27me3.shape[0] == 0:
	merged_df_27ac["27me3"] = "0"
	merged_df_27me3 = merged_df_27ac
else:
	overlaps_27me3["27me3"] = ones
	merged_df_27me3 = merged_df_27ac.merge(overlaps_27me3, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_27me3["27me3"] = merged_df_27me3["27me3"].fillna("0")

overlaps_4me1 = bed_genome.intersect(bed_4me1, u = True)
overlaps_4me1 = overlaps_4me1.sort()
overlaps_4me1 = pybedtools.BedTool.to_dataframe(overlaps_4me1)
if overlaps_4me1.shape[0] == 0:
	merged_df_27me3["4me1"] = "0"
	merged_df_4me1 = merged_df_27me3
else:
	overlaps_4me1["4me1"] = ones
	merged_df_4me1 = merged_df_27me3.merge(overlaps_4me1, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_4me1["4me1"] = merged_df_4me1["4me1"].fillna("0")

overlaps_4me3 = bed_genome.intersect(bed_4me3, u = True)
overlaps_4me3 = overlaps_4me3.sort()
overlaps_4me3 = pybedtools.BedTool.to_dataframe(overlaps_4me3)
if overlaps_4me3.shape[0] == 0:
	merged_df_4me1["4me3"] = "0"
	merged_df_4me3 = merged_df_4me1
else:
	overlaps_4me3["4me3"] = ones
	merged_df_4me3 = merged_df_4me1.merge(overlaps_4me3, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_4me3["4me3"] = merged_df_4me3["4me3"].fillna("0")

overlaps_atac = bed_genome.intersect(bed_atac, u = True)
overlaps_atac = overlaps_atac.sort()
overlaps_atac = pybedtools.BedTool.to_dataframe(overlaps_atac)
if overlaps_atac.shape[0] == 0:
	merged_df_4me3["atac"] = "0"
	merged_df_atac = merged_df_4me1
else:
	overlaps_atac["atac"] = ones
	merged_df_atac = merged_df_4me3.merge(overlaps_atac, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_atac["atac"] = merged_df_atac["atac"].fillna("0")

overlaps_ctcf = bed_genome.intersect(bed_ctcf, u = True)
overlaps_ctcf = overlaps_ctcf.sort()
overlaps_ctcf = pybedtools.BedTool.to_dataframe(overlaps_ctcf)
if overlaps_ctcf.shape[0] == 0:
	merged_df_atac["ctcf"] = "0"
	merged_df_ctcf = merged_df_atac
else:
	overlaps_ctcf["ctcf"] = ones
	merged_df_ctcf = merged_df_atac.merge(overlaps_ctcf, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_ctcf["ctcf"] = merged_df_ctcf["ctcf"].fillna("0")

overlaps_dnase = bed_genome.intersect(bed_dnase, u = True)
overlaps_dnase = overlaps_dnase.sort()
overlaps_dnase = pybedtools.BedTool.to_dataframe(overlaps_dnase)
if overlaps_dnase.shape[0] == 0:
	merged_df_ctcf["dnase"] = "0"
	merged_df_dnase = merged_df_ctcf
else:
	overlaps_dnase["dnase"] = ones
	merged_df_dnase = merged_df_ctcf.merge(overlaps_dnase, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_dnase["dnase"] = merged_df_dnase["dnase"].fillna("0")

overlaps_ep300 = bed_genome.intersect(bed_ep300, u = True)
overlaps_ep300 = overlaps_ep300.sort()
overlaps_ep300 = pybedtools.BedTool.to_dataframe(overlaps_ep300)
if overlaps_ep300.shape[0] == 0:
	merged_df_dnase["ep300"] = "0"
	merged_df_ep300 = merged_df_dnase
else:
	overlaps_ep300["ep300"] = ones
	merged_df_ep300 = merged_df_dnase.merge(overlaps_ep300, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_ep300["ep300"] = merged_df_ep300["ep300"].fillna("0")

overlaps_h2afz = bed_genome.intersect(bed_h2afz, u = True)
overlaps_h2afz = overlaps_h2afz.sort()
overlaps_h2afz = pybedtools.BedTool.to_dataframe(overlaps_h2afz)
if overlaps_h2afz.shape[0] == 0:
	merged_df_ep300["h2afz"] = "0"
	merged_df_h2afz = merged_df_ep300
else:
	overlaps_h2afz["h2afz"] = ones
	merged_df_h2afz = merged_df_ep300.merge(overlaps_h2afz, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h2afz["h2afz"] = merged_df_h2afz["h2afz"].fillna("0")

overlaps_h3k36me3 = bed_genome.intersect(bed_h3k36me3, u = True)
overlaps_h3k36me3 = overlaps_h3k36me3.sort()
overlaps_h3k36me3 = pybedtools.BedTool.to_dataframe(overlaps_h3k36me3)
if overlaps_h3k36me3.shape[0] == 0:
	merged_df_h2afz["h3k36me3"] = "0"
	merged_df_h3k36me3 = merged_df_h2afz
else:
	overlaps_h3k36me3["h3k36me3"] = ones
	merged_df_h3k36me3 = merged_df_h2afz.merge(overlaps_h3k36me3, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h3k36me3["h3k36me3"] = merged_df_h3k36me3["h3k36me3"].fillna("0")

overlaps_h3k4me2 = bed_genome.intersect(bed_h3k4me2, u = True)
overlaps_h3k4me2 = overlaps_h3k4me2.sort()
overlaps_h3k4me2 = pybedtools.BedTool.to_dataframe(overlaps_h3k4me2)
if overlaps_h3k4me2.shape[0] == 0:
	merged_df_h3k36me3["h3k4me2"] = "0"
	merged_df_h3k4me2 = merged_df_h3k36me3
else:
	overlaps_h3k4me2["h3k4me2"] = ones
	merged_df_h3k4me2 = merged_df_h3k36me3.merge(overlaps_h3k4me2, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h3k4me2["h3k4me2"] = merged_df_h3k4me2["h3k4me2"].fillna("0")

overlaps_h3k79me2 = bed_genome.intersect(bed_h3k79me2, u = True)
overlaps_h3k79me2 = overlaps_h3k79me2.sort()
overlaps_h3k79me2 = pybedtools.BedTool.to_dataframe(overlaps_h3k79me2)
if overlaps_h3k79me2.shape[0] == 0:
	merged_df_h3k4me2["h3k79me2"] = "0"
	merged_df_h3k79me2 = merged_df_h3k4me2
else:
	overlaps_h3k79me2["h3k79me2"] = ones
	merged_df_h3k79me2 = merged_df_h3k4me2.merge(overlaps_h3k79me2, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h3k79me2["h3k79me2"] = merged_df_h3k79me2["h3k79me2"].fillna("0")

overlaps_h3k9ac = bed_genome.intersect(bed_h3k9ac, u = True)
overlaps_h3k9ac = overlaps_h3k9ac.sort()
overlaps_h3k9ac = pybedtools.BedTool.to_dataframe(overlaps_h3k9ac)
if overlaps_h3k9ac.shape[0] == 0:
	merged_df_h3k79me2["h3k9ac"] = "0"
	merged_df_h3k9ac = merged_df_h3k79me2
else:
	overlaps_h3k9ac["h3k9ac"] = ones
	merged_df_h3k9ac = merged_df_h3k79me2.merge(overlaps_h3k9ac, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h3k9ac["h3k9ac"] = merged_df_h3k9ac["h3k9ac"].fillna("0")

overlaps_h3k9me3 = bed_genome.intersect(bed_h3k9me3, u = True)
overlaps_h3k9me3 = overlaps_h3k9me3.sort()
overlaps_h3k9me3 = pybedtools.BedTool.to_dataframe(overlaps_h3k9me3)
if overlaps_h3k9me3.shape[0] == 0:
	merged_df_h3k9ac["h3k9me3"] = "0"
	merged_df_h3k9me3 = merged_df_h3k9ac
else:
	overlaps_h3k9me3["h3k9me3"] = ones
	merged_df_h3k9me3 = merged_df_h3k9ac.merge(overlaps_h3k9me3, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h3k9me3["h3k9me3"] = merged_df_h3k9me3["h3k9me3"].fillna("0")

overlaps_h4k20me1 = bed_genome.intersect(bed_h4k20me1, u = True)
overlaps_h4k20me1 = overlaps_h4k20me1.sort()
overlaps_h4k20me1 = pybedtools.BedTool.to_dataframe(overlaps_h4k20me1)
if overlaps_h4k20me1.shape[0] == 0:
	merged_df_h3k9me3["h4k20me1"] = "0"
	merged_df_h4k20me1 = merged_df_h3k9me3
else:
	overlaps_h4k20me1["h4k20me1"] = ones
	merged_df_h4k20me1 = merged_df_h3k9me3.merge(overlaps_h4k20me1, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_h4k20me1["h4k20me1"] = merged_df_h4k20me1["h4k20me1"].fillna("0")

overlaps_polr2a = bed_genome.intersect(bed_polr2a, u = True)
overlaps_polr2a = overlaps_polr2a.sort()
overlaps_polr2a = pybedtools.BedTool.to_dataframe(overlaps_polr2a)
if overlaps_polr2a.shape[0] == 0:
	merged_df_h4k20me1["polr2a"] = "0"
	merged_df_polr2a = merged_df_h4k20me1
else:
	overlaps_polr2a["polr2a"] = ones
	merged_df_polr2a = merged_df_h4k20me1.merge(overlaps_polr2a, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_polr2a["polr2a"] = merged_df_polr2a["polr2a"].fillna("0")

overlaps_rad21 = bed_genome.intersect(bed_rad21, u = True)
overlaps_rad21 = overlaps_rad21.sort()
overlaps_rad21 = pybedtools.BedTool.to_dataframe(overlaps_rad21)
if overlaps_rad21.shape[0] == 0:
	merged_df_polr2a["rad21"] = "0"
	merged_df_rad21 = merged_df_polr2a
else:
	overlaps_rad21["rad21"] = ones
	merged_df_rad21 = merged_df_polr2a.merge(overlaps_rad21, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_rad21["rad21"] = merged_df_rad21["rad21"].fillna("0")

overlaps_smc3 = bed_genome.intersect(bed_smc3, u = True)
overlaps_smc3 = overlaps_smc3.sort()
overlaps_smc3 = pybedtools.BedTool.to_dataframe(overlaps_smc3)
if overlaps_smc3.shape[0] == 0:
	merged_df_rad21["smc3"] = "0"
	merged_df_smc3 = merged_df_rad21
else:
	overlaps_smc3["smc3"] = ones
	merged_df_smc3 = merged_df_rad21.merge(overlaps_smc3, how = "left", left_on = ["chrom", "start", "end"], right_on = ["chrom", "start", "end"])
	merged_df_smc3["smc3"] = merged_df_smc3["smc3"].fillna("0")

#Add the tissue type in the final column
merged_df_smc3["tissue"] = args.tissue

merged_df_smc3.to_csv(args.tissue + "_epigenome_hg38.txt", sep = "\t", header = True, index = False)