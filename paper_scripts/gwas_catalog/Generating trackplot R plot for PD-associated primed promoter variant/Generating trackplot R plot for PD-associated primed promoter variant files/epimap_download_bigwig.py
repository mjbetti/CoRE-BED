###Developed by Michael J Betti, April 2021, updated 8 July 2022
__author__ = "Michael J Betti"
__license__ = "MIT"
###Developed by Michael J Betti, April 2021, updated 25 February 2022
__author__ = "Michael J Betti"
__copyright__ = "Copyright 2021, Michael J Betti"
__license__ = "BSD"
__maintainer__ = "Michael J Betti"
__email__ = "mjbetti3@gmail.com"
__status__ = "Development"

import os, sys, requests, argparse, pybedtools, pandas as pd
from itertools import chain
import subprocess

parser = argparse.ArgumentParser(add_help = True)
#parser.add_argument("-i", "--input", type = str, required = True, help = "the input bed file (required)")
parser.add_argument("-s", "--separator", type = str, required = False, help = "the upstream boundary distance from a TSS (default: 5000 bp)", default = "\t")
parser.add_argument("-g", "--ref_genome", type = str, required = True, help = "the human or mouse reference genome build on which the input coordinates are based (required) (valid options: GRCh38/hg38, GRCh37/hg19, GRCm39/mm39, GRCm38/mm10, or GRCm37/mm9)")
parser.add_argument("-t", "--tissue", type = str, required = True, help = "the tissue of interest (required) (valid hg19 multi-sample options: Adipose, Blood, Bone, Brain, Cancer, Digestive, Endocrine, Endothelial, ESC_deriv, ESC, Eye, Heart, HSC_and_B_cell, iPSC, Kidney, Liver, Lung, Lymphoblastoid, Mesench, Muscle, Myosat, Neurosph, Other, Pancreas, Placenta_and_EMM, PNS, Reproductive, Sm_muscle, Spleen, Stromal, Thymus, Urinary, User_provided_files, User_provided_urls; valid hg19 single-sample options: adipose_adipocyte, adipose_34m, adipose_omental_fat_pad_51f, adipose_omental_fat_pad_53f, adipose_subcutaneous_25f, adipose_subcutaneous_41f, adipose_subcutaneous_49f, adipose_subcutaneous_59f, adipose_subcutaneous_81f, blood_cd4_alpha_beta_memory_t_cell, blood_cd4_alpha_beta_memory_t_cell_blood_origin, blood_cd4_alpha_beta_t_cell, blood_cd4_alpha_beta_t_cell_33f, blood_cd4_alpha_beta_t_cell_21m, blood_cd4_alpha_beta_t_cell_37m, blood_cd4_alpha_beta_t_cell_phorbol, blood_cd4_cd25_alpha_beta_t_cell, blood_cd8_alpha_beta_memory_t_cell, blood_cd8_alpha_beta_t_cell, blood_cd8_alpha_beta_t_cell_33f, blood_cd8_alpha_beta_t_cell_34f, blood_cd8_alpha_beta_t_cell_21m, blood_cd8_alpha_beta_t_cell_28m, blood_cd8_alpha_beta_t_cell_37m, blood_effector_memory_cd4_alpha_beta_t_cell, blood_mononuclear_cell_male, blood_naive_thymus_cd4_alpha_beta_t_cell, blood_naive_thymus_cd4_alpha_beta_t_cell_35f, blood_naive_thymus_cd4_alpha_beta_t_cell_26m, blood_peripheral_mononuclear_cell_28f, blood_peripheral_mononuclear_cell_27m, blood_peripheral_mononuclear_cell_28m, blood_peripheral_mononuclear_cell_32m, blood_peripheral_mononuclear_cell_39m, blood_regulatory_t_cell_35f, blood_regulatory_t_cell_28m, blood_regulatory_t_cell_blood_origin, blood_t_cell, blood_t_cell_21m, blood_t_cell_36m, blood_t_cell_37m, blood_t1_helper_cell, blood_t1_helper_cell_26f, blood_t1_helper_cell_33m, blood_t17_helper_cell, blood_t17_helper_cell_blood_origin, blood_t17_helper_cell_phorbol, blood_t2_helper_cell_26f, blood_t2_helper_cell_33m, bone_arm, bone_femur, bone_marrow_stroma, bone_leg, bone_osteoblast, brain_ammons_horn_84m, brain_angular_gyrus_75f, brain_angular_gyrus_81m, brain_astrocyte, brain_astrocyte_cerebellum, brain_astrocyte_hippocampus, brain_astrocyte_spinal_cord, brain_embryo_112_days, brain_embryo_56_58_days, brain_embryo_80_days, brain_embryo_105_days_f, brain_embryo_109_days_f, brain_embryo_117_days_f, brain_embryo_120_days_f, brain_embryo_142_days_f, brain_embryo_17_weeks_f, brain_embryo_85_days_f, brain_embryo_96_days_f, brain_embryo_101_days_m, brain_embryo_104_days_m, brain_embryo_105_days_m, brain_embryo_122_days_m, brain_embryo_72_76_days_m, brain_caudate_nucleus_75f, brain_caudate_nucleus_78m, brain_caudate_nucleus_81m, brain_cerebellar_cortex_78_81m, brain_cerebellum_27_35m, brain_cerebellum_53m, brain_cingulate_gyrus_75f, brain_cingulate_gyrus_81m, brain_frontal_cortex_67_80f, brain_frontal_cortex_27_35m, brain_germinal_matrix_20_weeks_m, brain_globus_pallidus_78_84m, brain_inferior_parietal_cortex_84m, brain_hippocampus_75f, brain_hippocampus_73m, brain_medulla_oblongata_78_84m, brain_midbrain_78_84m, brain_middle_frontal_area_75f, brain_middle_frontal_area_81m, brain_middle_frontal_gyrus_78m, brain_occipital_lobe_84m, brain_pons_78m, brain_putamen_78m, brain_substantia_nigra_75f, brain_substantia_nigra_81m, brain_superior_temporal_gyrus_84m, brain_temporal_lobe_75f, brain_temporal_lobe_81m, cancer_prostate_epithelial_22Rv1, cancer_pancreas_adenocarcinoma_8988T, cancer_glioblastoma_A172, cancer_lung_epithelial_carcinoma_A549, cancer_lung_epithelial_carcinoma_A549_treated_ethanol_1hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_1hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_10hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_10min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_12hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_15min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_2hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_20min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_25min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_3hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_30min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_4hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_5hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_5min, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_6hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_7hr, cancer_lung_epithelial_carcinoma_A549_treated_dexamethasone_8hr, cancer_muscle_ewing_sarcoma_A673, cancer_adenoid_cystic_carcinoma_ACC112, cancer_renal_cell_carcinoma_ACHN, cancer_neuroblastoma_BE2C, cancer_prostate_C4-2B, cancer_colorectal_adenocarcinoma_Caco-2, cancer_kidney_clear_cell_carcinoma_Caki2, cancer_myelogenous_leukemia_CMK, cancer_melanoma_COLO829, cancer_desmoplastic_medulloblastoma_Daoy, cancer_acute_lymphoblastic_leukemia_DND-41, cancer_b_cell_lymphoma_DOHH2, cancer_kidney_rhaboid_tumor_G401, cancer_neuroglioma_H4, cancer_glioblastoma_H54, cancer_haploid_myelogenous_leukemia_HAP-1, cancer_colorectal_adenocarcinoma_HCT116, cancer_cervix_adenocarcinoma_HeLa-S3, cancer_cervix_adenocarcinoma_HeLa-S3_G1b_phase, cancer_cervix_adenocarcinoma_HeLa-S3_treated_interferon_alpha_4hr, cancer_hepatocellular_carcinoma_HepG2, cancer_acute_promyelocytic_leukemia_HL-60, cancer_colorectal_adenocarcinoma_HT-29, cancer_fibrosarcoma_HT1080, cancer_hepatocellular_carcinoma_HuH-7.5, cancer_endometrial_adenocarcinoma_Ishikawa_treated_dmso_1hr, cancer_endometrial_adenocarcinoma_Ishikawa_treated_17b-estradiol_30min, cancer_endometrial_adenocarcinoma_Ishikawa_treated_afimoxifene_30min, cancer_myelogenous_leukemia_K562, cancer_myelogenous_leukemia_K562_G1_phase, cancer_myelogenous_leukemia_K562_G2_phase, cancer_myelogenous_leukemia_K562_treated_dmso_72hr, cancer_myelogenous_leukemia_K562_treated_vorinostat_72hr, cancer_myelogenous_leukemia_K562_treated_sodium_butyrate_72hr, cancer_b_cell_lymphoma_Karpas-422, cancer_myelogenous_leukemia_KBM-7, cancer_myeloma_KMS-11, cancer_acute_lymphoblastic_leukemia_KOPT-K1, cancer_prostate_adenocarcinoma_LNCaP_clone_FGC, cancer_prostate_adenocarcinoma_LNCaP_clone_FGC_treated_17b_12hr, cancer_acute_lymphoblastic_leukemia_Loucy, cancer_colorectal_adenocarcinoma_LoVo, cancer_glioblastoma_M059J, cancer_mammary_gland_adenocarcinoma_MCF-7, cancer_mammary_gland_adenocarcinoma_MCF-7_originated_from_MCF-7, cancer_mammary_gland_adenocarcinoma_MCF-7_treated_lactate_24hr, cancer_mammary_gland_adenocarcinoma_MCF-7_treated_estradiol_1hr, cancer_mammary_gland_adenocarcinoma_MCF-7_treated_ctcf_shrna_knockown, cancer_medulloblastoma, cancer_osteosarcoma_MG63, cancer_burkitt_lymphoma_NAMALWA, cancer_burkitt_lymphoma_NAMALWA_treated_sendai_virus_2hr, cancer_burkitt_lymphoma_NAMALWA_treated_sendai_virus_2hr, cancer_acute_promyelocytic_leukemia_NB4, cancer_squamous_cell_carcinoma_NCI-H226, cancer_large_cell_lung_NCI-H460, cancer_myeloma_NCI-H929, cancer_testicular_embryonal_carcinoma_NT2_D1, cancer_b_cell_lymphoma_OCI-LY1, cancer_b_cell_lymphoma_OCI-LY3, cancer_b_cell_lymphoma_OCI-LY7, cancer_pancreas_duct_epithelial_carcinoma_Panc1, cancer_parathyroid_adenoma_62m, cancer_prostate_adenocarcinoma_PC-3, cancer_lung_adenocarcinoma_PC-9, cancer_renal_cell_adenocarcinoma_RCC_7860, cancer_renal_cell_carcinoma, cancer_colon_carcinoma_RKO, cancer_melanoma_RPMI-7951, cancer_plasma_cell_myeloma_RPMI8226, cancer_rhabdomyosarcoma_SJCRH30, cancer_osteosarcoma_SJSA1, cancer_melanoma_SK-MEL-5, cancer_neuroblastoma_SK-N-DZ, cancer_neuroblastoma_SK-N-DZ_treated_dmso_72hr, cancer_neuroepithelioma_SK-N-MC, cancer_neuroblastoma_SK-N-SH, cancer_neuroblastoma_SK-N-SH_treated_retinoic_acid_48hr, cancer_b_cell_lymphoma_SU-DHL-6, cancer_colorectal_adenocarcinoma_SW480, cancer_mammary_gland_ductal_carcinoma_T47D, cancer_mammary_gland_ductal_carcinoma_T47D_treated_17b-estradiol_30min, cancer_prostate_epithelial_carcinoma_VCaP, cancer_eye_retinoblastoma_WERI-Rb-1, digestive_colon_mucosa_56f, digestive_colon_mucosa_73f, digestive_duodenum_mucosa_59m, digestive_duodenum_mucosa_76m, digestive_esophagus_30f, digestive_esophagus_34m, digestive_esophagus_muscularis_mucosa_53f, digestive_esophagus_muscularis_mucosa_37m, digestive_esophagus_squamous_epithelium_53f, digestive_esophagus_squamous_epithelium_37m, digestive_gastroesophageal_sphincter_53f, digestive_large_intestine_embryo_103_days_f, digestive_large_intestine_embryo_105_days_f, digestive_large_intestine_embryo_107_days_f, digestive_large_intestine_embryo_108_days_f, digestive_large_intestine_embryo_110_days_f, digestive_large_intestine_embryo_120_days_f, digestive_large_intestine_embryo_91_days_f, digestive_large_intestine_embryo_98_days_f, digestive_large_intestine_embryo_105_days_m, digestive_large_intestine_embryo_108_days_m, digestive_large_intestine_embryo_113_days_m, digestive_large_intestine_embryo_115_days_m, digestive_large_intestine_embryo_91_days_m, digestive_rectum_mucosa_50f, digestive_rectum_mucosa_61f, digestive_rectum_mucosa_59m, digestive_stomach_mucosa_59m, digestive_peyers_patch_53f, digestive_peyers_patch_37m, digestive_peyers_patch_54m, digestive_sigmoid_colon_51f, digestive_sigmoid_colon_53f, digestive_sigmoid_colon_34m, digestive_sigmoid_colon_54m, digestive_sigmoid_colon_3m, digestive_small_intestine_30f, digestive_small_intestine_105_days_f, digestive_small_intestine_107_days_f, digestive_small_intestine_108_days_f, digestive_small_intestine_110_days_f, digestive_small_intestine_120_days_f, digestive_small_intestine_91_days_f, digestive_small_intestine_98_days_f, digestive_small_intestine_34m, digestive_small_intestine_3m, digestive_small_intestine_105_days_m, digestive_small_intestine_108_days_m, digestive_small_intestine_115_days_m, digestive_small_intestine_87_days_m, digestive_small_intestine_91_days_m, digestive_stomach_101_days, digestive_stomach_30f, digestive_stomach_51f, digestive_stomach_53f, digestive_stomach_f, digestive_stomach_105_days_f, digestive_stomach_107_days_f, digestive_stomach_108_days_f, digestive_stomach_121_days_f, digestive_stomach_147_days_f, digestive_stomach_96_days_f, digestive_stomach_98_days_f, digestive_stomach_34m, digestive_stomach_54m, digestive_stomach_3m, digestive_stomach_108_days_m, digestive_stomach_127_days_m, digestive_stomach_58_76_days_m, digestive_stomach_91_days_m, digestive_transverse_colon_51f, digestive_transverse_colon_53f, digestive_transverse_colon_37m, digestive_transverse_colon_54m, endocrine_adrenal_gland_96_days, endocrine_adrenal_gland_30f, endocrine_adrenal_gland_51f, endocrine_adrenal_gland_53f, endocrine_adrenal_gland_108_days_f, endocrine_adrenal_gland_113_days_f, endocrine_adrenal_gland_85_days_f, endocrine_adrenal_gland_34m, endocrine_adrenal_gland_37m, endocrine_adrenal_gland_54m, endocrine_adrenal_gland_101_days_m, endocrine_adrenal_gland_108_days_m, endocrine_adrenal_gland_85_days_m, endocrine_adrenal_gland_97_days_m, endocrine_pancreas_59, endocrine_pancreas_m, endocrine_pancreas_45m, endocrine_pancreas_46m, endocrine_ovary_30f, endocrine_ovary_51f, endocrine_ovary_53f, endocrine_ovary_embryo_f, endocrine_testis_37m, endocrine_testis_54m, endocrine_testis_embryo_m, endocrine_thyroid_gland_51f, endocrine_thyroid_gland_53f, endocrine_thyroid_gland_37m, endocrine_thyroid_gland_54m, endothelial_brain_microvascular, endothelial_dermis_blood_vessel_adult_f, endothelial_dermis_blood_vessel_newborn_m, endothelial_dermis_microvascular_lymphatoc_vessel_f, endothelial_dermis_microvascular_lymphatoc_vessel_m, endothelial_umbilical_vein_newborn_m, endothelial_umbilical_vein_newborn, endothelial_glomerulus, endothelial_kidney_capillary_113_days_f, endothelial_lung_microvascular_f, endothelial_pulmonary_artery_f, epithelial_amnion, epithelial_bronchial, epithelial_bronchial_f_treated_retinoic_acid, epithelial_choroid_plexus, epithelial_colon, epithelial_esophagus, epithelial_prostate, epithelial_prostate_m, epithelial_proximal_tubule, epithelial_foreskin_keratinocyte_newborn_m, epithelial_foreskin_keratinocyte_newborn_2-4_days_m, epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_2.5_days, epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_5.5_days, epithelial_glomerulus_visceral_3, epithelial_proximal_tubule_HK-2, epithelial_pancreatic_duct_dpde6-e6e7, epithelial_bone_marrow_hs-27a, epithelial_iris_pigment, epithelial_keratinocyte_f, epithelial_keratinocyte_m, epithelial_kidney, epithelial_glomerulus_43_62_m, epithelial_tubule_80f_62m, epithelial_tubule_80f_treated_cisplatin, epithelial_skin_leg_53f, epithelial_skin_leg_37m, epithelial_mammary_luminal_33f, epithelial_mammary_f, epithelial_mammary_18f, epithelial_mammary_50f, epithelial_breast_mcf_10a, epithelial_breast_mcf_10a_treated_tamoxifen_24hr, epithelial_breast_mcf_10a_treated_tamoxifen_6hr, epithelial_mammary_myoepithelial_33f, epithelial_mammary_myoepithelial_36f, epithelial_non-pigmented_ciliary, epithelial_renal_cortical, epithelial_retinal, epithelial_prostate_rwpe1, epithelial_prostate_rwpe2, epithelial_skin_of_body_82_days_f, esc_deriv_bipolar_neuron_gm23338_origin_treated_doxycycline_4_days, esc_deriv_cardiac_mesoderm_h7-hesc_origin, esc_deriv_cardiac_muscle_rues2_origin, esc_deriv_neural_progenitor_h9_origin_1, esc_deriv_ectodermal, esc_deriv_endodermal, esc_deriv_endodermal_hues64_origin, esc_deriv_hepatocyte_h9_origin, esc_deriv_mesenchymal_stem_h1-hesc_origin, esc_deriv_mesendoderm_h1-hesc_origin, esc_deriv_mesodermal_hues64_origin, esc_deriv_neural_crest_h1-hesc_origin, esc_deriv_neural_progenitor_h9_origin_2, esc_deriv_neural_progenitor_h1-hesc_origin, esc_deriv_neural_progenitor_h9_origin_3, esc_deriv_neuron_h9_origin, esc_deriv_smooth_muscle_h9_origin, esc_deriv_trophoblast_h1-hesc_origin, esc_elf-1, esc_es-I3, esc_h1-hesc, esc_h7-hesc, esc_h9, esc_hues48, esc_hues6, esc_hues64, esc_ucsf-4, eye_56_days_76_days_m, eye_76_days_f, eye_125_days_103_days_m, eye_74_days_85_days, eye_89_days_f, heart_aorta_30f, heart_aorta_34m, heart_ascending_aorta_51f, heart_ascending_aorta_53f, heart_coronary_artery_51f, heart_coronary_artery_53f, heart_101_days, heart_59_days_76_days_f, heart_80_days, heart_96_days, heart_103_days_f, heart_105_days_f, heart_110_days_f, heart_116_days_98_days_f, heart_116_days_117_days_f, heart_147_days_f, heart_91_days_f, heart_left_ventricle_53f, heart_left_ventricle_101_days_103_days_f, heart_left_ventricle_136_days_f, heart_left_ventricle_34m, heart_left_ventricle_3m, heart_27m_35m, heart_3m, heart_105_days_m, heart_110_days_m, heart_120_days_m, heart_72_days_76_days_m, heart_91_days_m, heart_96_days_m, heart_right_ventricle_101_days_103_days_f, heart_right_ventricle_34m, heart_right_ventricle_3m, heart_left_atrium_101_days_f, heart_right_atrium_51f, heart_right_atrium_53f, heart_right_atrium_34m, heart_thoracic_aorta_37m, heart_thoracic_aorta_54m, heart_tibial_artery_53f, heart_tibial_artery_37m, hsc_and_b_cell_b_cell, hsc_and_b_cell_b_cell_27f, hsc_and_b_cell_b_cell_27f_43f, hsc_and_b_cell_b_cell_34f, hsc_and_b_cell_b_cell_43f, hsc_and_b_cell_b_cell_21m, hsc_and_b_cell_b_cell_37m, hsc_and_b_cell_cd14_monocyte_f, hsc_and_b_cell_cd14_monocyte_34f, hsc_and_b_cell_cd14_monocyte_21m, hsc_and_b_cell_cd14_monocyte_37m, hsc_and_b_cell_cd1c_myeloid_dendritic, hsc_and_b_cell_cmp_cd34, hsc_and_b_cell_cmp_cd34_f, hsc_and_b_cell_cmp_cd34_27f, hsc_and_b_cell_cmp_cd34_33f, hsc_and_b_cell_cmp_cd34_50f, hsc_and_b_cell_cmp_cd34_m, hsc_and_b_cell_cmp_cd34_adult_m, hsc_and_b_cell_cmp_cd34_23m, hsc_and_b_cell_cmp_cd34_36m, hsc_and_b_cell_cmp_cd34_37m, hsc_and_b_cell_cmp_cd34_42m, hsc_and_b_cell_cmp_cd34_43m, hsc_and_b_cell_cmp_cd34_49m, hsc_and_b_cell_germinal_center, hsc_and_b_cell_mpp, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_20_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_11_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_13_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_15_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_17_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_18_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_4_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_6_days, hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_8_days, hsc_and_b_cell_lymphocyte_jurkat_clone_e61, hsc_and_b_cell_lymphocyte_naive_b_cell, hsc_and_b_cell_nk_34f, hsc_and_b_cell_nk_21m, hsc_and_b_cell_nk_37m, hsc_and_b_cell_neutrophil, hsc_and_b_cell_neutrophil_m, ipsc_cwru1_m, ipsc_gm23338_53m_gm23248_origin, ipsc_ips_df_19.11_newborn_m, ipsc_ips_df_19.7_newborn_m, ipsc_ips_df_14.7_newborn_m, ipsc_ips_df_6.9_newborn_m, ipsc_ips-11a_36m, ipsc_ips-15b_48f, ipsc_ips-18a_48f, ipsc_ips-18c_48f, ipsc_ips-20b_55m, ipsc_ips-nihi11_71m_ag20443_origin, ipsc_ips-nihi7_85f_ag08395_origin, ipsc_l1-s8, ipsc_l1-s8r, kidney_hek293, kidney_hek293t, kidney_59_days_59_days_f, kidney_80_days, kidney_105_days_f, kidney_108_days_f, kidney_113_days_f, kidney_120_days_f, kidney_121_days_f, kidney_76_days_f_76_days_m, kidney_85_days_f, kidney_50m, kidney_67m, kidney_105_days_m, kidney_85_days_m, kidney_87_days_m, kidney_left_107_days_f, kidney_left_110_days_f, kidney_left_147_days_f, kidney_left_59_days_f_91_days_m, kidney_left_87_days_f, kidney_left_89_days_f, kidney_left_115_days_m, kidney_left_87_days_m, kidney_left_96_days_m, kidney_left_renal_cortex_interstitium_105_days_m, kidney_left_renal_cortex_interstitium_120_days_m, kidney_left_renal_pelvis_105_days_m, kidney_left_renal_pelvis_120_days_m, kidney_renal_cortex_interstitium_103_days_f, kidney_renal_cortex_interstitium_120_days_f, kidney_renal_cortex_interstitium_89_days_f, kidney_renal_cortex_interstitium_96_days_f, kidney_renal_cortex_interstitium_108_days_m, kidney_renal_cortex_interstitium_113_days_m, kidney_renal_cortex_interstitium_127_days_m, kidney_renal_cortex_interstitium_91_days_m, kidney_renal_cortex_interstitium_97_days_m, kidney_renal_pelvis_103_days_f, kidney_renal_pelvis_105_days_f, kidney_renal_pelvis_89_days_f, kidney_renal_pelvis_96_days_f, kidney_renal_pelvis_108_days_m, kidney_renal_pelvis_113_days_m, kidney_renal_pelvis_127_days_m, kidney_renal_pelvis_91_days_m, kidney_renal_pelvis_97_days_m, kidney_right_107_days_f, kidney_right_117_days_f, kidney_right_147_days_f, kidney_right_87_days_f, kidney_right_98_days_f, kidney_right_108_days_m, kidney_right_115_days_m, kidney_right_87_days_m, kidney_right_91_days_m, kidney_right_96_days_m, kidney_right_renal_cortex_interstitium_105_days_m, kidney_right_renal_cortex_interstitium_120_days_m, kidney_right_renal_pelvis_105_days_m, kidney_right_renal_pelvis_120_days_m, liver_hepatic_stellate_cell_59f, liver_hepatocyte, liver_59_days_80_days, liver_25f, liver_101_days_113_days_f, liver_31m, liver_32m, liver_78m, liver_right_lobe_53f, lung_left_105_days_f, lung_left_107_days_f, lung_left_108_days_f, lung_left_110_days_f, lung_left_117_days_f, lung_left_91_days_f, lung_left_98_days_f, lung_left_105_days_m, lung_left_113_days_m, lung_left_115_days_m, lung_left_87_days_m, lung_left_91_days_m, lung_left_96_days_m, lung_101_days, lung_112_days, lung_67_days, lung_80_days_76_days_m, lung_30f, lung_108_days_f, lung_120_days_f, lung_76_days_f, lung_82_days_f, lung_85_days_f, lung_96_days_f, lung_3m, lung_103_days_m, lung_108_days_m, lung_54_days_58_days_m, lung_82_days_m, lung_right_105_days_f, lung_right_107_days_f, lung_right_108_days_f, lung_right_110_days_f, lung_right_117_days_f, lung_right_91_days_f, lung_right_98_days_f, lung_right_105_days_m, lung_right_115_days_m, lung_right_87_days_m, lung_right_96_days_m, lung_left_upper_lobe_51f, lung_left_upper_lobe_53f, lung_left_upper_lobe_37m, lymphoblastoid_gm06990, lymphoblastoid_gm08714, lymphoblastoid_gm10248, lymphoblastoid_gm10266, lymphoblastoid_gm12864, lymphoblastoid_gm12865, lymphoblastoid_gm12875, lymphoblastoid_gm12878, lymphoblastoid_gm12891, lymphoblastoid_gm12892, lymphoblastoid_gm18507, lymphoblastoid_gm19238, lymphoblastoid_gm19239, lymphoblastoid_gm19240, mesench_adipocyte_msc_origin, mesench_msc_dedifferentiated_amniotic_fluid_origin, mesench_embryonic_facial_prominence_53_days_58_days, mesench_msc_adipose_origin, muscle_cardiac_myocyte, muscle_forelimb_108_days_f, muscle_gastrocnemius_medialis_51f, muscle_gastrocnemius_medialis_53f, muscle_gastrocnemius_medialis_37m, muscle_gastrocnemius_medialis_54m, muscle_leg_hindlimb_120_days_m, muscle_arm_101_days, muscle_arm_105_days_f, muscle_arm_115_days_f, muscle_arm_120_days_f, muscle_arm_85_days_f, muscle_arm_98_days_f, muscle_arm_101_days_m muscle_arm_104_days_m, muscle_arm_105_days_m, muscle_arm_113_days_m, muscle_arm_115_days_m, muscle_arm_120_days_m, muscle_arm_127_days_m, muscle_arm_96_days_m, muscle_arm_97_days_m, muscle_back_105_days_f, muscle_back_113_days_f, muscle_back_115_days_f, muscle_back_85_days_f, muscle_back_98_days_f, muscle_back_101_days_m, muscle_back_104_days_m, muscle_back_105_days_m, muscle_back_108_days_m, muscle_back_127_days_m, muscle_back_91_days_m, muscle_back_96_days_m, muscle_back_97_days_m, muscle_leg_105_days_f, muscle_leg_110_days_f, muscle_leg_113_days_f, muscle_leg_115_days_f, muscle_leg_85_days_f, muscle_leg_101_days_m, muscle_leg_104_days_m, muscle_leg_105_days_m, muscle_leg_115_days_m, muscle_leg_127_days_m, muscle_leg_96_days_m, muscle_leg_97_days_m, muscle_trunk_113_days_f, muscle_trunk_115_days_f, muscle_trunk_120_days_f, muscle_trunk_121_days_f, muscle_psoas_30f, muscle_psoas_27m_35m, muscle_psoas_34m, muscle_psoas_3m, muscle_skeletal_cell, muscle_skeletal_tissue, muscle_skeletal_tissue_72f, muscle_skeletal_tissue_54m, muscle_tongue_59_days_f_76_days_f, muscle_tongue_72_days_m, myosat_skeletal_muscle_myoblast_lhcn-m2, myosat_myocyte_lhcn-m2_origin, myosat_myotube_skeletal_muscle_myoblast_origin, myosat_skeletal_muscle_myoblast, myosat_skeletal_muscle_myoblast_22m, myosat_skeletal_muscle_satellite_cell_mesoderm_origin_f, neurosph_15_weeks_ganglionic_eminence_origin, neurosph_17_weeks_f_ganglionic_cortex_origin, neurosph_17_weeks_f_ganglionic_eminence_origin, neurosph_olfactory_cell_line, other_breast_epithelium_51f, other_breast_epithelium_53f, other_epidermal_melanocyte, other_foreskin_melanocyte_newborn_m, other_limb_embryo_53_days_56_days, other_limb_embryo_58_days_59_days, other_mammary_stem_cell, pancreas_body_51f, pancreas_body_53f, pancreas_body_37m, pancreas_body_54m, pancreas_islet_precursor_cell, pancreas_30f, pancreas_34m, placenta_and_eem_amnion_16_weeks_m, placenta_and_eem_amnion_stem_cell, placenta_and_eem_chorion, placenta_and_eem_chorion_40_weeks_f, placenta_and_eem_chorion_16_weeks_m, placenta_and_eem_chorionic_villus_16_weeks, placenta_and_eem_chorionic_villus_40_weeks_f, placenta_and_eem_chorionic_villus_16_weeks_m, placenta_and_eem_chorionic_villus_38_weeks_m, placenta_and_eem_trophoblast_htr-8_svneo, placenta_and_eem_placenta_102_days, placenta_and_eem_placenta_16_weeks, placenta_and_eem_placenta_53_days, placenta_and_eem_placenta_56_days_59_days, placenta_and_eem_placenta_101_days_f_105_days_m, placenta_and_eem_placenta_105_days_f, placenta_and_eem_placenta_108_days_f, placenta_and_eem_placenta_113_days_f, placenta_and_eem_placenta_85_days_f, placenta_and_eem_placenta_16_weeks_m, placenta_and_eem_placenta_85_days_m, placenta_and_eem_placenta_91_days_m, placenta_and_eem_placental_basal_plate_40_weeks_f, placenta_and_eem_placental_basal_plate_38_weeks_m, placenta_and_eem_trophoblast_cell_17_weeks_18_weeks, placenta_and_eem_trophoblast_cell_21_weeks, placenta_and_eem_trophoblast_cell_23_weeks, placenta_and_eem_trophoblast_cell_39_weeks_40_weeks, placenta_and_eem_trophoblast_20_weeks_f, placenta_and_eem_trophoblast_40_weeks_f, placenta_and_eem_umbilical_cord_59_days_76_days_m, pns_spinal_cord_108_days_f, pns_spinal_cord_113_days_f, pns_spinal_cord_59_days_f_72_days_m, pns_spinal_cord_87_days_f, pns_spinal_cord_89_days_f, pns_spinal_cord_105_days_m, pns_spinal_cord_96_days_m, pns_tibial_nerve_51f, pns_tibial_nerve_53f, pns_tibial_nerve_37m, reproductive_prostate_gland_37m, reproductive_prostate_gland_54m, reproductive_uterus_53f, reproductive_vagina_51f, reproductive_vagina_53f, sm_muscle_colon_56f, sm_muscle_colon_77f, sm_muscle_duodenum_59m, sm_muscle_duodenum_73m, sm_muscle_rectum_50f, sm_muscle_brain_vasculature_smooth_cell, sm_muscle_stomach_84f, sm_muscle_stomach_59m, spleen_112_days, spleen_30f, spleen_53f, spleen_34m, spleen_54m, spleen_3m, stromal_skin_fibroblast_ag04449, stromal_lung_fibroblast_ag04450, stromal_skin_fibroblast_ag08395, stromal_lung_fibroblast_ag08396, stromal_skin_fibroblast_ag08396, stromal_gingival_fibroblast_ag09319, stromal_skin_fibroblast_ag10803, stromal_skin_fibroblast_ag20443, stromal_skin_fibroblast_bj, stromal_brain_pericyte, stromal_cardiac_fibroblast, stromal_cardiac_fibroblast_f, stromal_cardiac_fibroblast_94_days_f_98_days_f, stromal_skin_fibroblast_eh, stromal_skin_fibroblast_el, stromal_skin_fibroblast_elr, stromal_breast_fibroblast_17f, stromal_breast_fibroblast_26f, stromal_dermis_fibroblast, stromal_dermis_fibroblast_f, stromal_dermis_fibroblast_none_f, stromal_gingival_fibroblast, stromal_lung_fibroblast, stromal_lung_fibroblast_11f_45m, stromal_lung_fibroblast_45m, stromal_mammary_fibroblast_f, stromal_peridontal_ligament_fibroblast_m, stromal_pulmonary_artery_fibroblast, stromal_skin_fibroblast_abdomen_97_days_m, stromal_aorta_fibroblast_f, stromal_conjunctiva_fibroblast, stromal_villous_mesenchyme_fibroblast, stromal_foreskin_fibroblast_newborn_m, stromal_skin_fibroblast_gm03348, stromal_skin_fibroblast_gm04503, stromal_skin_fibroblast_gm04504, stromal_skin_fibroblast_gm23248, stromal_foreskin_fibroblast_hff-myc_foreskin_fibroblast_origin, stromal_lung_fibroblast_imr-90, stromal_lung_fibroblast_97_days_m, stromal_bone_marrow_m, stromal_lung_fibroblast_wi38, thymus_embryo_f, thymus_105_days_f, thymus_110_days_f, thymus_113_days_f, thymus_147_days_f, thymus_98_days_f, thymus_3m, thymus_104_days_m, thymus_108_days_m, thymus_113_days_m, thymus_127_days_m, urinary_bladder_34m, urinary_bladder_76_days_m, urinary_urothelium_cell_line, all, adipose, blood, bone, brain, cancer, digestive, endocrine, endothelial, epithelial, esc_deriv, esc, eye, heart, hsc_and_b_cell, ipsc, kidney, liver, lung, lymphoblastoid, mesench, muscle, myosat, neurosph, other, pancreas, placenta_and_emm, pns, reproductive, sm_muscle, spleen, stromal, thymus, urinary, User_provided_files, User_provided_urls; valid hg38 multi-sample options: User_provided_files, User_provided_urls; valid mouse options: User_provided_files, User_provided_urls)")
parser.add_argument("-ud", "--tss_distance_upstream", type = int, required = False, help = "the upstream boundary distance from a TSS (default: 5000 bp)", default = 5000)
parser.add_argument("-dd", "--tss_distance_downstream", type = int, required = False, help = "the downstream boundary distance from a TSS (default: 1000 bp)", default = 1000)
parser.add_argument("-o", "--output", type = str, required = False, help = "the name of the output file", default = "out.bed")
parser.add_argument("-r", "--ref_dir", type = str, required = False, help = "the path of the reference file directory", default = "ref_files")
parser.add_argument("--no_multianno", required = False, help = "if a coordinate overlaps with multiple regions, keep the most significant occurance", action = "store_true")
parser.add_argument("--write_summary", required = False, help = "Write out a summary of regulatory counts as a .txt file", action = "store_true")
parser.add_argument("--write_anno_only", required = False, help = "Instead of the input file appended with an annotation column, write out only the annotation column", action = "store_true")
parser.add_argument("--bed_cols", type = str, required = False, help = "if the input is not in traditional UCSC BED format, specify the column numbers of chr, start, and end separated by commas", default = "1,2,3")
parser.add_argument("--input_header", required = False, help = "use if the input file has a header", action = "store_true")
parser.add_argument("--user_4me1", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K4me1 ChIP-seq peaks")
parser.add_argument("--user_4me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K4me3 ChIP-seq peaks")
parser.add_argument("--user_27ac", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K27ac ChIP-seq peaks")
parser.add_argument("--user_27me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K27me3 ChIP-seq peaks")
#parser.add_argument("--user_36me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K36me3 ChIP-seq peaks")
parser.add_argument("--user_ctcf", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided CTCF ChIP-seq peaks")
parser.add_argument("--user_dnase", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided DNase-seq peaks")
parser.add_argument("--user_tissue_names", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify the names of each corresponding tissue type")
parser.add_argument("-v", "--verbose", required = False, help = "return logging as terminal output", action = "store_true")
parser.add_argument("-h3k4me1", "--only_h3k4me1", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K4me1 overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k4me3", "--only_h3k4me3", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K4me3 overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k27ac", "--only_h3k27ac", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K27ac overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k27me3", "--only_h3k27me3", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K27me3 overlap in the target tissue", action = "store_true")
parser.add_argument("-dhs", "--only_dhs", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for DHS overlap in the target tissue", action = "store_true")
parser.add_argument("-ctcf", "--only_ctcf", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for CTCF overlap in the target tissue", action = "store_true")
args = parser.parse_args()

#Check that required arguments are specified
#assert args.input, "Must specify input file (-i, --input)"
assert args.ref_genome, "Must specify reference genome build (-g, --ref_genome)"
assert args.tissue, "Must specify tissue type (-t, --tissue)"
if args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39" or args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38" or args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "grcm37":
	assert args.tissue.lower() == "user_provided_files" or args.tissue.lower() == "user_provided_urls"
if args.tissue.lower() == "user_provided_files" or args.tissue.lower() == "user_provided_urls":
	assert args.user_4me1 and args.user_4me3 and args.user_27ac and args.user_27me3 and args.user_ctcf and args.user_dnase and args.user_tissue_names, "Must provide histone ChIP-seq, CTCF ChIP-seq, and DNase-seq files when using the User_provided_files tissue option or URLs when using the User_provided_urls option. Must also provide the tissue/cell types names of each specified set of files."

#If the reference genome is human, check that specified tissue type is one of the 30 valid options
###Create an array with all of the parsed in tissue types
tissue_array = args.tissue.lower().split(",")
adipose_array = ["adipose_adipocyte", "adipose_34m", "adipose_omental_fat_pad_51f", "adipose_omental_fat_pad_53f", "adipose_subcutaneous_25f", "adipose_subcutaneous_41f", "adipose_subcutaneous_49f", "adipose_subcutaneous_59f", "adipose_subcutaneous_81f"]
blood_array = ["blood_cd4_alpha_beta_memory_t_cell", "blood_cd4_alpha_beta_memory_t_cell_blood_origin", "blood_cd4_alpha_beta_t_cell", "blood_cd4_alpha_beta_t_cell_33f", "blood_cd4_alpha_beta_t_cell_21m", "blood_cd4_alpha_beta_t_cell_37m", "blood_cd4_alpha_beta_t_cell_phorbol", "blood_cd4_cd25_alpha_beta_t_cell", "blood_cd8_alpha_beta_memory_t_cell", "blood_cd8_alpha_beta_t_cell", "blood_cd8_alpha_beta_t_cell_33f", "blood_cd8_alpha_beta_t_cell_34f", "blood_cd8_alpha_beta_t_cell_21m", "blood_cd8_alpha_beta_t_cell_28m", "blood_cd8_alpha_beta_t_cell_37m", "blood_effector_memory_cd4_alpha_beta_t_cell", "blood_mononuclear_cell_male", "blood_naive_thymus_cd4_alpha_beta_t_cell", "blood_naive_thymus_cd4_alpha_beta_t_cell_35f", "blood_naive_thymus_cd4_alpha_beta_t_cell_26m", "blood_peripheral_mononuclear_cell_28f", "blood_peripheral_mononuclear_cell_27m", "blood_peripheral_mononuclear_cell_28m", "blood_peripheral_mononuclear_cell_32m", "blood_peripheral_mononuclear_cell_39m", "blood_regulatory_t_cell_35f", "blood_regulatory_t_cell_28m", "blood_regulatory_t_cell_blood_origin", "blood_t_cell", "blood_t_cell_21m", "blood_t_cell_36m", "blood_t_cell_37m", "blood_t1_helper_cell", "blood_t1_helper_cell_26f", "blood_t1_helper_cell_33m", "blood_t17_helper_cell", "blood_t17_helper_cell_blood_origin", "blood_t17_helper_cell_phorbol", "blood_t2_helper_cell_26f", "blood_t2_helper_cell_33m"]
bone_array = ["bone_arm", "bone_femur", "bone_marrow_stroma", "bone_leg", "bone_osteoblast"]
brain_array = ["brain_ammons_horn_84m", "brain_angular_gyrus_75f", "brain_angular_gyrus_81m", "brain_astrocyte", "brain_astrocyte_cerebellum", "brain_astrocyte_hippocampus", "brain_astrocyte_spinal_cord", "brain_embryo_112_days", "brain_embryo_56_58_days", "brain_embryo_80_days", "brain_embryo_105_days_f", "brain_embryo_109_days_f", "brain_embryo_117_days_f", "brain_embryo_120_days_f", "brain_embryo_142_days_f", "brain_embryo_17_weeks_f", "brain_embryo_85_days_f", "brain_embryo_96_days_f", "brain_embryo_101_days_m", "brain_embryo_104_days_m", "brain_embryo_105_days_m", "brain_embryo_122_days_m", "brain_embryo_72_76_days_m", "brain_caudate_nucleus_75f", "brain_caudate_nucleus_78m", "brain_caudate_nucleus_81m", "brain_cerebellar_cortex_78_81m", "brain_cerebellum_27_35m", "brain_cerebellum_53m", "brain_cingulate_gyrus_75f", "brain_cingulate_gyrus_81m", "brain_frontal_cortex_67_80f", "brain_frontal_cortex_27_35m", "brain_germinal_matrix_20_weeks_m", "brain_globus_pallidus_78_84m", "brain_inferior_parietal_cortex_84m", "brain_hippocampus_75f", "brain_hippocampus_73m", "brain_medulla_oblongata_78_84m", "brain_midbrain_78_84m", "brain_middle_frontal_area_75f", "brain_middle_frontal_area_81m", "brain_middle_frontal_gyrus_78m", "brain_occipital_lobe_84m", "brain_pons_78m", "brain_putamen_78m", "brain_substantia_nigra_75f", "brain_substantia_nigra_81m", "brain_superior_temporal_gyrus_84m", "brain_temporal_lobe_75f", "brain_temporal_lobe_81m"]
cancer_array = ["cancer_prostate_epithelial_22rv1", "cancer_pancreas_adenocarcinoma_8988t", "cancer_glioblastoma_a172", "cancer_lung_epithelial_carcinoma_a549", "cancer_lung_epithelial_carcinoma_a549_treated_ethanol_1hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_1hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_12hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_15min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_2hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_20min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_25min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_3hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_30min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_4hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5min", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_6hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_7hr", "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_8hr", "cancer_muscle_ewing_sarcoma_a673", "cancer_adenoid_cystic_carcinoma_acc112", "cancer_renal_cell_carcinoma_achn", "cancer_neuroblastoma_be2c", "cancer_prostate_c4-2b", "cancer_colorectal_adenocarcinoma_caco-2", "cancer_kidney_clear_cell_carcinoma_caki2", "cancer_myelogenous_leukemia_cmk", "cancer_melanoma_colo829", "cancer_desmoplastic_medulloblastoma_daoy", "cancer_acute_lymphoblastic_leukemia_dnd-41", "cancer_b_cell_lymphoma_dohh2", "cancer_kidney_rhaboid_tumor_g401", "cancer_neuroglioma_h4", "cancer_glioblastoma_h54", "cancer_haploid_myelogenous_leukemia_hap-1", "cancer_colorectal_adenocarcinoma_hct116", "cancer_cervix_adenocarcinoma_hela-s3", "cancer_cervix_adenocarcinoma_hela-s3_g1b_phase", "cancer_cervix_adenocarcinoma_hela-s3_treated_interferon_alpha_4hr", "cancer_hepatocellular_carcinoma_hepg2", "cancer_acute_promyelocytic_leukemia_hl-60", "cancer_colorectal_adenocarcinoma_ht-29", "cancer_fibrosarcoma_ht1080", "cancer_hepatocellular_carcinoma_huh-7.5", "cancer_endometrial_adenocarcinoma_ishikawa_treated_dmso_1hr", "cancer_endometrial_adenocarcinoma_ishikawa_treated_17b-estradiol_30min", "cancer_endometrial_adenocarcinoma_ishikawa_treated_afimoxifene_30min", "cancer_myelogenous_leukemia_k562", "cancer_myelogenous_leukemia_k562_g1_phase", "cancer_myelogenous_leukemia_k562_g2_phase", "cancer_myelogenous_leukemia_k562_treated_dmso_72hr", "cancer_myelogenous_leukemia_k562_treated_vorinostat_72hr", "cancer_myelogenous_leukemia_k562_treated_sodium_butyrate_72hr", "cancer_b_cell_lymphoma_karpas-422", "cancer_myelogenous_leukemia_kbm-7", "cancer_myeloma_kms-11", "cancer_acute_lymphoblastic_leukemia_kopt-k1", "cancer_prostate_adenocarcinoma_lncap_clone_fgc", "cancer_prostate_adenocarcinoma_lncap_clone_fgc_treated_17b_12hr", "cancer_acute_lymphoblastic_leukemia_loucy", "cancer_colorectal_adenocarcinoma_lovo", "cancer_glioblastoma_m059j", "cancer_mammary_gland_adenocarcinoma_mcf-7", "cancer_mammary_gland_adenocarcinoma_mcf-7_originated_from_mcf-7", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_lactate_24hr", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_estradiol_1hr", "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_ctcf_shrna_knockown", "cancer_medulloblastoma", "cancer_osteosarcoma_mg63", "cancer_burkitt_lymphoma_namalwa", "cancer_burkitt_lymphoma_namalwa_treated_sendai_virus_2hr", "cancer_acute_promyelocytic_leukemia_nb4", "cancer_squamous_cell_carcinoma_nci-h226", "cancer_large_cell_lung_nci-h460", "cancer_myeloma_nci-h929", "cancer_testicular_embryonal_carcinoma_nt2_d1", "cancer_b_cell_lymphoma_oci-ly1", "cancer_b_cell_lymphoma_oci-ly3", "cancer_b_cell_lymphoma_oci-ly7", "cancer_pancreas_duct_epithelial_carcinoma_panc1", "cancer_parathyroid_adenoma_62m", "cancer_prostate_adenocarcinoma_pc-3", "cancer_lung_adenocarcinoma_pc-9", "cancer_renal_cell_adenocarcinoma_rcc_7860", "cancer_renal_cell_carcinoma", "cancer_colon_carcinoma_rko", "cancer_melanoma_rpmi-7951", "cancer_plasma_cell_myeloma_rpmi8226", "cancer_rhabdomyosarcoma_sjcrh30", "cancer_osteosarcoma_sjsa1", "cancer_melanoma_sk-mel-5", "cancer_neuroblastoma_sk-n-dz", "cancer_neuroblastoma_sk-n-dz_treated_dmso_72hr", "cancer_neuroepithelioma_sk-n-mc", "cancer_neuroblastoma_sk-n-sh", "cancer_neuroblastoma_sk-n-sh_treated_retinoic_acid_48hr", "cancer_b_cell_lymphoma_su-dhl-6", "cancer_colorectal_adenocarcinoma_sw480", "cancer_mammary_gland_ductal_carcinoma_t47d", "cancer_mammary_gland_ductal_carcinoma_t47d_treated_17b-estradiol_30min", "cancer_prostate_epithelial_carcinoma_vcap", "cancer_eye_retinoblastoma_weri-rb-1"]
digestive_array = ["digestive_colon_mucosa_56f", "digestive_colon_mucosa_73f", "digestive_duodenum_mucosa_59m", "digestive_duodenum_mucosa_76m", "digestive_esophagus_30f", "digestive_esophagus_34m", "digestive_esophagus_muscularis_mucosa_53f", "digestive_esophagus_muscularis_mucosa_37m", "digestive_esophagus_squamous_epithelium_53f", "digestive_esophagus_squamous_epithelium_37m", "digestive_gastroesophageal_sphincter_53f", "digestive_large_intestine_embryo_103_days_f", "digestive_large_intestine_embryo_105_days_f", "digestive_large_intestine_embryo_107_days_f", "digestive_large_intestine_embryo_108_days_f", "digestive_large_intestine_embryo_110_days_f", "digestive_large_intestine_embryo_120_days_f", "digestive_large_intestine_embryo_91_days_f", "digestive_large_intestine_embryo_98_days_f", "digestive_large_intestine_embryo_105_days_m", "digestive_large_intestine_embryo_108_days_m", "digestive_large_intestine_embryo_113_days_m", "digestive_large_intestine_embryo_115_days_m", "digestive_large_intestine_embryo_91_days_m", "digestive_rectum_mucosa_50f", "digestive_rectum_mucosa_61f", "digestive_rectum_mucosa_59m", "digestive_stomach_mucosa_59m", "digestive_peyers_patch_53f", "digestive_peyers_patch_37m", "digestive_peyers_patch_54m", "digestive_sigmoid_colon_51f", "digestive_sigmoid_colon_53f", "digestive_sigmoid_colon_34m", "digestive_sigmoid_colon_54m", "digestive_sigmoid_colon_3m", "digestive_small_intestine_30f", "digestive_small_intestine_105_days_f", "digestive_small_intestine_107_days_f", "digestive_small_intestine_108_days_f", "digestive_small_intestine_110_days_f", "digestive_small_intestine_120_days_f", "digestive_small_intestine_91_days_f", "digestive_small_intestine_98_days_f", "digestive_small_intestine_34m", "digestive_small_intestine_3m", "digestive_small_intestine_105_days_m", "digestive_small_intestine_108_days_m", "digestive_small_intestine_115_days_m", "digestive_small_intestine_87_days_m", "digestive_small_intestine_91_days_m", "digestive_stomach_101_days", "digestive_stomach_30f", "digestive_stomach_51f", "digestive_stomach_53f", "digestive_stomach_f", "digestive_stomach_105_days_f", "digestive_stomach_107_days_f", "digestive_stomach_108_days_f", "digestive_stomach_121_days_f", "digestive_stomach_147_days_f", "digestive_stomach_96_days_f", "digestive_stomach_98_days_f", "digestive_stomach_34m", "digestive_stomach_54m", "digestive_stomach_3m", "digestive_stomach_108_days_m", "digestive_stomach_127_days_m", "digestive_stomach_58_76_days_m", "digestive_stomach_91_days_m", "digestive_transverse_colon_51f", "digestive_transverse_colon_53f", "digestive_transverse_colon_37m", "digestive_transverse_colon_54m"]
endocrine_array = ["endocrine_adrenal_gland_96_days", "endocrine_adrenal_gland_30f", "endocrine_adrenal_gland_51f", "endocrine_adrenal_gland_53f", "endocrine_adrenal_gland_108_days_f", "endocrine_adrenal_gland_113_days_f", "endocrine_adrenal_gland_85_days_f", "endocrine_adrenal_gland_34m", "endocrine_adrenal_gland_37m", "endocrine_adrenal_gland_54m", "endocrine_adrenal_gland_101_days_m", "endocrine_adrenal_gland_108_days_m", "endocrine_adrenal_gland_85_days_m", "endocrine_adrenal_gland_97_days_m", "endocrine_pancreas_59", "endocrine_pancreas_m", "endocrine_pancreas_45m", "endocrine_pancreas_46m", "endocrine_ovary_30f", "endocrine_ovary_51f", "endocrine_ovary_53f", "endocrine_ovary_embryo_f", "endocrine_testis_37m", "endocrine_testis_54m", "endocrine_testis_embryo_m", "endocrine_thyroid_gland_51f", "endocrine_thyroid_gland_53f", "endocrine_thyroid_gland_37m", "endocrine_thyroid_gland_54m"]
endothelial_array = ["endothelial_brain_microvascular", "endothelial_dermis_blood_vessel_adult_f", "endothelial_dermis_blood_vessel_newborn_m", "endothelial_dermis_microvascular_lymphatoc_vessel_f", "endothelial_dermis_microvascular_lymphatoc_vessel_m", "endothelial_umbilical_vein_newborn_m", "endothelial_umbilical_vein_newborn", "endothelial_glomerulus", "endothelial_kidney_capillary_113_days_f", "endothelial_lung_microvascular_f", "endothelial_pulmonary_artery_f"]
epithelial_array = ["epithelial_amnion", "epithelial_bronchial", "epithelial_bronchial_f_treated_retinoic_acid", "epithelial_choroid_plexus", "epithelial_colon", "epithelial_esophagus", "epithelial_prostate", "epithelial_prostate_m", "epithelial_proximal_tubule", "epithelial_foreskin_keratinocyte_newborn_m", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_2.5_days", "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_5.5_days", "epithelial_glomerulus_visceral_3", "epithelial_proximal_tubule_hk-2", "epithelial_pancreatic_duct_dpde6-e6e7", "epithelial_bone_marrow_hs-27a", "epithelial_iris_pigment", "epithelial_keratinocyte_f", "epithelial_keratinocyte_m", "epithelial_kidney", "epithelial_glomerulus_43_62_m", "epithelial_tubule_80f_62m", "epithelial_tubule_80f_treated_cisplatin", "epithelial_skin_leg_53f", "epithelial_skin_leg_37m", "epithelial_mammary_luminal_33f", "epithelial_mammary_f", "epithelial_mammary_18f", "epithelial_mammary_50f", "epithelial_breast_mcf_10a", "epithelial_breast_mcf_10a_treated_tamoxifen_24hr", "epithelial_breast_mcf_10a_treated_tamoxifen_6hr", "epithelial_mammary_myoepithelial_33f", "epithelial_mammary_myoepithelial_36f", "epithelial_non-pigmented_ciliary", "epithelial_renal_cortical", "epithelial_retinal", "epithelial_prostate_rwpe1", "epithelial_prostate_rwpe2", "epithelial_skin_of_body_82_days_f"]
esc_deriv_array = ["esc_deriv_bipolar_neuron_gm23338_origin_treated_doxycycline_4_days", "esc_deriv_cardiac_mesoderm_h7-hesc_origin", "esc_deriv_cardiac_muscle_rues2_origin", "esc_deriv_neural_progenitor_h9_origin_1", "esc_deriv_ectodermal", "esc_deriv_endodermal", "esc_deriv_endodermal_hues64_origin", "esc_deriv_hepatocyte_h9_origin", "esc_deriv_mesenchymal_stem_h1-hesc_origin", "esc_deriv_mesendoderm_h1-hesc_origin", "esc_deriv_mesodermal_hues64_origin", "esc_deriv_neural_crest_h1-hesc_origin", "esc_deriv_neural_progenitor_h9_origin_2", "esc_deriv_neural_progenitor_h1-hesc_origin", "esc_deriv_neural_progenitor_h9_origin_3", "esc_deriv_neuron_h9_origin", "esc_deriv_smooth_muscle_h9_origin", "esc_deriv_trophoblast_h1-hesc_origin"]
esc_array = ["esc_elf-1", "esc_es-i3", "esc_h1-hesc", "esc_h7-hesc", "esc_h9", "esc_hues48", "esc_hues6", "esc_hues64", "esc_ucsf-4"]
eye_array = ["eye_56_days_76_days_m", "eye_76_days_f", "eye_125_days_103_days_m", "eye_74_days_85_days", "eye_89_days_f"]
heart_array = ["heart_aorta_30f", "heart_aorta_34m", "heart_ascending_aorta_51f", "heart_ascending_aorta_53f", "heart_coronary_artery_51f", "heart_coronary_artery_53f", "heart_101_days", "heart_59_days_76_days_f", "heart_80_days", "heart_96_days", "heart_103_days_f", "heart_105_days_f", "heart_110_days_f", "heart_116_days_98_days_f", "heart_116_days_117_days_f", "heart_147_days_f", "heart_91_days_f", "heart_left_ventricle_53f", "heart_left_ventricle_101_days_103_days_f", "heart_left_ventricle_136_days_f", "heart_left_ventricle_34m", "heart_left_ventricle_3m", "heart_27m_35m", "heart_3m", "heart_105_days_m", "heart_110_days_m", "heart_120_days_m", "heart_72_days_76_days_m", "heart_91_days_m", "heart_96_days_m", "heart_right_ventricle_101_days_103_days_f", "heart_right_ventricle_34m", "heart_right_ventricle_3m", "heart_left_atrium_101_days_f", "heart_right_atrium_51f", "heart_right_atrium_53f", "heart_right_atrium_34m", "heart_thoracic_aorta_37m", "heart_thoracic_aorta_54m", "heart_tibial_artery_53f", "heart_tibial_artery_37m"]
hsc_and_b_cell_array = ["hsc_and_b_cell_b_cell", "hsc_and_b_cell_b_cell_27f", "hsc_and_b_cell_b_cell_27f_43f", "hsc_and_b_cell_b_cell_34f", "hsc_and_b_cell_b_cell_43f", "hsc_and_b_cell_b_cell_21m", "hsc_and_b_cell_b_cell_37m", "hsc_and_b_cell_cd14_monocyte_f", "hsc_and_b_cell_cd14_monocyte_34f", "hsc_and_b_cell_cd14_monocyte_21m", "hsc_and_b_cell_cd14_monocyte_37m", "hsc_and_b_cell_cd1c_myeloid_dendritic", "hsc_and_b_cell_cmp_cd34", "hsc_and_b_cell_cmp_cd34_f", "hsc_and_b_cell_cmp_cd34_27f", "hsc_and_b_cell_cmp_cd34_33f", "hsc_and_b_cell_cmp_cd34_50f", "hsc_and_b_cell_cmp_cd34_m", "hsc_and_b_cell_cmp_cd34_adult_m", "hsc_and_b_cell_cmp_cd34_23m", "hsc_and_b_cell_cmp_cd34_36m", "hsc_and_b_cell_cmp_cd34_37m", "hsc_and_b_cell_cmp_cd34_42m", "hsc_and_b_cell_cmp_cd34_43m", "hsc_and_b_cell_cmp_cd34_49m", "hsc_and_b_cell_germinal_center", "hsc_and_b_cell_mpp", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_20_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_11_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_13_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_15_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_17_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_18_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_4_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_6_days", "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_8_days", "hsc_and_b_cell_lymphocyte_jurkat_clone_e61", "hsc_and_b_cell_lymphocyte_naive_b_cell", "hsc_and_b_cell_nk_34f", "hsc_and_b_cell_nk_21m", "hsc_and_b_cell_nk_37m", "hsc_and_b_cell_neutrophil", "hsc_and_b_cell_neutrophil_m"]
ipsc_array = ["ipsc_cwru1_m", "ipsc_gm23338_53m_gm23248_origin", "ipsc_ips_df_19.11_newborn_m", "ipsc_ips_df_19.7_newborn_m", "ipsc_ips_df_14.7_newborn_m", "ipsc_ips_df_6.9_newborn_m", "ipsc_ips-11a_36m", "ipsc_ips-15b_48f", "ipsc_ips-18a_48f", "ipsc_ips-18c_48f", "ipsc_ips-20b_55m", "ipsc_ips-nihi11_71m_ag20443_origin", "ipsc_ips-nihi7_85f_ag08395_origin", "ipsc_l1-s8", "ipsc_l1-s8r"]
kidney_array = ["kidney_hek293", "kidney_hek293t", "kidney_59_days_59_days_f", "kidney_80_days", "kidney_105_days_f", "kidney_108_days_f", "kidney_113_days_f", "kidney_120_days_f", "kidney_121_days_f", "kidney_76_days_f_76_days_m", "kidney_85_days_f", "kidney_50m", "kidney_67m", "kidney_105_days_m", "kidney_85_days_m", "kidney_87_days_m", "kidney_left_107_days_f", "kidney_left_110_days_f", "kidney_left_147_days_f", "kidney_left_59_days_f_91_days_m", "kidney_left_87_days_f", "kidney_left_89_days_f", "kidney_left_115_days_m", "kidney_left_87_days_m", "kidney_left_96_days_m", "kidney_left_renal_cortex_interstitium_105_days_m", "kidney_left_renal_cortex_interstitium_120_days_m", "kidney_left_renal_pelvis_105_days_m", "kidney_left_renal_pelvis_120_days_m", "kidney_renal_cortex_interstitium_103_days_f", "kidney_renal_cortex_interstitium_120_days_f", "kidney_renal_cortex_interstitium_89_days_f", "kidney_renal_cortex_interstitium_96_days_f", "kidney_renal_cortex_interstitium_108_days_m", "kidney_renal_cortex_interstitium_113_days_m", "kidney_renal_cortex_interstitium_127_days_m", "kidney_renal_cortex_interstitium_91_days_m", "kidney_renal_cortex_interstitium_97_days_m", "kidney_renal_pelvis_103_days_f", "kidney_renal_pelvis_105_days_f", "kidney_renal_pelvis_89_days_f", "kidney_renal_pelvis_96_days_f", "kidney_renal_pelvis_108_days_m", "kidney_renal_pelvis_113_days_m", "kidney_renal_pelvis_127_days_m", "kidney_renal_pelvis_91_days_m", "kidney_renal_pelvis_97_days_m", "kidney_right_107_days_f", "kidney_right_117_days_f", "kidney_right_147_days_f", "kidney_right_87_days_f", "kidney_right_98_days_f", "kidney_right_108_days_m", "kidney_right_115_days_m", "kidney_right_87_days_m", "kidney_right_91_days_m", "kidney_right_96_days_m", "kidney_right_renal_cortex_interstitium_105_days_m", "kidney_right_renal_cortex_interstitium_120_days_m", "kidney_right_renal_pelvis_105_days_m", "kidney_right_renal_pelvis_120_days_m"]
liver_array = ["liver_hepatic_stellate_cell_59f", "liver_hepatocyte", "liver_59_days_80_days", "liver_25f", "liver_101_days_113_days_f", "liver_31m", "liver_32m", "liver_78m", "liver_right_lobe_53f"]
lung_array = ["lung_left_105_days_f", "lung_left_107_days_f", "lung_left_108_days_f", "lung_left_110_days_f", "lung_left_117_days_f", "lung_left_91_days_f", "lung_left_98_days_f", "lung_left_105_days_m", "lung_left_113_days_m", "lung_left_115_days_m", "lung_left_87_days_m", "lung_left_91_days_m", "lung_left_96_days_m", "lung_101_days", "lung_112_days", "lung_67_days", "lung_80_days_76_days_m", "lung_30f", "lung_108_days_f", "lung_120_days_f", "lung_76_days_f", "lung_82_days_f", "lung_85_days_f", "lung_96_days_f", "lung_3m", "lung_103_days_m", "lung_108_days_m", "lung_54_days_58_days_m", "lung_82_days_m", "lung_right_105_days_f", "lung_right_107_days_f", "lung_right_108_days_f", "lung_right_110_days_f", "lung_right_117_days_f", "lung_right_91_days_f", "lung_right_98_days_f", "lung_right_105_days_m", "lung_right_115_days_m", "lung_right_87_days_m", "lung_right_96_days_m", "lung_left_upper_lobe_51f", "lung_left_upper_lobe_53f", "lung_left_upper_lobe_37m"]
lymphoblastoid_array = ["lymphoblastoid_gm06990", "lymphoblastoid_gm08714", "lymphoblastoid_gm10248", "lymphoblastoid_gm10266", "lymphoblastoid_gm12864", "lymphoblastoid_gm12865", "lymphoblastoid_gm12875", "lymphoblastoid_gm12878", "lymphoblastoid_gm12891", "lymphoblastoid_gm12892", "lymphoblastoid_gm18507", "lymphoblastoid_gm19238", "lymphoblastoid_gm19239", "lymphoblastoid_gm19240"]
mesench_array = ["mesench_adipocyte_msc_origin", "mesench_msc_dedifferentiated_amniotic_fluid_origin", "mesench_embryonic_facial_prominence_53_days_58_days", "mesench_msc_adipose_origin"]
muscle_array = ["muscle_cardiac_myocyte", "muscle_forelimb_108_days_f", "muscle_gastrocnemius_medialis_51f", "muscle_gastrocnemius_medialis_53f", "muscle_gastrocnemius_medialis_37m", "muscle_gastrocnemius_medialis_54m", "muscle_leg_hindlimb_120_days_m", "muscle_arm_101_days", "muscle_arm_105_days_f", "muscle_arm_115_days_f", "muscle_arm_120_days_f", "muscle_arm_85_days_f", "muscle_arm_98_days_f", "muscle_arm_101_days_m", "muscle_arm_104_days_m", "muscle_arm_105_days_m", "muscle_arm_113_days_m", "muscle_arm_115_days_m", "muscle_arm_120_days_m", "muscle_arm_127_days_m", "muscle_arm_96_days_m", "muscle_arm_97_days_m", "muscle_back_105_days_f", "muscle_back_113_days_f", "muscle_back_115_days_f", "muscle_back_85_days_f", "muscle_back_98_days_f", "muscle_back_101_days_m", "muscle_back_104_days_m", "muscle_back_105_days_m", "muscle_back_108_days_m", "muscle_back_127_days_m", "muscle_back_91_days_m", "muscle_back_96_days_m", "muscle_back_97_days_m", "muscle_leg_105_days_f", "muscle_leg_110_days_f", "muscle_leg_113_days_f", "muscle_leg_115_days_f", "muscle_leg_85_days_f", "muscle_leg_101_days_m", "muscle_leg_104_days_m", "muscle_leg_105_days_m", "muscle_leg_115_days_m", "muscle_leg_127_days_m", "muscle_leg_96_days_m", "muscle_leg_97_days_m", "muscle_trunk_113_days_f", "muscle_trunk_115_days_f", "muscle_trunk_120_days_f", "muscle_trunk_121_days_f", "muscle_psoas_30f", "muscle_psoas_27m_35m", "muscle_psoas_34m", "muscle_psoas_3m", "muscle_skeletal_cell", "muscle_skeletal_tissue", "muscle_skeletal_tissue_72f", "muscle_skeletal_tissue_54m", "muscle_tongue_59_days_f_76_days_f", "muscle_tongue_72_days_m"]
myosat_array = ["myosat_skeletal_muscle_myoblast_lhcn-m2", "myosat_myocyte_lhcn-m2_origin", "myosat_myotube_skeletal_muscle_myoblast_origin", "myosat_skeletal_muscle_myoblast", "myosat_skeletal_muscle_myoblast_22m", "myosat_skeletal_muscle_satellite_cell_mesoderm_origin_f"]
neurosph_array = ["neurosph_15_weeks_ganglionic_eminence_origin", "neurosph_17_weeks_f_ganglionic_cortex_origin", "neurosph_17_weeks_f_ganglionic_eminence_origin", "neurosph_olfactory_cell_line"]
other_array = ["other_breast_epithelium_51f", "other_breast_epithelium_53f", "other_epidermal_melanocyte", "other_foreskin_melanocyte_newborn_m", "other_limb_embryo_53_days_56_days", "other_limb_embryo_58_days_59_days", "other_mammary_stem_cell"]
pancreas_array = ["pancreas_body_51f", "pancreas_body_53f", "pancreas_body_37m", "pancreas_body_54m", "pancreas_islet_precursor_cell", "pancreas_30f", "pancreas_34m"]
placenta_and_emm_array = ["placenta_and_eem_amnion_16_weeks_m", "placenta_and_eem_amnion_stem_cell", "placenta_and_eem_chorion", "placenta_and_eem_chorion_40_weeks_f", "placenta_and_eem_chorion_16_weeks_m", "placenta_and_eem_chorionic_villus_16_weeks", "placenta_and_eem_chorionic_villus_40_weeks_f", "placenta_and_eem_chorionic_villus_16_weeks_m", "placenta_and_eem_chorionic_villus_38_weeks_m", "placenta_and_eem_trophoblast_htr-8_svneo", "placenta_and_eem_placenta_102_days", "placenta_and_eem_placenta_16_weeks", "placenta_and_eem_placenta_53_days", "placenta_and_eem_placenta_56_days_59_days", "placenta_and_eem_placenta_101_days_f_105_days_m", "placenta_and_eem_placenta_105_days_f", "placenta_and_eem_placenta_108_days_f", "placenta_and_eem_placenta_113_days_f", "placenta_and_eem_placenta_85_days_f", "placenta_and_eem_placenta_16_weeks_m", "placenta_and_eem_placenta_85_days_m", "placenta_and_eem_placenta_91_days_m", "placenta_and_eem_placental_basal_plate_40_weeks_f", "placenta_and_eem_placental_basal_plate_38_weeks_m", "placenta_and_eem_trophoblast_cell_17_weeks_18_weeks", "placenta_and_eem_trophoblast_cell_21_weeks", "placenta_and_eem_trophoblast_cell_23_weeks", "placenta_and_eem_trophoblast_cell_39_weeks_40_weeks", "placenta_and_eem_trophoblast_20_weeks_f", "placenta_and_eem_trophoblast_40_weeks_f", "placenta_and_eem_umbilical_cord_59_days_76_days_m"]
pns_array = ["pns_spinal_cord_108_days_f", "pns_spinal_cord_113_days_f", "pns_spinal_cord_59_days_f_72_days_m", "pns_spinal_cord_87_days_f", "pns_spinal_cord_89_days_f", "pns_spinal_cord_105_days_m", "pns_spinal_cord_96_days_m", "pns_tibial_nerve_51f", "pns_tibial_nerve_53f", "pns_tibial_nerve_37m"]
reproductive_array = ["reproductive_prostate_gland_37m", "reproductive_prostate_gland_54m", "reproductive_uterus_53f", "reproductive_vagina_51f", "reproductive_vagina_53f"]
sm_muscle_array = ["sm_muscle_colon_56f", "sm_muscle_colon_77f", "sm_muscle_duodenum_59m", "sm_muscle_duodenum_73m", "sm_muscle_rectum_50f", "sm_muscle_brain_vasculature_smooth_cell", "sm_muscle_stomach_84f", "sm_muscle_stomach_59m"]
spleen_array = ["spleen_112_days", "spleen_30f", "spleen_53f", "spleen_34m", "spleen_54m", "spleen_3m"]
stromal_array = ["stromal_skin_fibroblast_ag04449", "stromal_lung_fibroblast_ag04450", "stromal_skin_fibroblast_ag08395", "stromal_lung_fibroblast_ag08396", "stromal_skin_fibroblast_ag08396", "stromal_gingival_fibroblast_ag09319", "stromal_skin_fibroblast_ag10803", "stromal_skin_fibroblast_ag20443", "stromal_skin_fibroblast_bj", "stromal_brain_pericyte", "stromal_cardiac_fibroblast", "stromal_cardiac_fibroblast_f", "stromal_cardiac_fibroblast_94_days_f_98_days_f", "stromal_skin_fibroblast_eh", "stromal_skin_fibroblast_el", "stromal_skin_fibroblast_elr", "stromal_breast_fibroblast_17f", "stromal_breast_fibroblast_26f", "stromal_dermis_fibroblast", "stromal_dermis_fibroblast_f", "stromal_dermis_fibroblast_none_f", "stromal_gingival_fibroblast", "stromal_lung_fibroblast", "stromal_lung_fibroblast_11f_45m", "stromal_lung_fibroblast_45m", "stromal_mammary_fibroblast_f", "stromal_peridontal_ligament_fibroblast_m", "stromal_pulmonary_artery_fibroblast", "stromal_skin_fibroblast_abdomen_97_days_m", "stromal_aorta_fibroblast_f", "stromal_conjunctiva_fibroblast", "stromal_villous_mesenchyme_fibroblast", "stromal_foreskin_fibroblast_newborn_m", "stromal_skin_fibroblast_gm03348", "stromal_skin_fibroblast_gm04503", "stromal_skin_fibroblast_gm04504", "stromal_skin_fibroblast_gm23248", "stromal_foreskin_fibroblast_hff-myc_foreskin_fibroblast_origin", "stromal_lung_fibroblast_imr-90", "stromal_lung_fibroblast_97_days_m", "stromal_bone_marrow_m", "stromal_lung_fibroblast_wi38"]
thymus_array = ["thymus_embryo_f", "thymus_105_days_f", "thymus_110_days_f", "thymus_113_days_f", "thymus_147_days_f", "thymus_98_days_f", "thymus_3m", "thymus_104_days_m", "thymus_108_days_m", "thymus_113_days_m", "thymus_127_days_m"]
urinary_array = ["urinary_bladder_34m", "urinary_bladder_76_days_m", "urinary_urothelium_cell_line"]

all_array = adipose_array + blood_array + bone_array + brain_array + cancer_array + digestive_array + endocrine_array + endothelial_array + epithelial_array + esc_deriv_array + esc_array + eye_array + heart_array + hsc_and_b_cell_array + ipsc_array + kidney_array + liver_array + lung_array + lymphoblastoid_array + mesench_array + muscle_array + myosat_array + neurosph_array + other_array + pancreas_array + placenta_and_emm_array + pns_array + reproductive_array + sm_muscle_array + spleen_array + stromal_array + thymus_array + urinary_array
tissue_array = list(chain.from_iterable(all_array if item == "all" else adipose_array if item == "adipose" else blood_array if item == "blood" else bone_array if item == "bone" else brain_array if item == "brain" else cancer_array if item == "cancer" else digestive_array if item == "digestive" else endocrine_array if item == "endocrine" else endothelial_array if item == "endothelial" else epithelial_array if item == "epithelial" else esc_deriv_array if item == "esc_deriv" else esc_array if item == "esc" else eye_array if item == "eye" else heart_array if item == "heart" else hsc_and_b_cell_array if item == "hsc_and_b_cell" else ipsc_array if item == "ipsc" else kidney_array if item == "kidney" else liver_array if item == "liver" else lung_array if item == "lung" else lymphoblastoid_array if item == "lymphoblastoid" else mesench_array if item == "mesench" else muscle_array if item == "muscle" else myosat_array if item == "myosat" else neurosph_array if item == "neurosph" else other_array if item == "other" else pancreas_array if item == "pancreas" else placenta_and_emm_array if item == "placenta_and_emm" else pns_array if item == "pns" else reproductive_array if item == "reproductive" else sm_muscle_array if item == "sm_muscle" else spleen_array if item == "spleen" else stromal_array if item == "stromal" else thymus_array if item == "thymus" else urinary_array if item == "urinary" else [item] for item in tissue_array))

###Create an array with all of the parsed in epigenetic arguments
if ("user_provided_files" in tissue_array or "user_provided_urls" in tissue_array):
	user_4me1_array = args.user_4me1.split(",")
	user_4me3_array = args.user_4me3.split(",")
	user_27ac_array = args.user_27ac.split(",")
	user_27me3_array = args.user_27me3.split(",")
	#user_36me3_array = args.user_36me3.split(",")
	user_ctcf_array = args.user_ctcf.split(",")
	user_dnase_array = args.user_dnase.split(",")
	user_tissue_names_array = args.user_tissue_names.lower().split(",")
	user_files_index = 0
	
elif args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39" or args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38" or args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "grcm37":
	assert args.tissue.lower() == "user_provided_files" or args.tissue.lower() == "user_provided_urls", "Tissue type must be one of the 2 valid mouse options (User_provided_files or User_provided_urls)"

#Download the appropriate reference files based on the specified genome build and tissue arguments
if not os.path.exists(args.ref_dir):
	os.mkdir(args.ref_dir)

def download_ref_genome(url_ref, out_name):
	out_path = os.path.join(args.ref_dir, out_name)
	r = requests.get(url_ref, allow_redirects = True)
	open(out_path, 'wb').write(r.content)

def download_ref(epimap_accession, out_name):
	#ATAC-seq
	if not os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " ATAC-seq...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_ATAC-seq_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_ATAC-seq.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_ATAC-seq_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_ATAC-seq_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_ATAC-seq_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_ATAC-seq_hg19" + ".bigWig")
	#EP300
	if not os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " EP300...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_EP300_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_EP300_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_EP300.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_EP300_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_EP300_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_EP300_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_EP300_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_EP300_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_EP300_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_EP300_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_EP300_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_EP300_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_EP300_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_EP300_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_EP300_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_EP300_hg19" + ".bigWig")
	#H2AFZ
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H2AFZ...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H2AFZ_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H2AFZ.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H2AFZ_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H2AFZ_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H2AFZ_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H2AFZ_hg19" + ".bigWig")
	#H3K4me1
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K4me1...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K4me1_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K4me1.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me1_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K4me1_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me1_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K4me1_hg19" + ".bigWig")
	#H3K4me2
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K4me2...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K4me2_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K4me2.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me2_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K4me2_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me2_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K4me2_hg19" + ".bigWig")
	#H3K4me3
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K4me3...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K4me3_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K4me3.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K4me3_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K4me3_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K4me3_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K4me3_hg19" + ".bigWig")
	#H3K9ac
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K9ac...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K9ac_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K9ac.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K9ac_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K9ac_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K9ac_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K9ac_hg19" + ".bigWig")
	#H3K9me3
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K9me3...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K9me3_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K9me3.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K9me3_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K9me3_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K9me3_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K9me3_hg19" + ".bigWig")				
	#H3K27ac
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K27ac...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K27ac_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K27ac.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K27ac_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K27ac_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K27ac_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K27ac_hg19" + ".bigWig")	
	#H3K27me3
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K27me3...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K27me3_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K27me3.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K27me3_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K27me3_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K27me3_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K27me3_hg19" + ".bigWig")	
	#DNase-seq
	if not os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " DNase-seq...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_DNase-seq_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_DNase-seq.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_DNase-seq_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_DNase-seq_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_DNase-seq_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_DNase-seq_hg19" + ".bigWig")	
	#CTCF
	if not os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " CTCF...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_CTCF_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_CTCF_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_CTCF.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_CTCF_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_CTCF_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_CTCF_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_CTCF_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_CTCF_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_CTCF_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_CTCF_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_CTCF_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_CTCF_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_CTCF_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_CTCF_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_CTCF_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_CTCF_hg19" + ".bigWig")
	#H3K36me3
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K36me3...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K36me3_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K36me3.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K36me3_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K36me3_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K36me3_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K36me3_hg19" + ".bigWig")
	#H3K79me2
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H3K79me2...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H3K79me2_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H3K79me2.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H3K79me2_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H3K79me2_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H3K79me2_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H3K79me2_hg19" + ".bigWig")
	#H4K20me1
	if not os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " H4K20me1...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_H4K20me1_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_H4K20me1.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_H4K20me1_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_H4K20me1_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_H4K20me1_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_H4K20me1_hg19" + ".bigWig")
	#POLR2A
	if not os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " POL2RA...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_POL2RA_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_POL2RA.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_POL2RA_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_POL2RA_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_POL2RA_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_POL2RA_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_POL2RA_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_POL2RA_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_POL2RA_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_POL2RA_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_POL2RA_hg19" + ".bigWig")
	#RAD21
	if not os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " RAD21...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_RAD21_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_RAD21_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_RAD21.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_RAD21_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_RAD21_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_RAD21_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_RAD21_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_RAD21_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_RAD21_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_RAD21_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_RAD21_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_RAD21_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_RAD21_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_RAD21_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_RAD21_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_RAD21_hg19" + ".bigWig")
	#SMC3
	if not os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19" + ".bed.gz"):
		if args.verbose:
			print("Downloading bigWig and compiling BED for " + out_name + " SMC3...")
		#Observed
		response_observed = requests.get("https://epigenome.wustl.edu/epimap/data/observed/FINAL_SMC3_" + epimap_accession + ".sub_VS_Uniform_BKG_CONTROL_36_50000000.pval.signal.bedgraph.gz.bigWig", allow_redirects = True)
		if response_observed.status_code == 200:
			out_path_observed = os.path.join(args.ref_dir, out_name + "_SMC3_hg19_observed" + ".bigWig")
			open(out_path_observed, 'wb').write(response_observed.content)
    	#Imputed
		response_imputed = requests.get("https://epigenome.wustl.edu/epimap/data/imputed/impute_" + epimap_accession + "_SMC3.bigWig", allow_redirects = True)
		if response_imputed.status_code == 200:
			out_path_imputed = os.path.join(args.ref_dir, out_name + "_SMC3_hg19_imputed" + ".bigWig")
			open(out_path_imputed, 'wb').write(response_imputed.content)
		#If both imputed and observed, merge files together into a bedGraph before peak calling. Otherwise, rename whatever invididual files have been donloaded and carry out the rest of the peak calling workflow
		if os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_observed" + ".bigWig") and os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_imputed" + ".bigWig"):
			subprocess.run(["bigWigMerge", os.path.join(args.ref_dir, out_name + "_SMC3_hg19_observed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_SMC3_hg19_imputed" + ".bigWig"), os.path.join(args.ref_dir, out_name + "_SMC3_hg19" + ".bedGraph")])
			subprocess.run(["bedGraphToBigWig", os.path.join(args.ref_dir, out_name + "_SMC3_hg19" + ".bedGraph"), os.path.join("/home/bettimj/reference_genomes", args.ref_genome + ".chrom.sizes"), os.path.join(args.ref_dir, out_name + "_SMC3_hg19" + ".bigWig")])
			os.remove(os.path.join(args.ref_dir, out_name + "_SMC3_hg19_observed" + ".bigWig"))
			os.remove(os.path.join(args.ref_dir, out_name + "_SMC3_hg19_imputed" + ".bigWig"))
		elif os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_observed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_imputed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_SMC3_hg19_observed" + ".bigWig", args.ref_dir + out_name + "_SMC3_hg19" + ".bigWig")
		elif os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_imputed" + ".bigWig") and not os.path.exists(args.ref_dir + "/" + out_name + "_SMC3_hg19_observed" + ".bigWig"):
			os.rename(args.ref_dir + "/" + out_name + "_SMC3_hg19_imputed" + ".bigWig", args.ref_dir + "/" + out_name + "_SMC3_hg19" + ".bigWig")
		
if not os.path.exists(args.ref_dir):
	os.mkdir(args.ref_dir)

if args.verbose:
	print("Downloading {ref_genome} TSS coordinates and {tissue_type} histone ChIP-seq bed files...".format(ref_genome = args.ref_genome, tissue_type = args.tissue))
	print("\n")

###Human###
#GRCh38/hg38
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	if not os.path.exists(args.ref_dir + "/refTSS_v3.1_human_coordinate.hg38.bed.gz"):
		download_ref_genome("http://reftss.clst.riken.jp/datafiles/3.1/human/refTSS_v3.1_human_coordinate.hg38.bed.gz", "refTSS_v3.1_human_coordinate.hg38.bed.gz")
	
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_hg38.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_hg38.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_hg38.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_hg38.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_hg38.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_hg38.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_hg38.bed.gz"))
			user_files_index += 1
			

###################################################################################
#GRCh37/hg19
if args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	#This TSS dataset was manually lifted over from the original hg38 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		###Adipose
		#Adipocyte (adipocyte)
		if tissue == "adipose_adipocyte":
			download_ref("BSS00038", tissue)
		#Adipose tissue (adipose tissue male adult (34 years))
		if tissue == "adipose_34m":
			download_ref("BSS00043", tissue)
		#Omental fat pad (omental fat pad female adult (51 year))
		if tissue == "adipose_omental_fat_pad_51f":
			download_ref("BSS01393", tissue)
		#Omental fat pad (omental fat pad female adult (53 year))
		if tissue == "adipose_omental_fat_pad_53f":
			download_ref("BSS01394", tissue)
		#Adipose (subcutaneous abdominal adipose tissue nuclear fraction female adult (25 years))
		if tissue == "adipose_subcutaneous_25f":
			download_ref("BSS01665", tissue)
		#Adipose (subcutaneous abdominal adipose tissue nuclear fraction female adult (41 year))
		if tissue == "adipose_subcutaneous_41f":
			download_ref("BSS01666", tissue)
		#Adipose (subcutaneous abdominal adipose tissue nuclear fraction female adult (49 years))
		if tissue == "adipose_subcutaneous_49f":
			download_ref("BSS01667", tissue)
		#Adipose (subcutaneous abdominal adipose tissue nuclear fraction female adult (59 years))
		if tissue == "adipose_subcutaneous_59f":
			download_ref("BSS01668", tissue)
		#Adipose (subcutaneous abdominal adipose tissue nuclear fraction female adult (81 year))
		if tissue == "adipose_subcutaneous_81f":
			download_ref("BSS01669", tissue)
			
		###Blood
		#CD4 T-cell (CD4-positive, alpha-beta memory T cell)
		if tissue == "blood_cd4_alpha_beta_memory_t_cell":
			download_ref("BSS00183", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta memory T cell originated from blood cell)
		if tissue == "blood_cd4_alpha_beta_memory_t_cell_blood_origin":
			download_ref("BSS00185", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta T cell)
		if tissue == "blood_cd4_alpha_beta_t_cell":
			download_ref("BSS00186", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta T cell female adult (33 years))
		if tissue == "blood_cd4_alpha_beta_t_cell_33f":
			download_ref("BSS00188", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta T cell male adult (21 year))
		if tissue == "blood_cd4_alpha_beta_t_cell_21m":
			download_ref("BSS00189", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta T cell male adult (37 years))
		if tissue == "blood_cd4_alpha_beta_t_cell_37m":
			download_ref("BSS00190", tissue)
		#CD4 T-cell (CD4-positive, alpha-beta T cell treated with phorbol 13-acetate 12-myristate , ionomycin)
		if tissue == "blood_cd4_alpha_beta_t_cell_phorbol":
			download_ref("BSS00191", tissue)
		#CD4 T-cell (CD4-positive, CD25-positive, alpha-beta regulatory T cell)
		if tissue == "blood_cd4_cd25_alpha_beta_t_cell":
			download_ref("BSS00192", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta memory T cell)
		if tissue == "blood_cd8_alpha_beta_memory_t_cell":
			download_ref("BSS00193", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell)
		if tissue == "blood_cd8_alpha_beta_t_cell":
			download_ref("BSS00194", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell female adult (33 years))
		if tissue == "blood_cd8_alpha_beta_t_cell_33f":
			download_ref("BSS00195", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell female adult (34 years))
		if tissue == "blood_cd8_alpha_beta_t_cell_34f":
			download_ref("BSS00196", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell male adult (21 year))
		if tissue == "blood_cd8_alpha_beta_t_cell_21m":
			download_ref("BSS00197", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell male adult (28 years))
		if tissue == "blood_cd8_alpha_beta_t_cell_28m":
			download_ref("BSS00198", tissue)
		#CD8 T-cell (CD8-positive, alpha-beta T cell male adult (37 years))
		if tissue == "blood_cd8_alpha_beta_t_cell_37m":
			download_ref("BSS00200", tissue)
		#CD4 T-cell (effector memory CD4-positive, alpha-beta T cell)
		if tissue == "blood_effector_memory_cd4_alpha_beta_t_cell":
			download_ref("BSS00274", tissue)
		#Mononuclear cell (mononuclear cell male)
		if tissue == "blood_mononuclear_cell_male":
			download_ref("BSS01279", tissue)
		#Naive T-cell (naive thymus-derived CD4-positive, alpha-beta T cell)
		if tissue == "blood_naive_thymus_cd4_alpha_beta_t_cell":
			download_ref("BSS01346", tissue)
		#Naive T-cell (naive thymus-derived CD4-positive, alpha-beta T cell female adult (35 years))
		if tissue == "blood_naive_thymus_cd4_alpha_beta_t_cell_35f":
			download_ref("BSS01347", tissue)
		#Naive T-cell (naive thymus-derived CD4-positive, alpha-beta T cell male adult (26 years))
		if tissue == "blood_naive_thymus_cd4_alpha_beta_t_cell_26m":
			download_ref("BSS01348", tissue)
		#Mononuclear cell (peripheral blood mononuclear cell female adult (28 years))
		if tissue == "blood_peripheral_mononuclear_cell_28f":
			download_ref("BSS01419", tissue)
		#Mononuclear cell (peripheral blood mononuclear cell male adult (27 years))
		if tissue == "blood_peripheral_mononuclear_cell_27m":
			download_ref("BSS01420", tissue)
		#Mononuclear cell (peripheral blood mononuclear cell male adult (28 years))
		if tissue == "blood_peripheral_mononuclear_cell_28m":
			download_ref("BSS01421", tissue)
		#Mononuclear cell (peripheral blood mononuclear cell male adult (32 years))
		if tissue == "blood_peripheral_mononuclear_cell_32m":
			download_ref("BSS01423", tissue)
		#Mononuclear cell (peripheral blood mononuclear cell male adult (39 years))
		if tissue == "blood_peripheral_mononuclear_cell_39m":
			download_ref("BSS01424", tissue)
		#Regulatory T-cell (regulatory T cell female adult (35 years))
		if tissue == "blood_regulatory_t_cell_35f":
			download_ref("BSS01478", tissue)
		#Regulatory T-cell (regulatory T cell male adult (28 years))
		if tissue == "blood_regulatory_t_cell_28m":
			download_ref("BSS01479", tissue)
		#Regulatory T-cell (regulatory T cell originated from blood cell)
		if tissue == "blood_regulatory_t_cell_blood_origin":
			download_ref("BSS01480", tissue)
		#T-cell (T-cell)
		if tissue == "blood_t_cell":
			download_ref("BSS01684", tissue)
		#T-cell (T-cell male adult (21 year))
		if tissue == "blood_t_cell_21m":
			download_ref("BSS01687", tissue)
		#T-cell (T-cell male adult (36 years))
		if tissue == "blood_t_cell_36m":
			download_ref("BSS01688", tissue)
		#T-cell (T-cell male adult (37 years))
		if tissue == "blood_t_cell_37m":
			download_ref("BSS01689", tissue)
		#T1 helper cell (T-helper 1 cell)
		if tissue == "blood_t1_helper_cell":
			download_ref("BSS01690", tissue)
		#T1 helper cell (T-helper 1 cell female adult (26 years))
		if tissue == "blood_t1_helper_cell_26f":
			download_ref("BSS01691", tissue)
		#T1 helper cell (T-helper 1 cell male adult (33 years))
		if tissue == "blood_t1_helper_cell_33m":
			download_ref("BSS01692", tissue)
		#T17 helper cell (T-helper 17 cell)
		if tissue == "blood_t17_helper_cell":
			download_ref("BSS01693", tissue)
		#T17 helper cell (T-helper 17 cell originated from blood cell)
		if tissue == "blood_t17_helper_cell_blood_origin":
			download_ref("BSS01694", tissue)
		#T17 helper cell (T-helper 17 cell treated with phorbol 13-acetate 12-myristate , ionomycin)
		if tissue == "blood_t17_helper_cell_phorbol":
			download_ref("BSS01695", tissue)
		#T2 helper cell (T-helper 2 cell female adult (26 years))
		if tissue == "blood_t2_helper_cell_26f":
			download_ref("BSS01697", tissue)
		#T2 helper cell (T-helper 2 cell male adult (33 years))
		if tissue == "blood_t2_helper_cell_33m":
			download_ref("BSS01698", tissue)
		
		###Bone
		#Bone arm (arm bone male embryo (81 day))
		if tissue == "bone_arm":
			download_ref("BSS00084", tissue)
		#Bone femur (femur female embryo (98 days))
		if tissue == "bone_femur":
			download_ref("BSS00330", tissue)
		#Bone marrow stroma (HS-5)
		if tissue == "bone_marrow_stroma":
			download_ref("BSS00705", tissue)
		#Bone leg (leg bone male embryo (81 day))
		if tissue == "bone_leg":
			download_ref("BSS01154", tissue)
		#Osteoblast (osteoblast - primary cell)
		if tissue == "bone_osteoblast":
			download_ref("BSS01397", tissue)
		
		##Brain
		#Ammons horn (Ammon s horn male adult (84 years))
		if tissue == "brain_ammons_horn_84m":
			download_ref("BSS00071", tissue)
		#Angular gyrus (angular gyrus female adult (75 years))
		if tissue == "brain_angular_gyrus_75f":
			download_ref("BSS00077", tissue)
		#Angular gyrus (angular gyrus male adult (81 year))
		if tissue == "brain_angular_gyrus_81m":
			download_ref("BSS00078", tissue)
		#Astrocyte (astrocyte)
		if tissue == "brain_astrocyte":
			download_ref("BSS00089", tissue)
		#Astrocyte cerebellum (astrocyte of the cerebellum)
		if tissue == "brain_astrocyte_cerebellum":
			download_ref("BSS00090", tissue)
		#Astrocyte hippocampus (astrocyte of the hippocampus)
		if tissue == "brain_astrocyte_hippocampus":
			download_ref("BSS00091", tissue)
		#Astrocyte spinal cord (astrocyte of the spinal cord)
		if tissue == "brain_astrocyte_spinal_cord":
			download_ref("BSS00092", tissue)
		#Brain (brain embryo (112 days))
		if tissue == "brain_embryo_112_days":
			download_ref("BSS00125", tissue)
		#Brain (brain embryo (56 days) and male embryo (58 days))
		if tissue == "brain_embryo_56_58_days":
			download_ref("BSS00126", tissue)
		#Brain (brain embryo (80 days))
		if tissue == "brain_embryo_80_days":
			download_ref("BSS00127", tissue)
		#Brain (brain female embryo (105 days))
		if tissue == "brain_embryo_105_days_f":
			download_ref("BSS00127", tissue)
		#Brain (brain female embryo (109 days))
		if tissue == "brain_embryo_109_days_f":
			download_ref("BSS00130", tissue)
		#Brain (brain female embryo (117 days))
		if tissue == "brain_embryo_117_days_f":
			download_ref("BSS00131", tissue)
		#Brain (brain female embryo (120 days))
		if tissue == "brain_embryo_120_days_f":
			download_ref("BSS00132", tissue)
		#Brain (brain female embryo (142 days))
		if tissue == "brain_embryo_142_days_f":
			download_ref("BSS00133", tissue)
		#Brain (brain female embryo (17 weeks))
		if tissue == "brain_embryo_17_weeks_f":
			download_ref("BSS00134", tissue)
		#Brain (brain female embryo (85 days))
		if tissue == "brain_embryo_85_days_f":
			download_ref("BSS00135", tissue)
		#Brain (brain female embryo (96 days))
		if tissue == "brain_embryo_96_days_f":
			download_ref("BSS00136", tissue)
		#Brain (brain male embryo (101 day))
		if tissue == "brain_embryo_101_days_m":
			download_ref("BSS00138", tissue)
		#Brain (brain male embryo (104 days))
		if tissue == "brain_embryo_104_days_m":
			download_ref("BSS00139", tissue)
		#Brain (brain male embryo (105 days))
		if tissue == "brain_embryo_105_days_m":
			download_ref("BSS00140", tissue)
		#Brain (brain male embryo (122 days))
		if tissue == "brain_embryo_122_days_m":
			download_ref("BSS00141", tissue)
		#Brain (brain male embryo (72 days) and male embryo (76 days))
		if tissue == "brain_embryo_72_76_days_m":
			download_ref("BSS00142", tissue)
		#Brain caudate nucleus (caudate nucleus female adult (75 years))
		if tissue == "brain_caudate_nucleus_75f":
			download_ref("BSS00173", tissue)
		#Brain caudate nucleus (caudate nucleus male adult (78 years))
		if tissue == "brain_caudate_nucleus_78m":
			download_ref("BSS00174", tissue)
		#Brain caudate nucleus (caudate nucleus male adult (81 year))
		if tissue == "brain_caudate_nucleus_81m":
			download_ref("BSS00175", tissue)
		#Brain cerebellar cortex (cerebellar cortex male adult (78 years) and male adult (84 years))
		if tissue == "brain_cerebellar_cortex_78_81m":
			download_ref("BSS00201", tissue)
		#Brain cerebellum (cerebellum male adult (27 years) and male adult (35 years))
		if tissue == "brain_cerebellum_27_35m":
			download_ref("BSS00206", tissue)
		#Brain cerebellum (cerebellum male adult (53 years))
		if tissue == "brain_cerebellum_53m":
			download_ref("BSS00207", tissue)
		#Brain cingulate gyrus (cingulate gyrus female adult (75 years))
		if tissue == "brain_cingulate_gyrus_75f":
			download_ref("BSS00219", tissue)
		#Brain cingulate gyrus (cingulate gyrus male adult (81 year))
		if tissue == "brain_cingulate_gyrus_81m":
			download_ref("BSS00220", tissue)
		#Brain frontal cortex (frontal cortex female adult (67 years) and female adult (80 years))
		if tissue == "brain_frontal_cortex_67_80f":
			download_ref("BSS00369", tissue)
		#Brain frontal cortex (frontal cortex male adult (27 years) and male adult (35 years))
		if tissue == "brain_frontal_cortex_27_35m":
			download_ref("BSS00369", tissue)
		#Brain germinal matrix (germinal matrix male embryo (20 weeks))
		if tissue == "brain_germinal_matrix_20_weeks_m":
			download_ref("BSS00385", tissue)
		#Brain globus pallidus (globus pallidus male adult (78 years) and male adult (84 years))
		if tissue == "brain_globus_pallidus_78_84m":
			download_ref("BSS00386", tissue)
		#Brain inferior parietal cortex (inferior parietal cortex male adult (84 years))
		if tissue == "brain_inferior_parietal_cortex_84m":
			download_ref("BSS00729", tissue)
		#Brain hippocampus (layer of hippocampus female adult (75 years))
		if tissue == "brain_hippocampus_75f":
			download_ref("BSS01124", tissue)
		#Brain hippocampus (layer of hippocampus male adult (73 years))
		if tissue == "brain_hippocampus_73m":
			download_ref("BSS01124", tissue)
		#Brain medulla oblongata (medulla oblongata male adult (78 years) and male adult (84 years))
		if tissue == "brain_medulla_oblongata_78_84m":
			download_ref("BSS01250", tissue)
		#Brain midbrain (midbrain male adult (78 years) and male adult (84 years))
		if tissue == "brain_midbrain_78_84m":
			download_ref("BSS01270", tissue)
		#Brain middle frontal area (middle frontal area 46 female adult (75 years))
		if tissue == "brain_middle_frontal_area_75f":
			download_ref("BSS01271", tissue)
		#Brain middle frontal area (middle frontal area 46 male adult (81 year))
		if tissue == "brain_middle_frontal_area_81m":
			download_ref("BSS01272", tissue)
		#Brain middle frontal gyrus (middle frontal gyrus male adult (78 years))
		if tissue == "brain_middle_frontal_gyrus_78m":
			download_ref("BSS01273", tissue)
		#Brain occipital lobe (occipital lobe male adult (84 years))
		if tissue == "brain_occipital_lobe_84m":
			download_ref("BSS01388", tissue)
		#Brain pons (pons male adult (78 years))
		if tissue == "brain_pons_78m":
			download_ref("BSS01451", tissue)
		#Brain putamen (putamen male adult (78 years))
		if tissue == "brain_putamen_78m":
			download_ref("BSS01469", tissue)
		#Brain substantia nigra (substantia nigra female adult (75 years))
		if tissue == "brain_substantia_nigra_75f":
			download_ref("BSS01675", tissue)
		#Brain substantia nigra (substantia nigra male adult (81 year))
		if tissue == "brain_substantia_nigra_81m":
			download_ref("BSS01676", tissue)
		#Brain superior temporal gyrus (superior temporal gyrus male adult (84 years))
		if tissue == "brain_superior_temporal_gyrus_84m":
			download_ref("BSS01677", tissue)
		#Brain temporal lobe (temporal lobe female adult (75 years))
		if tissue == "brain_temporal_lobe_75f":
			download_ref("BSS01712", tissue)
		#Brain temporal lobe (temporal lobe male adult (81 year))
		if tissue == "brain_temporal_lobe_81m":
			download_ref("BSS01714", tissue)
				
		##Cancer
		#Prostate epithelial carcinoma (22Rv1)
		if tissue == "cancer_prostate_epithelial_22rv1":
			download_ref("BSS00001", tissue)
		#Pancreas adenocarcinoma (8988T)
		if tissue == "cancer_pancreas_adenocarcinoma_8988t":
			download_ref("BSS00003", tissue)
		#Glioblastoma (A172)
		if tissue == "cancer_glioblastoma_a172":
			download_ref("BSS00004", tissue)
		#Lung epithelial carcinoma (A549)
		if tissue == "cancer_lung_epithelial_carcinoma_a549":
			download_ref("BSS00007", tissue)
		#Lung epithelial carcinoma (A549 treated with 0.02pct ethanol for 1 hour)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_ethanol_1hr":
			download_ref("BSS00013", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 1 hour)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_1hr":
			download_ref("BSS00015", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 10 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10hr":
			download_ref("BSS00016", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 10 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_10min":
			download_ref("BSS00017", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 12 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_12hr":
			download_ref("BSS00018", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 15 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_15min":
			download_ref("BSS00019", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 2 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_2hr":
			download_ref("BSS00020", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 20 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_20min":
			download_ref("BSS00021", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 25 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_25min":
			download_ref("BSS00022", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 3 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_3hr":
			download_ref("BSS00023", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 30 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_30min":
			download_ref("BSS00024", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 4 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_4hr":
			download_ref("BSS00025", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 5 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5hr":
			download_ref("BSS00026", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 5 minutes)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_5min":
			download_ref("BSS00027", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 6 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_6hr":
			download_ref("BSS00028", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 7 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_7hr":
			download_ref("BSS00029", tissue)
		#Lung epithelial carcinoma (A549 treated with 100 nM dexamethasone for 8 hours)
		if tissue == "cancer_lung_epithelial_carcinoma_a549_treated_dexamethasone_8hr":
			download_ref("BSS00030", tissue)
		#Muscle Ewing sarcoma (A673)
		if tissue == "cancer_muscle_ewing_sarcoma_a673":
			download_ref("BSS00035", tissue)
		#Adenoid cystic sarcoma (ACC112)
		if tissue == "cancer_adenoid_cystic_carcinoma_acc112":
			download_ref("BSS00036", tissue)
		#Renal cell carcinoma (ACHN)
		if tissue == "cancer_renal_cell_carcinoma_achn":
			download_ref("BSS00037", tissue)
		#Neuroblastoma (BE2C)
		if tissue == "cancer_neuroblastoma_be2c":
			download_ref("BSS00102", tissue)
		#Prostate cancer (C4-2B)
		if tissue == "cancer_prostate_c4-2b":
			download_ref("BSS00157", tissue)
		#Colorectal adenocarcinoma (Caco-2)
		if tissue == "cancer_colorectal_adenocarcinoma_caco-2":
			download_ref("BSS00159", tissue)
		#Kidney clear cell carcinoma (Caki2)
		if tissue == "cancer_kidney_clear_cell_carcinoma_caki2":
			download_ref("BSS00160", tissue)
		#Myelogenous leukemia (CMK)
		if tissue == "cancer_myelogenous_leukemia_cmk":
			download_ref("BSS00221", tissue)
		#Melanoma (COLO829)
		if tissue == "cancer_melanoma_colo829":
			download_ref("BSS00222", tissue)
		#Desmoplastic medulloblastoma (Daoy)
		if tissue == "cancer_desmoplastic_medulloblastoma_daoy":
			download_ref("BSS00246", tissue)
		#Acute lymphoblastic leukemia (DND-41)
		if tissue == "cancer_acute_lymphoblastic_leukemia_dnd-41":
			download_ref("BSS00267", tissue)
		#B cell lymphoma (DOHH2)
		if tissue == "cancer_b_cell_lymphoma_dohh2":
			download_ref("BSS00268", tissue)
		#Kidney rhaboid tumor (G401)
		if tissue == "cancer_kidney_rhaboid_tumor_g401":
			download_ref("BSS00372", tissue)
		#Neuroglioma (H4)
		if tissue == "cancer_neuroglioma_h4":
			download_ref("BSS00481", tissue)
		#Glioblastoma (H54)
		if tissue == "cancer_glioblastoma_h54":
			download_ref("BSS00481", tissue)
		#Haploid myelogenous leukemia (HAP-1)
		if tissue == "cancer_haploid_myelogenous_leukemia_hap-1":
			download_ref("BSS00491", tissue)
		#Colorectal adenocarcinoma (HCT116)
		if tissue == "cancer_colorectal_adenocarcinoma_hct116":
			download_ref("BSS00492", tissue)
		#Cervix adenocarcinoma (HeLa-S3)
		if tissue == "cancer_cervix_adenocarcinoma_hela-s3":
			download_ref("BSS00529", tissue)
		#Cervix adenocarcinoma (HeLa-S3 G1b phase)
		if tissue == "cancer_cervix_adenocarcinoma_hela-s3_g1b_phase":
			download_ref("BSS00531", tissue)
		#Cervix adenocarcinoma (HeLa-S3 treated with interferon alpha for 4 hours)
		if tissue == "cancer_cervix_adenocarcinoma_hela-s3_treated_interferon_alpha_4hr":
			download_ref("BSS00531", tissue)
		#Hepatocellular carcinoma (HepG2)
		if tissue == "cancer_hepatocellular_carcinoma_hepg2":
			download_ref("BSS00558", tissue)
		#Acute promyelocytic leukemia (HL-60)
		if tissue == "cancer_acute_promyelocytic_leukemia_hl-60":
			download_ref("BSS00702", tissue)
		#Colorectal adenocarcinoma (HT-29)
		if tissue == "cancer_colorectal_adenocarcinoma_ht-29":
			download_ref("BSS00708", tissue)
		#Fibrosarcoma (HT1080)
		if tissue == "cancer_fibrosarcoma_ht1080":
			download_ref("BSS00709", tissue)
		#Hepatocellular carcinoma (HuH-7.5)
		if tissue == "cancer_hepatocellular_carcinoma_huh-7.5":
			download_ref("BSS00719", tissue)
		#Endometrial ademocarcinoma (Ishikawa treated with 0.02pct dimethyl sulfoxide for 1 hour)
		if tissue == "cancer_endometrial_adenocarcinoma_ishikawa_treated_dmso_1hr":
			download_ref("BSS00745", tissue)
		#Endometrial ademocarcinoma (Ishikawa treated with 10 nM 17b-estradiol for 30 minutes)
		if tissue == "cancer_endometrial_adenocarcinoma_ishikawa_treated_17b-estradiol_30min":
			download_ref("BSS00748", tissue)
		#Endometrial ademocarcinoma (Ishikawa treated with 600 nM afimoxifene for 30 minutes)
		if tissue == "cancer_endometrial_adenocarcinoma_ishikawa_treated_afimoxifene_30min":
			download_ref("BSS00756", tissue)
		#Myelogenous leukemia (K562)
		if tissue == "cancer_myelogenous_leukemia_k562":
			download_ref("BSS00762", tissue)
		#Myelogenous leukemia (K562 G1 phase)
		if tissue == "cancer_myelogenous_leukemia_k562_g1_phase":
			download_ref("BSS01038", tissue)
		#Myelogenous leukemia (K562 G2 phase)
		if tissue == "cancer_myelogenous_leukemia_k562_g2_phase":
			download_ref("BSS01039", tissue)
		#Myelogenous leukemia (K562 treated with 0.05pct dimethyl sulfoxide for 72 hours)
		if tissue == "cancer_myelogenous_leukemia_k562_treated_dmso_72hr":
			download_ref("BSS01056", tissue)
		#Myelogenous leukemia (K562 treated with 1 mM vorinostat for 72 hours)
		if tissue == "cancer_myelogenous_leukemia_k562_treated_vorinostat_72hr":
			download_ref("BSS01057", tissue)
		#Myelogenous leukemia (K562 treated with 500 mM sodium butyrate for 72 hours)
		if tissue == "cancer_myelogenous_leukemia_k562_treated_sodium_butyrate_72hr":
			download_ref("BSS01059", tissue)
		#B cell lymphoma (Karpas-422)
		if tissue == "cancer_b_cell_lymphoma_karpas-422":
			download_ref("BSS01065", tissue)
		#Myelogenous leukemia (KBM-7)
		if tissue == "cancer_myelogenous_leukemia_kbm-7":
			download_ref("BSS01066", tissue)
		#Myeloma (KMS-11)
		if tissue == "cancer_myeloma_kms-11":
			download_ref("BSS01104", tissue)
		#Acute lymphoblastic leukemia (KOPT-K1)
		if tissue == "cancer_acute_lymphoblastic_leukemia_kopt-k1":
			download_ref("BSS01105", tissue)
		#Prostate adenocarcinoma (LNCaP clone FGC)
		if tissue == "cancer_prostate_adenocarcinoma_lncap_clone_fgc":
			download_ref("BSS01173", tissue)
		#Prostate adenocarcinoma (LNCaP clone FGC treated with 1 nM 17b-hydroxy-17-methylestra-4,9,11-trien-3-one for 12 hours)
		if tissue == "cancer_prostate_adenocarcinoma_lncap_clone_fgc_treated_17b_12hr":
			download_ref("BSS01174", tissue)
		#Acute lymphoblastic leukemia (Loucy)
		if tissue == "cancer_acute_lymphoblastic_leukemia_loucy":
			download_ref("BSS01178", tissue)
		#Colorectal adenocarcinoma (LoVo)
		if tissue == "cancer_colorectal_adenocarcinoma_lovo":
			download_ref("BSS01179", tissue)
		#Glioblastoma (M059J)
		if tissue == "cancer_glioblastoma_m059j":
			download_ref("BSS01208", tissue)
		#Mammary gland adenocarcinoma (MCF-7)
		if tissue == "cancer_mammary_gland_adenocarcinoma_mcf-7":
			download_ref("BSS01226", tissue)
		#Mammary gland adenocarcinoma (MCF-7 originated from MCF-7)
		if tissue == "cancer_mammary_gland_adenocarcinoma_mcf-7_originated_from_mcf-7":
			download_ref("BSS01235", tissue)
		#Mammary gland adenocarcinoma (MCF-7 treated with 10 mM lactate for 24 hours)
		if tissue == "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_lactate_24hr":
			download_ref("BSS01240", tissue)
		#Mammary gland adenocarcinoma (MCF-7 treated with 100 nM estradiol for 1 hour)
		if tissue == "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_estradiol_1hr":
			download_ref("BSS01243", tissue)
		#Mammary gland adenocarcinoma (MCF-7 treated with CTCF shRNA knockdown)
		if tissue == "cancer_mammary_gland_adenocarcinoma_mcf-7_treated_ctcf_shrna_knockown":
			download_ref("BSS01244", tissue)
		#Medulloblastoma (medulloblastoma)
		if tissue == "cancer_medulloblastoma":
			download_ref("BSS01251", tissue)
		#Osteosarcoma (MG63)
		if tissue == "cancer_osteosarcoma_mg63":
			download_ref("BSS01267", tissue)
		#Burkitt lymphoma (NAMALWA)
		if tissue == "cancer_burkitt_lymphoma_namalwa":
			download_ref("BSS01350", tissue)
		#Burkitt lymphoma (NAMALWA treated with Sendai virus for 2 hours)
		if tissue == "cancer_burkitt_lymphoma_namalwa_treated_sendai_virus_2hr":
			download_ref("BSS01351", tissue)
		#Acute promyelocytic leukemia (NB4)
		if tissue == "cancer_acute_promyelocytic_leukemia_nb4":
			download_ref("BSS01356", tissue)
		#Squamous cell carcinoma (NCI-H226)
		if tissue == "cancer_squamous_cell_carcinoma_nci-h226":
			download_ref("BSS01359", tissue)
		#Large cell lung cancer (NCI-H460)
		if tissue == "cancer_large_cell_lung_nci-h460":
			download_ref("BSS01360", tissue)
		#Myeloma (NCI-H929)
		if tissue == "cancer_myeloma_nci-h929":
			download_ref("BSS01365", tissue)
		#Testicular embryonal carcinoma (NT2 D1)
		if tissue == "cancer_testicular_embryonal_carcinoma_nt2_d1":
			download_ref("BSS01386", tissue)
		#B cell lymphoma (OCI-LY1)
		if tissue == "cancer_b_cell_lymphoma_oci-ly1":
			download_ref("BSS01389", tissue)
		#B cell lymphoma (OCI-LY3)
		if tissue == "cancer_b_cell_lymphoma_oci-ly3":
			download_ref("BSS01390", tissue)
		#B cell lymphoma (OCI-LY7)
		if tissue == "cancer_b_cell_lymphoma_oci-ly7":
			download_ref("BSS01391", tissue)
		#Pancreas duct epithelial carcinoma (Panc1)
		if tissue == "cancer_pancreas_duct_epithelial_carcinoma_panc1":
			download_ref("BSS01405", tissue)
		#Parathyroid adenoma (Parathyroid adenoma male adult (62 years))
		if tissue == "cancer_parathyroid_adenoma_62m":
			download_ref("BSS01411", tissue)
		#Prostate adenocarcinoma (PC-3)
		if tissue == "cancer_prostate_adenocarcinoma_pc-3":
			download_ref("BSS01414", tissue)
		#Lung adenocarcinoma (PC-9)
		if tissue == "cancer_lung_adenocarcinoma_pc-9":
			download_ref("BSS01415", tissue)
		#Renal cell adenocarcinoma (RCC 7860)
		if tissue == "cancer_renal_cell_adenocarcinoma_rcc_7860":
			download_ref("BSS01474", tissue)
		#Renal cell carcinoma (renal cell carcinoma)
		if tissue == "cancer_renal_cell_carcinoma":
			download_ref("BSS01481", tissue)
		#Colon carcinoma (RKO)
		if tissue == "cancer_colon_carcinoma_rko":
			download_ref("BSS01535", tissue)
		#Melanoma (RPMI-7951)
		if tissue == "cancer_melanoma_rpmi-7951":
			download_ref("BSS01536", tissue)
		#Plasma cell myeloma (RPMI8226)
		if tissue == "cancer_plasma_cell_myeloma_rpmi8226":
			download_ref("BSS01537", tissue)
		#Rhabdomyosarcoma (SJCRH30)
		if tissue == "cancer_rhabdomyosarcoma_sjcrh30":
			download_ref("BSS01549", tissue)
		#Osteosarcoma (SJSA1)
		if tissue == "cancer_osteosarcoma_sjsa1":
			download_ref("BSS01550", tissue)
		#Melanoma (SK-MEL-5)
		if tissue == "cancer_melanoma_sk-mel-5":
			download_ref("BSS01551", tissue)
		#Neuroblastoma (SK-N-DZ)
		if tissue == "cancer_neuroblastoma_sk-n-dz":
			download_ref("BSS01554", tissue)
		#Neuroblastoma (SK-N-DZ treated with dimethyl sulfoxide for 72 hours)
		if tissue == "cancer_neuroblastoma_sk-n-dz_treated_dmso_72hr":
			download_ref("BSS01558", tissue)
		#Neuroepithelioma (SK-N-MC)
		if tissue == "cancer_neuroepithelioma_sk-n-mc":
			download_ref("BSS01559", tissue)
		#Neuroblastoma (SK-N-SH)
		if tissue == "cancer_neuroblastoma_sk-n-sh":
			download_ref("BSS01562", tissue)
		#Neuroblastoma (SK-N-SH treated with 6 mM all-trans-retinoic acid for 48 hours)
		if tissue == "cancer_neuroblastoma_sk-n-sh_treated_retinoic_acid_48hr":
			download_ref("BSS01571", tissue)
		#B cell lymphoma (SU-DHL-6)
		if tissue == "cancer_b_cell_lymphoma_su-dhl-6":
			download_ref("BSS01664", tissue)
		#Colorectal adenocarcinoma (SW480)
		if tissue == "cancer_colorectal_adenocarcinoma_sw480":
			download_ref("BSS01682", tissue)
		#Mammary gland ductal carcinoma (T47D)
		if tissue == "cancer_mammary_gland_ductal_carcinoma_t47d":
			download_ref("BSS01699", tissue)
		#Mammary gland ductal carcinoma (T47D treated with 10 nM 17b-estradiol for 30 minutes)
		if tissue == "cancer_mammary_gland_ductal_carcinoma_t47d_treated_17b-estradiol_30min":
			download_ref("BSS01705", tissue)
		#Prostate epithelial carcinoma (VCaP)
		if tissue == "cancer_prostate_epithelial_carcinoma_vcap":
			download_ref("BSS01888", tissue)
		#Eye retinoblastoma (WERI-Rb-1)
		if tissue == "cancer_eye_retinoblastoma_weri-rb-1":
			download_ref("BSS01890", tissue)
			
		##Digestive
		#Colon mucosa (colonic mucosa female adult (56 years))
		if tissue == "digestive_colon_mucosa_56f":
			download_ref("BSS00227", tissue)
		#Colon mucosa (colonic mucosa female adult (73 years))
		if tissue == "digestive_colon_mucosa_73f":
			download_ref("BSS00228", tissue)
		#Duodenum mucosa (duodenal mucosa male adult (59 years))
		if tissue == "digestive_duodenum_mucosa_59m":
			download_ref("BSS00270", tissue)
		#Duodenum mucosa (duodenal mucosa male adult (76 years))
		if tissue == "digestive_duodenum_mucosa_76m":
			download_ref("BSS00271", tissue)
		#Esophagus (esophagus female adult (30 years))
		if tissue == "digestive_esophagus_30f":
			download_ref("BSS00316", tissue)
		#Esophagus (esophagus male adult (34 years))
		if tissue == "digestive_esophagus_34m":
			download_ref("BSS00318", tissue)
		#Esophagus muscularis mucosa (esophagus muscularis mucosa female adult (53 years))
		if tissue == "digestive_esophagus_muscularis_mucosa_53f":
			download_ref("BSS00321", tissue)
		#Esophagus muscularis mucosa (esophagus muscularis mucosa male adult (37 years))
		if tissue == "digestive_esophagus_muscularis_mucosa_37m":
			download_ref("BSS00322", tissue)
		#Esophagus squamous epithelium (esophagus squamous epithelium female adult (53 years))
		if tissue == "digestive_esophagus_squamous_epithelium_53f":
			download_ref("BSS00325", tissue)
		#Esophagus squamous epithelium (esophagus squamous epithelium male adult (37 years))
		if tissue == "digestive_esophagus_squamous_epithelium_37m":
			download_ref("BSS00326", tissue)
		#Gastroesophageal sphincter (gastroesophageal sphincter female adult (53 years))
		if tissue == "digestive_gastroesophageal_sphincter_53f":
			download_ref("BSS00381", tissue)
		#Large intestine (large intestine female embryo (103 days))
		if tissue == "digestive_large_intestine_embryo_103_days_f":
			download_ref("BSS01109", tissue)
		#Large intestine (large intestine female embryo (105 days))
		if tissue == "digestive_large_intestine_embryo_105_days_f":
			download_ref("BSS01110", tissue)
		#Large intestine (large intestine female embryo (107 days))
		if tissue == "digestive_large_intestine_embryo_107_days_f":
			download_ref("BSS01111", tissue)
		#Large intestine (large intestine female embryo (108 days))
		if tissue == "digestive_large_intestine_embryo_108_days_f":
			download_ref("BSS01112", tissue)
		#Large intestine (large intestine female embryo (110 days))
		if tissue == "digestive_large_intestine_embryo_110_days_f":
			download_ref("BSS01113", tissue)
		#Large intestine (large intestine female embryo (120 days))
		if tissue == "digestive_large_intestine_embryo_120_days_f":
			download_ref("BSS01114", tissue)
		#Large intestine (large intestine female embryo (91 day))
		if tissue == "digestive_large_intestine_embryo_91_days_f":
			download_ref("BSS01116", tissue)
		#Large intestine (large intestine female embryo (98 days))
		if tissue == "digestive_large_intestine_embryo_98_days_f":
			download_ref("BSS01117", tissue)
		#Large intestine (large intestine male embryo (105 days))
		if tissue == "digestive_large_intestine_embryo_105_days_m":
			download_ref("BSS01118", tissue)
		#Large intestine (large intestine male embryo (108 days))
		if tissue == "digestive_large_intestine_embryo_108_days_m":
			download_ref("BSS01119", tissue)
		#Large intestine (large intestine male embryo (113 days))
		if tissue == "digestive_large_intestine_embryo_113_days_m":
			download_ref("BSS01120", tissue)
		#Large intestine (large intestine male embryo (115 days))
		if tissue == "digestive_large_intestine_embryo_115_days_m":
			download_ref("BSS01121", tissue)
		#Large intestine (large intestine male embryo (91 day))
		if tissue == "digestive_large_intestine_embryo_91_days_m":
			download_ref("BSS01122", tissue)
		#Rectum mucosa (mucosa of rectum female adult (50 years))
		if tissue == "digestive_rectum_mucosa_50f":
			download_ref("BSS01282", tissue)
		#Rectum mucosa (mucosa of rectum female adult (61 year))
		if tissue == "digestive_rectum_mucosa_61f":
			download_ref("BSS01283", tissue)
		#Rectum mucosa (mucosa of stomach male adult (59 years))
		if tissue == "digestive_rectum_mucosa_59m":
			download_ref("BSS01284", tissue)
		#Stomach mucosa (mucosa of stomach male adult (59 years))
		if tissue == "digestive_stomach_mucosa_59m":
			download_ref("BSS01284", tissue)
		#Peyer's patch (Peyer's patch female adult (53 years))
		if tissue == "digestive_peyers_patch_53f":
			download_ref("BSS01426", tissue)
		#Peyer's patch (Peyer's patch male adult (37 years))
		if tissue == "digestive_peyers_patch_37m":
			download_ref("BSS01427", tissue)
		#Peyer's patch (Peyer's patch male adult (54 years))
		if tissue == "digestive_peyers_patch_54m":
			download_ref("BSS01428", tissue)
		#Sigmoid colon (sigmoid colon female adult (51 year))
		if tissue == "digestive_sigmoid_colon_51f":
			download_ref("BSS01542", tissue)
		#Sigmoid colon (sigmoid colon female adult (53 years))
		if tissue == "digestive_sigmoid_colon_53f":
			download_ref("BSS01543", tissue)
		#Sigmoid colon (sigmoid colon male adult (34 years))
		if tissue == "digestive_sigmoid_colon_34m":
			download_ref("BSS01545", tissue)
		#Sigmoid colon (sigmoid colon male adult (54 years))
		if tissue == "digestive_sigmoid_colon_54m":
			download_ref("BSS01547", tissue)
		#Sigmoid colon (sigmoid colon male child (3 years))
		if tissue == "digestive_sigmoid_colon_3m":
			download_ref("BSS01548", tissue)
		#Small intestine (small intestine female adult (30 years))
		if tissue == "digestive_small_intestine_30f":
			download_ref("BSS01588", tissue)
		#Small intestine (small intestine female embryo (105 days))
		if tissue == "digestive_small_intestine_105_days_f":
			download_ref("BSS01590", tissue)
		#Small intestine (small intestine female embryo (107 days))
		if tissue == "digestive_small_intestine_107_days_f":
			download_ref("BSS01591", tissue)
		#Small intestine (small intestine female embryo (108 days))
		if tissue == "digestive_small_intestine_108_days_f":
			download_ref("BSS01592", tissue)
		#Small intestine (small intestine female embryo (110 days))
		if tissue == "digestive_small_intestine_110_days_f":
			download_ref("BSS01593", tissue)
		#Small intestine (small intestine female embryo (120 days))
		if tissue == "digestive_small_intestine_120_days_f":
			download_ref("BSS01594", tissue)
		#Small intestine (small intestine female embryo (91 day))
		if tissue == "digestive_small_intestine_91_days_f":
			download_ref("BSS01595", tissue)
		#Small intestine (small intestine female embryo (98 days))
		if tissue == "digestive_small_intestine_98_days_f":
			download_ref("BSS01596", tissue)
		#Small intestine (small intestine male adult (34 years))
		if tissue == "digestive_small_intestine_34m":
			download_ref("BSS01597", tissue)
		#Small intestine (small intestine male child (3 years))
		if tissue == "digestive_small_intestine_3m":
			download_ref("BSS01599", tissue)
		#Small intestine (small intestine male embryo (105 days))
		if tissue == "digestive_small_intestine_105_days_m":
			download_ref("BSS01600", tissue)
		#Small intestine (small intestine male embryo (108 days))
		if tissue == "digestive_small_intestine_108_days_m":
			download_ref("BSS01601", tissue)
		#Small intestine (small intestine male embryo (115 days))
		if tissue == "digestive_small_intestine_115_days_m":
			download_ref("BSS01602", tissue)
		#Small intestine (small intestine male embryo (87 days))
		if tissue == "digestive_small_intestine_87_days_m":
			download_ref("BSS01603", tissue)
		#Small intestine (small intestine male embryo (91 days))
		if tissue == "digestive_small_intestine_91_days_m":
			download_ref("BSS01604", tissue)
		#Stomach (stomach embryo (101 day))
		if tissue == "digestive_stomach_101_days":
			download_ref("BSS01636", tissue)
		#Stomach (stomach female adult (30 years))
		if tissue == "digestive_stomach_30f":
			download_ref("BSS01637", tissue)
		#Stomach (stomach female adult (51 year))
		if tissue == "digestive_stomach_51f":
			download_ref("BSS01638", tissue)
		#Stomach (stomach female adult (53 years))
		if tissue == "digestive_stomach_53f":
			download_ref("BSS01639", tissue)
		#Stomach (stomach female embryo)
		if tissue == "digestive_stomach_f":
			download_ref("BSS01641", tissue)
		#Stomach (stomach female embryo (105 days))
		if tissue == "digestive_stomach_105_days_f":
			download_ref("BSS01642", tissue)
		#Stomach (stomach female embryo (107 days))
		if tissue == "digestive_stomach_107_days_f":
			download_ref("BSS01643", tissue)
		#Stomach (stomach female embryo (108 days))
		if tissue == "digestive_stomach_108_days_f":
			download_ref("BSS01644", tissue)
		#Stomach (stomach female embryo (121 day))
		if tissue == "digestive_stomach_121_days_f":
			download_ref("BSS01646", tissue)
		#Stomach (stomach female embryo (147 days))
		if tissue == "digestive_stomach_147_days_f":
			download_ref("BSS01647", tissue)
		#Stomach (stomach female embryo (96 days))
		if tissue == "digestive_stomach_96_days_f":
			download_ref("BSS01649", tissue)
		#Stomach (stomach female embryo (98 days))
		if tissue == "digestive_stomach_98_days_f":
			download_ref("BSS01650", tissue)
		#Stomach (stomach male adult (34 years))
		if tissue == "digestive_stomach_34m":
			download_ref("BSS01651", tissue)
		#Stomach (stomach male adult (54 years))
		if tissue == "digestive_stomach_54m":
			download_ref("BSS01653", tissue)
		#Stomach (stomach male child (3 years))
		if tissue == "digestive_stomach_3m":
			download_ref("BSS01654", tissue)
		#Stomach (stomach male embryo (108 days))
		if tissue == "digestive_stomach_108_days_m":
			download_ref("BSS01655", tissue)
		#Stomach (stomach male embryo (127 days))
		if tissue == "digestive_stomach_127_days_m":
			download_ref("BSS01656", tissue)
		#Stomach (stomach male embryo (58 days) and male embryo (76 days))
		if tissue == "digestive_stomach_58_76_days_m":
			download_ref("BSS01657", tissue)
		#Stomach (stomach male embryo (91 day))
		if tissue == "digestive_stomach_91_days_m":
			download_ref("BSS01658", tissue)
		#Transverse colon (transverse colon female adult (51 year))
		if tissue == "digestive_transverse_colon_51f":
			download_ref("BSS01848", tissue)
		#Transverse colon (transverse colon female adult (53 years))
		if tissue == "digestive_transverse_colon_53f":
			download_ref("BSS01849", tissue)
		#Transverse colon (transverse colon male adult (37 years))
		if tissue == "digestive_transverse_colon_37m":
			download_ref("BSS01850", tissue)
		#Transverse colon (transverse colon male adult (54 years))
		if tissue == "digestive_transverse_colon_54m":
			download_ref("BSS01851", tissue)
			
		##Endocrine
		#Adrenal gland (adrenal gland embryo (96 days))
		if tissue == "endocrine_adrenal_gland_96_days":
			download_ref("BSS00045", tissue)
		#Adrenal gland (adrenal gland female adult (30 years))
		if tissue == "endocrine_adrenal_gland_30f":
			download_ref("BSS00046", tissue)
		#Adrenal gland (adrenal gland female adult (51 year))
		if tissue == "endocrine_adrenal_gland_51f":
			download_ref("BSS00047", tissue)
		#Adrenal gland (adrenal gland female adult (53 years))
		if tissue == "endocrine_adrenal_gland_53f":
			download_ref("BSS00048", tissue)
		#Adrenal gland (adrenal gland female embryo (108 days))
		if tissue == "endocrine_adrenal_gland_108_days_f":
			download_ref("BSS00050", tissue)
		#Adrenal gland (adrenal gland female embryo (113 days))
		if tissue == "endocrine_adrenal_gland_113_days_f":
			download_ref("BSS00051", tissue)
		#Adrenal gland (adrenal gland female embryo (85 days))
		if tissue == "endocrine_adrenal_gland_85_days_f":
			download_ref("BSS00052", tissue)
		#Adrenal gland (adrenal gland male adult (34 years))
		if tissue == "endocrine_adrenal_gland_34m":
			download_ref("BSS00054", tissue)
		#Adrenal gland (adrenal gland male adult (37 years))
		if tissue == "endocrine_adrenal_gland_37m":
			download_ref("BSS00055", tissue)
		#Adrenal gland (adrenal gland male adult (54 years))
		if tissue == "endocrine_adrenal_gland_54m":
			download_ref("BSS00056", tissue)
		#Adrenal gland (adrenal gland male embryo (101 day))
		if tissue == "endocrine_adrenal_gland_101_days_m":
			download_ref("BSS00057", tissue)
		#Adrenal gland (adrenal gland male embryo (108 days))
		if tissue == "endocrine_adrenal_gland_108_days_m":
			download_ref("BSS00058", tissue)
		#Adrenal gland (adrenal gland male embryo (85 days))
		if tissue == "endocrine_adrenal_gland_85_days_m":
			download_ref("BSS00059", tissue)
		#Adrenal gland (adrenal gland male embryo (97 days))
		if tissue == "endocrine_adrenal_gland_97_days_m":
			download_ref("BSS00060", tissue)
		#Pancreas (endocrine pancreas adult (59 years))
		if tissue == "endocrine_pancreas_59":
			download_ref("BSS00281", tissue)
		#Pancreas (endocrine pancreas male)
		if tissue == "endocrine_pancreas_m":
			download_ref("BSS00282", tissue)
		#Pancreas (endocrine pancreas male adult (45 years))
		if tissue == "endocrine_pancreas_45m":
			download_ref("BSS00283", tissue)
		#Pancreas (endocrine pancreas male adult (46 years))
		if tissue == "endocrine_pancreas_46m":
			download_ref("BSS00284", tissue)
		#Ovary (ovary female adult (30 years))
		if tissue == "endocrine_ovary_30f":
			download_ref("BSS01399", tissue)
		#Ovary (ovary female adult (51 year))
		if tissue == "endocrine_ovary_51f":
			download_ref("BSS01401", tissue)
		#Ovary (ovary female adult (53 years))
		if tissue == "endocrine_ovary_53f":
			download_ref("BSS01402", tissue)
		#Ovary (ovary female embryo)
		if tissue == "endocrine_ovary_embryo_f":
			download_ref("BSS01403", tissue)
		#Testis (testis male adult (37 years))
		if tissue == "endocrine_testis_37m":
			download_ref("BSS01715", tissue)
		#Testis (testis male adult (54 years))
		if tissue == "endocrine_testis_54m":
			download_ref("BSS01718", tissue)
		#Testis (testis male embryo)
		if tissue == "endocrine_testis_embryo_m":
			download_ref("BSS01719", tissue)
		#Thyroid gland (thyroid gland female adult (51 year))
		if tissue == "endocrine_thyroid_gland_51f":
			download_ref("BSS01831", tissue)
		#Thyroid gland (thyroid gland female adult (53 years))
		if tissue == "endocrine_thyroid_gland_53f":
			download_ref("BSS01832", tissue)
		#Thyroid gland (thyroid gland male adult (37 years))
		if tissue == "endocrine_thyroid_gland_37m":
			download_ref("BSS01834", tissue)
		#Thyroid gland (thyroid gland male adult (54 years))
		if tissue == "endocrine_thyroid_gland_54m":
			download_ref("BSS01835", tissue)
			
		##Endothelial
		#Brain microvascular endothelial cell (brain microvascular endothelial cell)
		if tissue == "endothelial_brain_microvascular":
			download_ref("BSS00143", tissue)
		#Dermis blood vessel endothelial cell (dermis blood vessel endothelial cell female adult)
		if tissue == "endothelial_dermis_blood_vessel_adult_f":
			download_ref("BSS00258", tissue)
		#Dermis blood vessel endothelial cell (dermis blood vessel endothelial cell male newborn)
		if tissue == "endothelial_dermis_blood_vessel_newborn_m":
			download_ref("BSS00260", tissue)
		#Dermis blood vessel endothelial cell (dermis microvascular lymphatic vessel endothelial cell female)
		if tissue == "endothelial_dermis_microvascular_lymphatoc_vessel_f":
			download_ref("BSS00262", tissue)
		#Dermis blood vessel endothelial cell (dermis microvascular lymphatic vessel endothelial cell male)
		if tissue == "endothelial_dermis_microvascular_lymphatoc_vessel_m":
			download_ref("BSS00264", tissue)
		#Umbilical vein endothelial cell (endothelial cell of umbilical vein male newborn)
		if tissue == "endothelial_umbilical_vein_newborn_m":
			download_ref("BSS00296", tissue)
		#Umbilical vein endothelial cell (endothelial cell of umbilical vein newborn)
		if tissue == "endothelial_umbilical_vein_newborn":
			download_ref("BSS00298", tissue)
		#Glomerulus endothelial cell (glomerular endothelial cell)
		if tissue == "endothelial_glomerulus":
			download_ref("BSS00387", tissue)
		#Kidney capillary endothelial cell (kidney capillary endothelial cell female embryo (113 days))
		if tissue == "endothelial_kidney_capillary_113_days_f":
			download_ref("BSS01077", tissue)
		#Lung microvascular endothelial cell (lung microvascular endothelial cell female)
		if tissue == "endothelial_lung_microvascular_f":
			download_ref("BSS01206", tissue)
		#Pulmonary artery endothelial cell (pulmonary artery endothelial cell female)
		if tissue == "endothelial_pulmonary_artery_f":
			download_ref("BSS01465", tissue)	
			
		##Epithelial
		#Amnion epithelial cell (amniotic epithelial cell)
		if tissue == "epithelial_amnion":
			download_ref("BSS00075", tissue)
		#Bronchial epithelial cell (bronchial epithelial cell)
		if tissue == "epithelial_bronchial":
			download_ref("BSS00150", tissue)
		#Bronchial epithelial cell (bronchial epithelial cell female treated with retinoic acid)
		if tissue == "epithelial_bronchial_f_treated_retinoic_acid":
			download_ref("BSS00153", tissue)
		#Choroid plexus epithelial cell (choroid plexus epithelial cell)
		if tissue == "epithelial_choroid_plexus":
			download_ref("BSS00218", tissue)
		#Colon epithelial cell (colon epithelial cell line)
		if tissue == "epithelial_colon":
			download_ref("BSS00223", tissue)
		#Esophagus epithelial cell (epithelial cell of esophagus)
		if tissue == "epithelial_esophagus":
			download_ref("BSS00307", tissue)
		#Prostate epithelial cell (epithelial cell of prostate - primary cell - ENCODE)
		if tissue == "epithelial_prostate":
			download_ref("BSS00308", tissue)
		#Prostate epithelial cell (epithelial cell of prostate male)
		if tissue == "epithelial_prostate_m":
			download_ref("BSS00309", tissue)
		#Proximal tubule epithelial cell (epithelial cell of proximal tubule)
		if tissue == "epithelial_proximal_tubule":
			download_ref("BSS00310", tissue)
		#Foreskin keratinocyte (foreskin keratinocyte male newborn)
		if tissue == "epithelial_foreskin_keratinocyte_newborn_m":
			download_ref("BSS00354", tissue)
		#Foreskin keratinocyte (foreskin keratinocyte male newborn (2-4 days))
		if tissue == "epithelial_foreskin_keratinocyte_newborn_2-4_days_m":
			download_ref("BSS00355", tissue)
		#Foreskin keratinocyte (foreskin keratinocyte male newborn (2-4 days) treated with 1.2 mM calcium for 2.5 days)
		if tissue == "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_2.5_days":
			download_ref("BSS00361", tissue)
		#Foreskin keratinocyte (foreskin keratinocyte male newborn (2-4 days) treated with 1.2 mM calcium for 5.5 days)
		if tissue == "epithelial_foreskin_keratinocyte_newborn_2-4_days_m_treated_calcium_5.5_days":
			download_ref("BSS00367", tissue)
		#Glomerulus visceral epithelial cell (glomerular visceral epithelial cell child (3 years))
		if tissue == "epithelial_glomerulus_visceral_3":
			download_ref("BSS00389", tissue)
		#Proximal tubule epithelial cell (HK-2)
		if tissue == "epithelial_proximal_tubule_hk-2":
			download_ref("BSS00701", tissue)
		#Pancreatic duct epithelial cell (HPDE6-E6E7)
		if tissue == "epithelial_pancreatic_duct_dpde6-e6e7":
			download_ref("BSS00703", tissue)
		#Bone marrow epithelial cell (HS-27A)
		if tissue == "epithelial_bone_marrow_hs-27a":
			download_ref("BSS00704", tissue)
		#Iris pigment epithelial cell (iris pigment epithelial cell)
		if tissue == "epithelial_iris_pigment":
			download_ref("BSS00743", tissue)
		#Keratinocyte (keratinocyte female)
		if tissue == "epithelial_keratinocyte_f":
			download_ref("BSS01068", tissue)
		#Keratinocyte (keratinocyte male)
		if tissue == "epithelial_keratinocyte_m":
			download_ref("BSS01071", tissue)
		#Kidney epithelial cell (kidney epithelial cell)
		if tissue == "epithelial_kidney":
			download_ref("BSS01080", tissue)
		#Glomerulus epithelial cell (kidney glomerular epithelial cell male adult (43 years) and male adult (62 years))
		if tissue == "epithelial_glomerulus_43_62_m":
			download_ref("BSS01092", tissue)
		#Tubule cell (kidney tubule cell female adult (80 years) and male adult (62 years))
		if tissue == "epithelial_tubule_80f_62m":
			download_ref("BSS01102", tissue)
		#Tubule cell (kidney tubule cell female adult (80 years) treated with 5 mM cisplatin)
		if tissue == "epithelial_tubule_80f_treated_cisplatin":
			download_ref("BSS01103", tissue)
		#Skin leg (lower leg skin female adult (53 years))
		if tissue == "epithelial_skin_leg_53f":
			download_ref("BSS01181", tissue)
		#Skin leg (lower leg skin male adult (37 years))
		if tissue == "epithelial_skin_leg_37m":
			download_ref("BSS01182", tissue)
		#Mammary luminal epithelial cell (luminal epithelial cell of mammary gland female adult (33 years))
		if tissue == "epithelial_mammary_luminal_33f":
			download_ref("BSS01185", tissue)
		#Mammary epithelial cell (mammary epithelial cell female)
		if tissue == "epithelial_mammary_f":
			download_ref("BSS01209", tissue)
		#Mammary epithelial cell (mammary epithelial cell female adult (18 years))
		if tissue == "epithelial_mammary_18f":
			download_ref("BSS01211", tissue)
		#Mammary epithelial cell (mammary epithelial cell female adult (50 years))
		if tissue == "epithelial_mammary_50f":
			download_ref("BSS01213", tissue)
		#Breast epithelial cell (MCF 10A)
		if tissue == "epithelial_breast_mcf_10a":
			download_ref("BSS01217", tissue)
		#Breast epithelial cell (MCF 10A treated with 1 mM tamoxifen for 24 hours)
		if tissue == "epithelial_breast_mcf_10a_treated_tamoxifen_24hr":
			download_ref("BSS01224", tissue)
		#Breast epithelial cell (MCF 10A treated with 1 mM tamoxifen for 6 hours)
		if tissue == "epithelial_breast_mcf_10a_treated_tamoxifen_6hr":
			download_ref("BSS01225", tissue)
		#Mammary myoepithelial cell (myoepithelial cell of mammary gland female adult (33 years))
		if tissue == "epithelial_mammary_myoepithelial_33f":
			download_ref("BSS01340", tissue)
		#Mammary myoepithelial cell (myoepithelial cell of mammary gland female adult (36 years))
		if tissue == "epithelial_mammary_myoepithelial_36f":
			download_ref("BSS01341", tissue)
		#Non-pigmented ciliary epithelial cell (non-pigmented ciliary epithelial cell)
		if tissue == "epithelial_non-pigmented_ciliary":
			download_ref("BSS01385", tissue)
		#Renal cortical epithelial cell (renal cortical epithelial cell)
		if tissue == "epithelial_renal_cortical":
			download_ref("BSS01491", tissue)
		#Retinal epithelial cell (retinal pigment epithelial cell)
		if tissue == "epithelial_retinal":
			download_ref("BSS01505", tissue)
		#Prostate epithelial cell (RWPE1)
		if tissue == "epithelial_prostate_rwpe1":
			download_ref("BSS01538", tissue)
		#Prostate epithelial cell (RWPE2)
		if tissue == "epithelial_prostate_rwpe2":
			download_ref("BSS01539", tissue)
		#Skin of body (skin of body female embryo (82 days))
		if tissue == "epithelial_skin_of_body_82_days_f":
			download_ref("BSS01587", tissue)
			
		##ESC-deriv
		#Bipolar neuron deriv (bipolar neuron originated from GM23338 treated with 0.5 mg mL doxycycline hyclate for 4 days)
		if tissue == "esc_deriv_bipolar_neuron_gm23338_origin_treated_doxycycline_4_days":
			download_ref("BSS00112", tissue)
		#Cardiac mesoderm deriv (cardiac mesoderm originated from H7-hESC)
		if tissue == "esc_deriv_cardiac_mesoderm_h7-hesc_origin":
			download_ref("BSS00169", tissue)
		#Cardiac muscle deriv (cardiac muscle cell originated from RUES2)
		if tissue == "esc_deriv_cardiac_muscle_rues2_origin":
			download_ref("BSS00171", tissue)
		#Neural progenitor deriv (ecto neural progenitor cell originated from H9)
		if tissue == "esc_deriv_neural_progenitor_h9_origin_1":
			download_ref("BSS00272", tissue)
		#Ectodermal deriv (ectodermal cell originated from embryonic stem cell)
		if tissue == "esc_deriv_ectodermal":
			download_ref("BSS00273", tissue)
		#Endodermal cell (endodermal cell)
		if tissue == "esc_deriv_endodermal":
			download_ref("BSS00285", tissue)
		#Endodermal cell (endodermal cell originated from HUES64)
		if tissue == "esc_deriv_endodermal_hues64_origin":
			download_ref("BSS00287", tissue)
		#Hepatocyte deriv (hepatocyte originated from H9)
		if tissue == "esc_deriv_hepatocyte_h9_origin":
			download_ref("BSS00556", tissue)
		#Mesenchymal stem deriv (mesenchymal stem cell originated from H1-hESC)
		if tissue == "esc_deriv_mesenchymal_stem_h1-hesc_origin":
			download_ref("BSS01261", tissue)
		#Mesendoderm deriv (mesendoderm originated from H1-hESC)
		if tissue == "esc_deriv_mesendoderm_h1-hesc_origin":
			download_ref("BSS01263", tissue)
		#Mesodermal deriv (mesodermal cell originated from HUES64)
		if tissue == "esc_deriv_mesodermal_hues64_origin":
			download_ref("BSS01264", tissue)
		#Neural deriv (neural cell originated from H1-hESC)
		if tissue == "esc_deriv_neural_crest_h1-hesc_origin":
			download_ref("BSS01366", tissue)
		#Neural progenitor deriv (neural progenitor cell originated from H9)
		if tissue == "esc_deriv_neural_progenitor_h9_origin_2":
			download_ref("BSS01370", tissue)
		#Neural progenitor deriv (neural stem progenitor cell originated from H1-hESC)
		if tissue == "esc_deriv_neural_progenitor_h1-hesc_origin":
			download_ref("BSS01371", tissue)
		#Neural progenitor deriv (neural stem progenitor cell originated from H9)
		if tissue == "esc_deriv_neural_progenitor_h9_origin_3":
			download_ref("BSS01372", tissue)
		#Neuron deriv (neuron originated from H9)
		if tissue == "esc_deriv_neuron_h9_origin":
			download_ref("BSS01375", tissue)
		#Smooth muscle deriv (smooth muscle cell originated from H9)
		if tissue == "esc_deriv_smooth_muscle_h9_origin":
			download_ref("BSS01612", tissue)
		#Trophoblast deriv (trophoblast cell originated from H1-hESC)
		if tissue == "esc_deriv_trophoblast_h1-hesc_origin":
			download_ref("BSS01857", tissue)
			
		##ESC
		#ESC (ELF-1)
		if tissue == "esc_elf-1":
			download_ref("BSS00277", tissue)
		#ESC (ES-I3)
		if tissue == "esc_es-i3":
			download_ref("BSS00315", tissue)
		#ESC (H1-hESC)
		if tissue == "esc_h1-hesc":
			download_ref("BSS00478", tissue)
		#ESC (H7-hESC)
		if tissue == "esc_h7-hesc":
			download_ref("BSS00483", tissue)
		#ESC (H9)
		if tissue == "esc_h9":
			download_ref("BSS00484", tissue)
		#ESC (HUES48)
		if tissue == "esc_hues48":
			download_ref("BSS00715", tissue)
		#ESC (HUES6)
		if tissue == "esc_hues6":
			download_ref("BSS00716", tissue)
		#ESC (HUES64)
		if tissue == "esc_hues64":
			download_ref("BSS00717", tissue)
		#ESC (UCSF-4)
		if tissue == "esc_ucsf-4":
			download_ref("BSS01866", tissue)
			
		##Eye
		#Eye (eye embryo (56 days) and male embryo (76 days))
		if tissue == "eye_56_days_76_days_m":
			download_ref("BSS00328", tissue)
		#Eye (eye female embryo (76 days))
		if tissue == "eye_76_days_f":
			download_ref("BSS00329", tissue)
		#Eye (retina embryo (125 days) and male embryo (103 days))
		if tissue == "eye_125_days_103_days_m":
			download_ref("BSS01502", tissue)
		#Eye (retina embryo (74 days) and embryo (85 days))
		if tissue == "eye_74_days_85_days":
			download_ref("BSS01503", tissue)
		#Eye (retina female embryo (89 days))
		if tissue == "eye_89_days_f":
			download_ref("BSS01504", tissue)
		
		##Heart
		#Aorta (aorta female adult (30 years))
		if tissue == "heart_aorta_30f":
			download_ref("BSS00079", tissue)
		#Aorta (aorta male adult (34 years))
		if tissue == "heart_aorta_34m":
			download_ref("BSS00080", tissue)
		#Ascending aorta (ascending aorta female adult (51 year))
		if tissue == "heart_ascending_aorta_51f":
			download_ref("BSS00087", tissue)
		#Ascending aorta (ascending aorta female adult (53 years))
		if tissue == "heart_ascending_aorta_53f":
			download_ref("BSS00088", tissue)
		#Coronary artery (coronary artery female adult (51 year))
		if tissue == "heart_coronary_artery_51f":
			download_ref("BSS00242", tissue)
		#Coronary artery (coronary artery female adult (53 years))
		if tissue == "heart_coronary_artery_53f":
			download_ref("BSS00243", tissue)
		#Heart (heart embryo (101 day))
		if tissue == "heart_101_days":
			download_ref("BSS00493", tissue)
		#Heart (heart embryo (59 days) and female embryo (76 days))
		if tissue == "heart_59_days_76_days_f":
			download_ref("BSS00494", tissue)
		#Heart (heart embryo (80 days))
		if tissue == "heart_80_days":
			download_ref("BSS00495", tissue)
		#Heart (heart embryo (96 days))
		if tissue == "heart_96_days":
			download_ref("BSS00496", tissue)
		#Heart (heart female embryo (103 days))
		if tissue == "heart_103_days_f":
			download_ref("BSS00498", tissue)
		#Heart (heart female embryo (105 days))
		if tissue == "heart_105_days_f":
			download_ref("BSS00499", tissue)
		#Heart (heart female embryo (110 days))
		if tissue == "heart_110_days_f":
			download_ref("BSS00500", tissue)
		#Heart (heart female embryo (116 days) and female embryo (98 days))
		if tissue == "heart_116_days_98_days_f":
			download_ref("BSS00501", tissue)
		#Heart (heart female embryo (117 days))
		if tissue == "heart_116_days_117_days_f":
			download_ref("BSS00502", tissue)
		#Heart (heart female embryo (147 days))
		if tissue == "heart_147_days_f":
			download_ref("BSS00503", tissue)
		#Heart (heart female embryo (91 day))
		if tissue == "heart_91_days_f":
			download_ref("BSS00505", tissue)
		#Heart left ventricle (heart left ventricle female adult (53 years))
		if tissue == "heart_left_ventricle_53f":
			download_ref("BSS00507", tissue)
		#Heart left ventricle (heart left ventricle female embryo (101 day) and female embryo (103 days))
		if tissue == "heart_left_ventricle_101_days_103_days_f":
			download_ref("BSS00508", tissue)
		#Heart left ventricle (heart left ventricle female embryo (136 days))
		if tissue == "heart_left_ventricle_136_days_f":
			download_ref("BSS00509", tissue)
		#Heart left ventricle (heart left ventricle male adult (34 years))
		if tissue == "heart_left_ventricle_34m":
			download_ref("BSS00512", tissue)
		#Heart left ventricle (heart left ventricle male child (3 years))
		if tissue == "heart_left_ventricle_3m":
			download_ref("BSS00513", tissue)
		#Heart (heart male adult (27 years) and male adult (35 years))
		if tissue == "heart_27m_35m":
			download_ref("BSS00514", tissue)
		#Heart (heart male child (3 years))
		if tissue == "heart_3m":
			download_ref("BSS00516", tissue)	
		#Heart (heart male embryo (105 days))
		if tissue == "heart_105_days_m":
			download_ref("BSS00517", tissue)
		#Heart (heart male embryo (110 days))
		if tissue == "heart_110_days_m":
			download_ref("BSS00518", tissue)
		#Heart (heart male embryo (120 days))
		if tissue == "heart_120_days_m":
			download_ref("BSS00519", tissue)
		#Heart (heart male embryo (72 days) and male embryo (76 days))
		if tissue == "heart_72_days_76_days_m":
			download_ref("BSS00520", tissue)
		#Heart (heart male embryo (91 day))
		if tissue == "heart_91_days_m":
			download_ref("BSS00521", tissue)
		#Heart (heart male embryo (96 days))
		if tissue == "heart_96_days_m":
			download_ref("BSS00522", tissue)
		#Heart right ventricle (heart right ventricle female embryo (101 day) and female embryo (103 days))
		if tissue == "heart_right_ventricle_101_days_103_days_f":
			download_ref("BSS00523", tissue)
		#Heart right ventricle (heart right ventricle male adult (34 years))
		if tissue == "heart_right_ventricle_34m":
			download_ref("BSS00524", tissue)
		#Heart right ventricle (heart right ventricle male child (3 years))
		if tissue == "heart_right_ventricle_3m":
			download_ref("BSS00525", tissue)
		#Heart left atrium (left cardiac atrium female embryo (101 day))
		if tissue == "heart_left_atrium_101_days_f":
			download_ref("BSS01127", tissue)
		#Heart right atrium (right atrium auricular region female adult (51 year))
		if tissue == "heart_right_atrium_51f":
			download_ref("BSS01506", tissue)
		#Heart right atrium (right atrium auricular region female adult (53 years))
		if tissue == "heart_right_atrium_53f":
			download_ref("BSS01507", tissue)
		#Heart right atrium (right cardiac atrium male adult (34 years))
		if tissue == "heart_right_atrium_34m":
			download_ref("BSS01508", tissue)
		#Thoracic aorta (thoracic aorta male adult (37 years))
		if tissue == "heart_thoracic_aorta_37m":
			download_ref("BSS01814", tissue)
		#Thoracic aorta (thoracic aorta male adult (54 years))
		if tissue == "heart_thoracic_aorta_54m":
			download_ref("BSS01815", tissue)
		#Tibial artery (tibial artery female adult (53 years))
		if tissue == "heart_tibial_artery_53f":
			download_ref("BSS01837", tissue)
		#Tibial artery (tibial artery male adult (37 years))
		if tissue == "heart_tibial_artery_37m":
			download_ref("BSS01838", tissue)
			
		##HSC & B-cell
		#B cell (B cell)
		if tissue == "hsc_and_b_cell_b_cell":
			download_ref("BSS00093", tissue)
		#B cell (B cell female adult (27 years))
		if tissue == "hsc_and_b_cell_b_cell_27f":
			download_ref("BSS00095", tissue)
		#B cell (B cell female adult (27 years) and female adult (43 years))
		if tissue == "hsc_and_b_cell_b_cell_27f_43f":
			download_ref("BSS00096", tissue)
		#B cell (B cell female adult (34 years))
		if tissue == "hsc_and_b_cell_b_cell_34f":
			download_ref("BSS00097", tissue)
		#B cell (B cell female adult (43 years))
		if tissue == "hsc_and_b_cell_b_cell_43f":
			download_ref("BSS00098", tissue)
		#B cell (B cell male adult (21 year))
		if tissue == "hsc_and_b_cell_b_cell_21m":
			download_ref("BSS00100", tissue)
		#B cell (B cell male adult (37 years))
		if tissue == "hsc_and_b_cell_b_cell_37m":
			download_ref("BSS00101", tissue)
		#CD14 monocyte (CD14-positive monocyte female)
		if tissue == "hsc_and_b_cell_cd14_monocyte_f":
			download_ref("BSS00178", tissue)
		#CD14 monocyte (CD14-positive monocyte female adult (34 years))
		if tissue == "hsc_and_b_cell_cd14_monocyte_34f":
			download_ref("BSS00179", tissue)
		#CD14 monocyte (CD14-positive monocyte male adult (21 year))
		if tissue == "hsc_and_b_cell_cd14_monocyte_21m":
			download_ref("BSS00180", tissue)
		#CD14 monocyte (CD14-positive monocyte male adult (37 years))
		if tissue == "hsc_and_b_cell_cd14_monocyte_37m":
			download_ref("BSS00181", tissue)
		#CD1c myeloid dendritic cell (CD1c-positive myeloid dendritic cell)
		if tissue == "hsc_and_b_cell_cd1c_myeloid_dendritic":
			download_ref("BSS00182", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive)
		if tissue == "hsc_and_b_cell_cmp_cd34":
			download_ref("BSS00229", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive female)
		if tissue == "hsc_and_b_cell_cmp_cd34_f":
			download_ref("BSS00230", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive female adult (27 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_27f":
			download_ref("BSS00231", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive female adult (33 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_33f":
			download_ref("BSS00232", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive female adult (50 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_50f":
			download_ref("BSS00233", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male)
		if tissue == "hsc_and_b_cell_cmp_cd34_m":
			download_ref("BSS00234", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult)
		if tissue == "hsc_and_b_cell_cmp_cd34_adult_m":
			download_ref("BSS00235", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (23 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_23m":
			download_ref("BSS00236", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (36 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_36m":
			download_ref("BSS00237", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (37 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_37m":
			download_ref("BSS00238", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (42 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_42m":
			download_ref("BSS00239", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (43 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_43m":
			download_ref("BSS00240", tissue)
		#CD34 CMP (common myeloid progenitor, CD34-positive male adult (49 years))
		if tissue == "hsc_and_b_cell_cmp_cd34_49m":
			download_ref("BSS00241", tissue)
		#Germinal center (germinal center)
		if tissue == "hsc_and_b_cell_germinal_center":
			download_ref("BSS00384", tissue)
		#MPP (hematopoietic multipotent progenitor cell)
		if tissue == "hsc_and_b_cell_mpp":
			download_ref("BSS00543", tissue)
		#MPP (hematopoietic multipotent progenitor cell male adult (25 years) treated with erythropoietin for 20 days, hydrocortisone succinate for 20 days, kit ligand for 20 days, interleukin-3 for 20 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_20_days":
			download_ref("BSS00544", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 11 day, kit ligand for 11 day, hydrocortisone succinate for 11 day, erythropoietin for 11 day)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_11_days":
			download_ref("BSS00545", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 13 days, kit ligand for 13 days, hydrocortisone succinate for 13 days, erythropoietin for 13 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_13_days":
			download_ref("BSS00546", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 15 days, kit ligand for 15 days, hydrocortisone succinate for 15 days, erythropoietin for 15 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_15_days":
			download_ref("BSS00547", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 17 days, kit ligand for 17 days, hydrocortisone succinate for 17 days, erythropoietin for 17 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_17_days":
			download_ref("BSS00548", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 18 days, kit ligand for 18 days, hydrocortisone succinate for 18 days, erythropoietin for 18 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_18_days":
			download_ref("BSS00549", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 4 days, kit ligand for 4 days, hydrocortisone succinate for 4 days, erythropoietin for 4 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_4_days":
			download_ref("BSS00550", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 6 days, kit ligand for 6 days, hydrocortisone succinate for 6 days, erythropoietin for 6 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_6_days":
			download_ref("BSS00551", tissue)
		#MPP (hematopoietic multipotent progenitor cell treated with interleukin-3 for 8 days, kit ligand for 8 days, hydrocortisone succinate for 8 days, erythropoietin for 8 days)
		if tissue == "hsc_and_b_cell_mpp_25m_treated_erythropoietin_hydrocortisone_succinate_ligand_il3_8_days":
			download_ref("BSS00552", tissue)
		#Lymphocyte (Jurkat clone E61)
		if tissue == "hsc_and_b_cell_lymphocyte_jurkat_clone_e61":
			download_ref("BSS00760", tissue)
		#Naive B cell (naive B cell)
		if tissue == "hsc_and_b_cell_lymphocyte_naive_b_cell":
			download_ref("BSS01345", tissue)
		#NK cell (natural killer cell female adult (34 years))
		if tissue == "hsc_and_b_cell_nk_34f":
			download_ref("BSS01353", tissue)
		#NK cell (natural killer cell male adult (21 year))
		if tissue == "hsc_and_b_cell_nk_21m":
			download_ref("BSS01354", tissue)
		#NK cell (natural killer cell male adult (37 years))
		if tissue == "hsc_and_b_cell_nk_37m":
			download_ref("BSS01355", tissue)
		#Neutrophil (neutrophil)
		if tissue == "hsc_and_b_cell_neutrophil":
			download_ref("BSS01380", tissue)
		#Neutrophil (neutrophil male)
		if tissue == "hsc_and_b_cell_neutrophil_m":
			download_ref("BSS01381", tissue)
			
		##iPSC
		#iPSC (CWRU1 male)
		if tissue == "ipsc_cwru1_m":
			download_ref("BSS00244", tissue)
		#iPSC (GM23338 male adult (53 years) originated from GM23248)
		if tissue == "ipsc_gm23338_53m_gm23248_origin":
			download_ref("BSS00477", tissue)
		#iPSC (iPS DF 19.11 male newborn)
		if tissue == "ipsc_ips_df_19.11_newborn_m":
			download_ref("BSS00731", tissue)
		#iPSC (iPS DF 19.7 male newborn)
		if tissue == "ipsc_ips_df_19.7_newborn_m":
			download_ref("BSS00732", tissue)
		#iPSC (iPS DF 4.7 male newborn)
		if tissue == "ipsc_ips_df_14.7_newborn_m":
			download_ref("BSS00733", tissue)
		#iPSC (iPS DF 6.9 male newborn)
		if tissue == "ipsc_ips_df_6.9_newborn_m":
			download_ref("BSS00734", tissue)
		#iPSC (iPS-11a male adult (36 years))
		if tissue == "ipsc_ips-11a_36m":
			download_ref("BSS00735", tissue)
		#iPSC (iPS-15b female adult (48 years))
		if tissue == "ipsc_ips-15b_48f":
			download_ref("BSS00736", tissue)
		#iPSC (iPS-18a female adult (48 years))
		if tissue == "ipsc_ips-18a_48f":
			download_ref("BSS00737", tissue)
		#iPSC (iPS-18c female adult (48 years))
		if tissue == "ipsc_ips-18c_48f":
			download_ref("BSS00738", tissue)
		#iPSC (iPS-20b male adult (55 years))
		if tissue == "ipsc_ips-20b_55m":
			download_ref("BSS00738", tissue)
		#iPSC (iPS-NIHi11 male adult (71 year) originated from AG20443)
		if tissue == "ipsc_ips-nihi11_71m_ag20443_origin":
			download_ref("BSS00741", tissue)
		#iPSC (iPS-NIHi7 female adult (85 years) originated from AG08395)
		if tissue == "ipsc_ips-nihi7_85f_ag08395_origin":
			download_ref("BSS00742", tissue)
		#iPSC (L1-S8)
		if tissue == "ipsc_l1-s8":
			download_ref("BSS01107", tissue)
		#iPSC (L1-S8R)
		if tissue == "ipsc_l1-s8r":
			download_ref("BSS01108", tissue)
			
		##Kidney
		#Kidney cell (HEK293)
		if tissue == "kidney_hek293":
			download_ref("BSS00526", tissue)
		#Kidney cell (HEK293T)
		if tissue == "kidney_hek293t":
			download_ref("BSS00528", tissue)
		#Kidney cell (kidney embryo (59 days) and female embryo (59 days))
		if tissue == "kidney_59_days_59_days_f":
			download_ref("BSS01078", tissue)
		#Kidney (kidney embryo (80 days))
		if tissue == "kidney_80_days":
			download_ref("BSS01079", tissue)
		#Kidney (kidney female embryo (105 days))
		if tissue == "kidney_105_days_f":
			download_ref("BSS01084", tissue)
		#Kidney (kidney female embryo (108 days))
		if tissue == "kidney_108_days_f":
			download_ref("BSS01085", tissue)
		#Kidney (kidney female embryo (113 days))
		if tissue == "kidney_113_days_f":
			download_ref("BSS01086", tissue)
		#Kidney (kidney female embryo (120 days))
		if tissue == "kidney_120_days_f":
			download_ref("BSS01088", tissue)
		#Kidney (kidney female embryo (121 day))
		if tissue == "kidney_121_days_f":
			download_ref("BSS01089", tissue)
		#Kidney (kidney female embryo (76 days) and male embryo (76 days))
		if tissue == "kidney_76_days_f_76_days_m":
			download_ref("BSS01090", tissue)
		#Kidney (kidney female embryo (85 days))
		if tissue == "kidney_85_days_f":
			download_ref("BSS01091", tissue)
		#Kidney (kidney male adult (50 years))
		if tissue == "kidney_50m":
			download_ref("BSS01096", tissue)
		#Kidney (kidney male adult (67 years))
		if tissue == "kidney_67m":
			download_ref("BSS01097", tissue)
		#Kidney (kidney male embryo (105 days))
		if tissue == "kidney_105_days_m":
			download_ref("BSS01099", tissue)
		#Kidney (kidney male embryo (85 days))
		if tissue == "kidney_85_days_m":
			download_ref("BSS01100", tissue)
		#Kidney (kidney male embryo (87 days))
		if tissue == "kidney_87_days_m":
			download_ref("BSS01101", tissue)
		#Kidney (left kidney female embryo (107 days))
		if tissue == "kidney_left_107_days_f":
			download_ref("BSS01128", tissue)
		#Kidney (left kidney female embryo (110 days))
		if tissue == "kidney_left_110_days_f":
			download_ref("BSS01129", tissue)
		#Kidney (left kidney female embryo (147 days))
		if tissue == "kidney_left_147_days_f":
			download_ref("BSS01130", tissue)
		#Kidney (left kidney female embryo (59 days) and male embryo (91 day))
		if tissue == "kidney_left_59_days_f_91_days_m":
			download_ref("BSS01131", tissue)
		#Kidney (left kidney female embryo (87 days))
		if tissue == "kidney_left_87_days_f":
			download_ref("BSS01132", tissue)
		#Kidney (left kidney female embryo (98 days))
		if tissue == "kidney_left_89_days_f":
			download_ref("BSS01133", tissue)
		#Kidney (left kidney male embryo (115 days))
		if tissue == "kidney_left_115_days_m":
			download_ref("BSS01134", tissue)
		#Kidney (left kidney male embryo (87 days))
		if tissue == "kidney_left_87_days_m":
			download_ref("BSS01135", tissue)
		#Kidney (left kidney male embryo (96 days))
		if tissue == "kidney_left_96_days_m":
			download_ref("BSS01136", tissue)
		#Renal cortex interstitium (left renal cortex interstitium male embryo (105 days))
		if tissue == "kidney_left_renal_cortex_interstitium_105_days_m":
			download_ref("BSS01150", tissue)
		#Renal cortex interstitium (left renal cortex interstitium male embryo (120 days))
		if tissue == "kidney_left_renal_cortex_interstitium_120_days_m":
			download_ref("BSS01151", tissue)
		#Renal pelvis (left renal pelvis male embryo (105 days))
		if tissue == "kidney_left_renal_pelvis_105_days_m":
			download_ref("BSS01152", tissue)
		#Renal pelvis (left renal pelvis male embryo (120 days))
		if tissue == "kidney_left_renal_pelvis_120_days_m":
			download_ref("BSS01153", tissue)
		#Renal cortex interstitium (renal cortex interstitium female embryo (103 days))
		if tissue == "kidney_renal_cortex_interstitium_103_days_f":
			download_ref("BSS01482", tissue)
		#Renal cortex interstitium (renal cortex interstitium female embryo (120 days))
		if tissue == "kidney_renal_cortex_interstitium_120_days_f":
			download_ref("BSS01483", tissue)
		#Renal cortex interstitium (renal cortex interstitium female embryo (89 days))
		if tissue == "kidney_renal_cortex_interstitium_89_days_f":
			download_ref("BSS01484", tissue)
		#Renal cortex interstitium (renal cortex interstitium female embryo (96 days))
		if tissue == "kidney_renal_cortex_interstitium_96_days_f":
			download_ref("BSS01485", tissue)
		#Renal cortex interstitium (renal cortex interstitium male embryo (108 days))
		if tissue == "kidney_renal_cortex_interstitium_108_days_m":
			download_ref("BSS01486", tissue)
		#Renal cortex interstitium (renal cortex interstitium male embryo (113 days))
		if tissue == "kidney_renal_cortex_interstitium_113_days_m":
			download_ref("BSS01487", tissue)
		#Renal cortex interstitium (renal cortex interstitium male embryo (127 days))
		if tissue == "kidney_renal_cortex_interstitium_127_days_m":
			download_ref("BSS01488", tissue)
		#Renal cortex interstitium (renal cortex interstitium male embryo (91 day))
		if tissue == "kidney_renal_cortex_interstitium_91_days_m":
			download_ref("BSS01489", tissue)
		#Renal cortex interstitium (renal cortex interstitium male embryo (97 days))
		if tissue == "kidney_renal_cortex_interstitium_97_days_m":
			download_ref("BSS01490", tissue)
		#Renal pelvis (renal pelvis female embryo (103 days))
		if tissue == "kidney_renal_pelvis_103_days_f":
			download_ref("BSS01493", tissue)
		#Renal pelvis (renal pelvis female embryo (105 days))
		if tissue == "kidney_renal_pelvis_105_days_f":
			download_ref("BSS01494", tissue)
		#Renal pelvis (renal pelvis female embryo (89 days))
		if tissue == "kidney_renal_pelvis_89_days_f":
			download_ref("BSS01495", tissue)
		#Renal pelvis (renal pelvis female embryo (96 days))
		if tissue == "kidney_renal_pelvis_96_days_f":
			download_ref("BSS01496", tissue)
		#Renal pelvis (renal pelvis male embryo (108 days))
		if tissue == "kidney_renal_pelvis_108_days_m":
			download_ref("BSS01497", tissue)
		#Renal pelvis (renal pelvis male embryo (113 days))
		if tissue == "kidney_renal_pelvis_113_days_m":
			download_ref("BSS01498", tissue)
		#Renal pelvis (renal pelvis male embryo (127 days))
		if tissue == "kidney_renal_pelvis_127_days_m":
			download_ref("BSS01499", tissue)
		#Renal pelvis (renal pelvis male embryo (91 day))
		if tissue == "kidney_renal_pelvis_91_days_m":
			download_ref("BSS01500", tissue)
		#Renal pelvis (renal pelvis male embryo (97 days))
		if tissue == "kidney_renal_pelvis_97_days_m":
			download_ref("BSS01501", tissue)
		#Kidney (right kidney female embryo (107 days))
		if tissue == "kidney_right_107_days_f":
			download_ref("BSS01509", tissue)
		#Kidney (right kidney female embryo (117 days))
		if tissue == "kidney_right_117_days_f":
			download_ref("BSS01510", tissue)
		#Kidney (right kidney female embryo (147 days))
		if tissue == "kidney_right_147_days_f":
			download_ref("BSS01511", tissue)
		#Kidney (right kidney female embryo (87 days))
		if tissue == "kidney_right_87_days_f":
			download_ref("BSS01512", tissue)
		#Kidney (right kidney female embryo (98 days))
		if tissue == "kidney_right_98_days_f":
			download_ref("BSS01513", tissue)
		#Kidney (right kidney male embryo (108 days))
		if tissue == "kidney_right_108_days_m":
			download_ref("BSS01514", tissue)
		#Kidney (right kidney male embryo (115 days))
		if tissue == "kidney_right_115_days_m":
			download_ref("BSS01515", tissue)
		#Kidney (right kidney male embryo (87 days))
		if tissue == "kidney_right_87_days_m":
			download_ref("BSS01516", tissue)
		#Kidney (right kidney male embryo (91 day))
		if tissue == "kidney_right_91_days_m":
			download_ref("BSS01517", tissue)
		#Kidney (right kidney male embryo (96 days))
		if tissue == "kidney_right_96_days_m":
			download_ref("BSS01518", tissue)
		#Renal cortex interstitium (right renal cortex interstitium male embryo (105 days))
		if tissue == "kidney_right_renal_cortex_interstitium_105_days_m":
			download_ref("BSS01531", tissue)
		#Renal cortex interstitium (right renal cortex interstitium male embryo (120 days))
		if tissue == "kidney_right_renal_cortex_interstitium_120_days_m":
			download_ref("BSS01532", tissue)
		#Renal pelvis (right renal pelvis male embryo (105 days))
		if tissue == "kidney_right_renal_pelvis_105_days_m":
			download_ref("BSS01533", tissue)
		#Renal pelvis (right renal pelvis male embryo (120 days))
		if tissue == "kidney_right_renal_pelvis_120_days_m":
			download_ref("BSS01534", tissue)
			
		##Liver
		#Hepatic stellate cell (hepatic stellate cell female adult (59 years))
		if tissue == "liver_hepatic_stellate_cell_59f":
			download_ref("BSS00553", tissue)
		#Hepatocyte (hepatocyte)
		if tissue == "liver_hepatocyte":
			download_ref("BSS00554", tissue)
		#Liver (liver embryo (59 days) and embryo (80 days))
		if tissue == "liver_59_days_80_days":
			download_ref("BSS01158", tissue)
		#Liver (liver female adult (25 years))
		if tissue == "liver_25f":
			download_ref("BSS01159", tissue)
		#Liver (liver female embryo (101 day) and female embryo (113 days))
		if tissue == "liver_101_days_113_days_f":
			download_ref("BSS01164", tissue)
		#Liver (liver male adult (31 year))
		if tissue == "liver_31m":
			download_ref("BSS01168", tissue)
		#Liver (liver male adult (32 years))
		if tissue == "liver_32m":
			download_ref("BSS01169", tissue)
		#Liver (liver male adult (78 years))
		if tissue == "liver_78m":
			download_ref("BSS01170", tissue)
		#Liver (right lobe of liver female adult (53 years))
		if tissue == "liver_right_lobe_53f":
			download_ref("BSS01519", tissue)
			
		##Lung
		#Lung (left lung female embryo (105 days))
		if tissue == "lung_left_105_days_f":
			download_ref("BSS01137", tissue)
		#Lung (left lung female embryo (107 days))
		if tissue == "lung_left_107_days_f":
			download_ref("BSS01138", tissue)
		#Lung (left lung female embryo (108 days))
		if tissue == "lung_left_108_days_f":
			download_ref("BSS01139", tissue)
		#Lung (left lung female embryo (110 days))
		if tissue == "lung_left_110_days_f":
			download_ref("BSS01140", tissue)
		#Lung (left lung female embryo (117 days))
		if tissue == "lung_left_117_days_f":
			download_ref("BSS01141", tissue)
		#Lung (left lung female embryo (91 day))
		if tissue == "lung_left_91_days_f":
			download_ref("BSS01142", tissue)
		#Lung (left lung female embryo (98 days))
		if tissue == "lung_left_98_days_f":
			download_ref("BSS01143", tissue)
		#Lung (left lung male embryo (105 days))
		if tissue == "lung_left_105_days_m":
			download_ref("BSS01144", tissue)
		#Lung (left lung male embryo (113 days))
		if tissue == "lung_left_113_days_m":
			download_ref("BSS01145", tissue)
		#Lung (left lung male embryo (115 days))
		if tissue == "lung_left_115_days_m":
			download_ref("BSS01146", tissue)
		#Lung (left lung male embryo (87 days))
		if tissue == "lung_left_87_days_m":
			download_ref("BSS01147", tissue)
		#Lung (left lung male embryo (91 day))
		if tissue == "lung_left_91_days_m":
			download_ref("BSS01148", tissue)
		#Lung (left lung male embryo (96 days))
		if tissue == "lung_left_96_days_m":
			download_ref("BSS01149", tissue)
		#Lung (lung embryo (101 day))
		if tissue == "lung_101_days":
			download_ref("BSS01186", tissue)
		#Lung (lung embryo (112 days))
		if tissue == "lung_112_days":
			download_ref("BSS01187", tissue)
		#Lung (lung embryo (67 days))
		if tissue == "lung_67_days":
			download_ref("BSS01188", tissue)
		#Lung (lung embryo (80 days) and male embryo (76 days))
		if tissue == "lung_80_days_76_days_m":
			download_ref("BSS01189", tissue)
		#Lung (lung female adult (30 years))
		if tissue == "lung_30f":
			download_ref("BSS01190", tissue)
		#Lung (lung female embryo (108 days))
		if tissue == "lung_108_days_f":
			download_ref("BSS01192", tissue)
		#Lung (lung female embryo (120 days))
		if tissue == "lung_120_days_f":
			download_ref("BSS01193", tissue)
		#Lung (lung female embryo (76 days))
		if tissue == "lung_76_days_f":
			download_ref("BSS01195", tissue)
		#Lung (lung female embryo (82 days))
		if tissue == "lung_82_days_f":
			download_ref("BSS01196", tissue)
		#Lung (lung female embryo (85 days))
		if tissue == "lung_85_days_f":
			download_ref("BSS01197", tissue)
		#Lung (lung female embryo (96 days))
		if tissue == "lung_96_days_f":
			download_ref("BSS01198", tissue)
		#Lung (lung male child (3 years))
		if tissue == "lung_3m":
			download_ref("BSS01201", tissue)
		#Lung (lung male embryo (103 days))
		if tissue == "lung_103_days_m":
			download_ref("BSS01202", tissue)
		#Lung (lung male embryo (108 days))
		if tissue == "lung_108_days_m":
			download_ref("BSS01203", tissue)
		#Lung (lung male embryo (54 days) and male embryo (58 days))
		if tissue == "lung_54_days_58_days_m":
			download_ref("BSS01204", tissue)
		#Lung (lung male embryo (82 days))
		if tissue == "lung_82_days_m":
			download_ref("BSS01205", tissue)
		#Lung (right lung female embryo (105 days))
		if tissue == "lung_right_105_days_f":
			download_ref("BSS01520", tissue)
		#Lung (right lung female embryo (107 days))
		if tissue == "lung_right_107_days_f":
			download_ref("BSS01521", tissue)
		#Lung (right lung female embryo (108 days))
		if tissue == "lung_right_108_days_f":
			download_ref("BSS01522", tissue)
		#Lung (right lung female embryo (110 days))
		if tissue == "lung_right_110_days_f":
			download_ref("BSS01523", tissue)
		#Lung (right lung female embryo (117 days))
		if tissue == "lung_right_117_days_f":
			download_ref("BSS01524", tissue)
		#Lung (right lung female embryo (91 day))
		if tissue == "lung_right_91_days_f":
			download_ref("BSS01525", tissue)
		#Lung (right lung female embryo (98 days))
		if tissue == "lung_right_98_days_f":
			download_ref("BSS01526", tissue)
		#Lung (right lung male embryo (105 days))
		if tissue == "lung_right_105_days_m":
			download_ref("BSS01527", tissue)
		#Lung (right lung male embryo (115 days))
		if tissue == "lung_right_115_days_m":
			download_ref("BSS01528", tissue)
		#Lung (right lung male embryo (87 days))
		if tissue == "lung_right_87_days_m":
			download_ref("BSS01529", tissue)
		#Lung (right lung male embryo (96 days))
		if tissue == "lung_right_96_days_m":
			download_ref("BSS01530", tissue)
		#Lung (upper lobe of left lung female adult (51 year))
		if tissue == "lung_left_upper_lobe_51f":
			download_ref("BSS01869", tissue)
		#Lung (upper lobe of left lung female adult (53 years))
		if tissue == "lung_left_upper_lobe_53f":
			download_ref("BSS01870", tissue)
		#Lung (upper lobe of left lung male adult (37 years))
		if tissue == "lung_left_upper_lobe_37m":
			download_ref("BSS01871", tissue)
			
		##Lymphoblastoid
		#Lymphoblastoid cell line (GM06990)
		if tissue == "lymphoblastoid_gm06990":
			download_ref("BSS00395", tissue)
		#Lymphoblastoid cell line (GM08714)
		if tissue == "lymphoblastoid_gm08714":
			download_ref("BSS00403", tissue)
		#Lymphoblastoid cell line (GM10248)
		if tissue == "lymphoblastoid_gm10248":
			download_ref("BSS00404", tissue)
		#Lymphoblastoid cell line (GM10266)
		if tissue == "lymphoblastoid_gm10266":
			download_ref("BSS00405", tissue)
		#Lymphoblastoid cell line (GM12864)
		if tissue == "lymphoblastoid_gm12864":
			download_ref("BSS00427", tissue)
		#Lymphoblastoid cell line (GM12865)
		if tissue == "lymphoblastoid_gm12865":
			download_ref("BSS00428", tissue)
		#Lymphoblastoid cell line (GM12875)
		if tissue == "lymphoblastoid_gm12875":
			download_ref("BSS00438", tissue)
		#Lymphoblastoid cell line (GM12878)
		if tissue == "lymphoblastoid_gm12878":
			download_ref("BSS00439", tissue)
		#Lymphoblastoid cell line (GM12891)
		if tissue == "lymphoblastoid_gm12891":
			download_ref("BSS00452", tissue)
		#Lymphoblastoid cell line (GM12892)
		if tissue == "lymphoblastoid_gm12892":
			download_ref("BSS00454", tissue)
		#Lymphoblastoid cell line (GM18507)
		if tissue == "lymphoblastoid_gm18507":
			download_ref("BSS00462", tissue)
		#Lymphoblastoid cell line (GM19238)
		if tissue == "lymphoblastoid_gm19238":
			download_ref("BSS00471", tissue)
		#Lymphoblastoid cell line (GM19239)
		if tissue == "lymphoblastoid_gm19239":
			download_ref("BSS00472", tissue)
		#Lymphoblastoid cell line (GM19240)
		if tissue == "lymphoblastoid_gm19240":
			download_ref("BSS00473", tissue)
			
		##Mesench
		#Adipocyte from MSC (adipocyte originated from mesenchymal stem cell)
		if tissue == "mesench_adipocyte_msc_origin":
			download_ref("BSS00039", tissue)
		#Amniotic fluid from MSC (dedifferentiated amniotic fluid mesenchymal stem cell)
		if tissue == "mesench_msc_dedifferentiated_amniotic_fluid_origin":
			download_ref("BSS00250", tissue)
		#Embryonic facial prominence (embryonic facial prominence embryo (53 days) and embryo (58 days))
		if tissue == "mesench_embryonic_facial_prominence_53_days_58_days":
			download_ref("BSS00279", tissue)
		#Mesenchymal stem cell (mesenchymal stem cell originated from adipose tissue)
		if tissue == "mesench_msc_adipose_origin":
			download_ref("BSS01260", tissue)
			
		##Muscle
		#Cardiac myocyte (cardiac muscle cell)
		if tissue == "muscle_cardiac_myocyte":
			download_ref("BSS00170", tissue)
		#Arm muscle (forelimb muscle female embryo (108 days))
		if tissue == "muscle_forelimb_108_days_f":
			download_ref("BSS00352", tissue)
		#Gastrocnemius medialis (gastrocnemius medialis female adult (51 year))
		if tissue == "muscle_gastrocnemius_medialis_51f":
			download_ref("BSS00376", tissue)
		#Gastrocnemius medialis (gastrocnemius medialis female adult (53 years))
		if tissue == "muscle_gastrocnemius_medialis_53f":
			download_ref("BSS00377", tissue)
		#Gastrocnemius medialis (gastrocnemius medialis male adult (37 years))
		if tissue == "muscle_gastrocnemius_medialis_37m":
			download_ref("BSS00378", tissue)
		#Gastrocnemius medialis (gastrocnemius medialis male adult (54 years))
		if tissue == "muscle_gastrocnemius_medialis_54m":
			download_ref("BSS00379", tissue)
		#Leg muscle (hindlimb muscle male embryo (120 days))
		if tissue == "muscle_leg_hindlimb_120_days_m":
			download_ref("BSS00700", tissue)
		#Arm muscle (muscle of arm embryo (101 day))
		if tissue == "muscle_arm_101_days":
			download_ref("BSS01289", tissue)
		#Arm muscle (muscle of arm female embryo (105 days))
		if tissue == "muscle_arm_105_days_f":
			download_ref("BSS01290", tissue)
		#Arm muscle (muscle of arm female embryo (115 days))
		if tissue == "muscle_arm_115_days_f":
			download_ref("BSS01291", tissue)
		#Arm muscle (muscle of arm female embryo (120 days))
		if tissue == "muscle_arm_120_days_f":
			download_ref("BSS01292", tissue)
		#Arm muscle (muscle of arm female embryo (85 days))
		if tissue == "muscle_arm_85_days_f":
			download_ref("BSS01293", tissue)
		#Arm muscle (muscle of arm female embryo (98 days))
		if tissue == "muscle_arm_98_days_f":
			download_ref("BSS01294", tissue)
		#Arm muscle (muscle of arm male embryo (101 day))
		if tissue == "muscle_arm_101_days_m":
			download_ref("BSS01295", tissue)
		#Arm muscle (muscle of arm male embryo (104 days))
		if tissue == "muscle_arm_104_days_m":
			download_ref("BSS01296", tissue)
		#Arm muscle (muscle of arm male embryo (105 days))
		if tissue == "muscle_arm_105_days_m":
			download_ref("BSS01297", tissue)
		#Arm muscle (muscle of arm male embryo (113 days))
		if tissue == "muscle_arm_113_days_m":
			download_ref("BSS01298", tissue)
		#Arm muscle (muscle of arm male embryo (115 days))
		if tissue == "muscle_arm_115_days_m":
			download_ref("BSS01299", tissue)
		#Arm muscle (muscle of arm male embryo (120 days))
		if tissue == "muscle_arm_120_days_m":
			download_ref("BSS01300", tissue)
		#Arm muscle (muscle of arm male embryo (127 days))
		if tissue == "muscle_arm_127_days_m":
			download_ref("BSS01301", tissue)
		#Arm muscle (muscle of arm male embryo (96 days))
		if tissue == "muscle_arm_96_days_m":
			download_ref("BSS01303", tissue)
		#Arm muscle (muscle of arm male embryo (97 days))
		if tissue == "muscle_arm_97_days_m":
			download_ref("BSS01304", tissue)
		#Back muscle (muscle of back female embryo (105 days))
		if tissue == "muscle_back_105_days_f":
			download_ref("BSS01305", tissue)
		#Back muscle (muscle of back female embryo (113 days))
		if tissue == "muscle_back_113_days_f":
			download_ref("BSS01306", tissue)
		#Back muscle (muscle of back female embryo (115 days))
		if tissue == "muscle_back_115_days_f":
			download_ref("BSS01307", tissue)
		#Back muscle (muscle of back female embryo (85 days))
		if tissue == "muscle_back_85_days_f":
			download_ref("BSS01308", tissue)
		#Back muscle (muscle of back female embryo (98 days))
		if tissue == "muscle_back_98_days_f":
			download_ref("BSS01309", tissue)
		#Back muscle (muscle of back male embryo (101 day))
		if tissue == "muscle_back_101_days_m":
			download_ref("BSS01310", tissue)
		#Back muscle (muscle of back male embryo (104 days))
		if tissue == "muscle_back_104_days_m":
			download_ref("BSS01311", tissue)
		#Back muscle (muscle of back male embryo (105 days))
		if tissue == "muscle_back_105_days_m":
			download_ref("BSS01312", tissue)
		#Back muscle (muscle of back male embryo (108 days))
		if tissue == "muscle_back_108_days_m":
			download_ref("BSS01313", tissue)
		#Back muscle (muscle of back male embryo (127 days))
		if tissue == "muscle_back_127_days_m":
			download_ref("BSS01314", tissue)
		#Back muscle (muscle of back male embryo (91 day))
		if tissue == "muscle_back_91_days_m":
			download_ref("BSS01315", tissue)
		#Back muscle (muscle of back male embryo (96 days))
		if tissue == "muscle_back_96_days_m":
			download_ref("BSS01316", tissue)
		#Back muscle (muscle of back male embryo (97 days))
		if tissue == "muscle_back_97_days_m":
			download_ref("BSS01317", tissue)
		#Leg muscle (muscle of leg female embryo (105 days))
		if tissue == "muscle_leg_105_days_f":
			download_ref("BSS01318", tissue)
		#Leg muscle (muscle of leg female embryo (110 days))
		if tissue == "muscle_leg_110_days_f":
			download_ref("BSS01319", tissue)
		#Leg muscle (muscle of leg female embryo (113 days))
		if tissue == "muscle_leg_113_days_f":
			download_ref("BSS01320", tissue)
		#Leg muscle (muscle of leg female embryo (115 days))
		if tissue == "muscle_leg_115_days_f":
			download_ref("BSS01321", tissue)
		#Leg muscle (muscle of leg female embryo (85 days))
		if tissue == "muscle_leg_85_days_f":
			download_ref("BSS01322", tissue)
		#Leg muscle (muscle of leg male embryo (101 day))
		if tissue == "muscle_leg_101_days_m":
			download_ref("BSS01323", tissue)
		#Leg muscle (muscle of leg male embryo (104 days))
		if tissue == "muscle_leg_104_days_m":
			download_ref("BSS01324", tissue)
		#Leg muscle (muscle of leg male embryo (105 days))
		if tissue == "muscle_leg_105_days_m":
			download_ref("BSS01325", tissue)
		#Leg muscle (muscle of leg male embryo (115 days))
		if tissue == "muscle_leg_115_days_m":
			download_ref("BSS01327", tissue)
		#Leg muscle (muscle of leg male embryo (127 days))
		if tissue == "muscle_leg_127_days_m":
			download_ref("BSS01328", tissue)
		#Leg muscle (muscle of leg male embryo (96 days))
		if tissue == "muscle_leg_96_days_m":
			download_ref("BSS01329", tissue)
		#Leg muscle (muscle of leg male embryo (97 days))
		if tissue == "muscle_leg_97_days_m":
			download_ref("BSS01330", tissue)
		#Trunk muscle (muscle of trunk female embryo (113 days))
		if tissue == "muscle_trunk_113_days_f":
			download_ref("BSS01331", tissue)
		#Trunk muscle (muscle of trunk female embryo (115 days))
		if tissue == "muscle_trunk_115_days_f":
			download_ref("BSS01332", tissue)
		#Trunk muscle (muscle of trunk female embryo (120 days))
		if tissue == "muscle_trunk_120_days_f":
			download_ref("BSS01333", tissue)
		#Trunk muscle (muscle of trunk female embryo (121 day))
		if tissue == "muscle_trunk_121_days_f":
			download_ref("BSS01334", tissue)
		#Psoas muscle (psoas muscle female adult (30 years))
		if tissue == "muscle_psoas_30f":
			download_ref("BSS01460", tissue)
		#Psoas muscle (psoas muscle male adult (27 years) and male adult (35 years))
		if tissue == "muscle_psoas_27m_35m":
			download_ref("BSS01461", tissue)
		#Psoas muscle (psoas muscle male adult (34 years))
		if tissue == "muscle_psoas_34m":
			download_ref("BSS01462", tissue)
		#Psoas muscle (psoas muscle male child (3 years))
		if tissue == "muscle_psoas_3m":
			download_ref("BSS01463", tissue)
		#Skeletal muscle cell (skeletal muscle cell)
		if tissue == "muscle_skeletal_cell":
			download_ref("BSS01572", tissue)
		#Skeletal muscle (skeletal muscle tissue)
		if tissue == "muscle_skeletal_tissue":
			download_ref("BSS01577", tissue)
		#Skeletal muscle (skeletal muscle tissue female adult (72 years))
		if tissue == "muscle_skeletal_tissue_72f":
			download_ref("BSS01578", tissue)
		#Skeletal muscle (skeletal muscle tissue male adult (54 years))
		if tissue == "muscle_skeletal_tissue_54m":
			download_ref("BSS01581", tissue)
		#Tongue (tongue female embryo (59 days) and female embryo (76 days))
		if tissue == "muscle_tongue_59_days_f_76_days_f":
			download_ref("BSS01845", tissue)
		#Tongue (tongue male embryo (72 days))
		if tissue == "muscle_tongue_72_days_m":
			download_ref("BSS01846", tissue)
			
		##Myosat
		#Skeletal muscle myoblast (LHCN-M2)
		if tissue == "myosat_skeletal_muscle_myoblast_lhcn-m2":
			download_ref("BSS01155", tissue)
		#Myocyte (myocyte originated from LHCN-M2)
		if tissue == "myosat_myocyte_lhcn-m2_origin":
			download_ref("BSS01338", tissue)
		#Myotube (myotube originated from skeletal muscle myoblast)
		if tissue == "myosat_myotube_skeletal_muscle_myoblast_origin":
			download_ref("BSS01344", tissue)
		#Skeletal muscle myoblast (skeletal muscle myoblast)
		if tissue == "myosat_skeletal_muscle_myoblast":
			download_ref("BSS01573", tissue)
		#Skeletal muscle myoblast (skeletal muscle myoblast male adult (22 years))
		if tissue == "myosat_skeletal_muscle_myoblast_22m":
			download_ref("BSS01574", tissue)
		#Skeletal muscle satellite cell (skeletal muscle satellite cell female adult originated from mesodermal cell)
		if tissue == "myosat_skeletal_muscle_satellite_cell_mesoderm_origin_f":
			download_ref("BSS01576", tissue)
			
		##Neurosph
		#Neurosphere (neurosphere embryo (15 weeks) originated from ganglionic eminence)
		if tissue == "neurosph_15_weeks_ganglionic_eminence_origin":
			download_ref("BSS01377", tissue)
		#Neurosphere (neurosphere female embryo (17 weeks) originated from cortex)
		if tissue == "neurosph_17_weeks_f_ganglionic_cortex_origin":
			download_ref("BSS01378", tissue)
		#Neurosphere (neurosphere female embryo (17 weeks) originated from ganglionic eminence)
		if tissue == "neurosph_17_weeks_f_ganglionic_eminence_origin":
			download_ref("BSS01379", tissue)
		#Olfactory neurosphere (olfactory neurosphere cell line)
		if tissue == "neurosph_olfactory_cell_line":
			download_ref("BSS01392", tissue)
		
		##Other
		#Breast epithelium (breast epithelium female adult (51 year))
		if tissue == "other_breast_epithelium_51f":
			download_ref("BSS00145", tissue)
		#Breast epithelium (breast epithelium female adult (53 years))
		if tissue == "other_breast_epithelium_53f":
			download_ref("BSS00146", tissue)
		#Epidermal melanocyte (epidermal melanocyte)
		if tissue == "other_epidermal_melanocyte":
			download_ref("BSS00304", tissue)
		#Foreskin melanocyte (foreskin melanocyte male newborn)
		if tissue == "other_foreskin_melanocyte_newborn_m":
			download_ref("BSS00368", tissue)
		#Limb embryo (limb embryo (53 days) and embryo (56 days))
		if tissue == "other_limb_embryo_53_days_56_days":
			download_ref("BSS01156", tissue)
		#Limb embryo (limb embryo (58 days) and embryo (59 days))
		if tissue == "other_limb_embryo_58_days_59_days":
			download_ref("BSS01157", tissue)
		#Mammary stem cell (mammary stem cell)
		if tissue == "other_mammary_stem_cell":
			download_ref("BSS01216", tissue)
			
		##Pancreas
		#Body of pancreas (body of pancreas female adult (51 year))
		if tissue == "pancreas_body_51f":
			download_ref("BSS00121", tissue)
		#Body of pancreas (body of pancreas female adult (53 years))
		if tissue == "pancreas_body_53f":
			download_ref("BSS00122", tissue)
		#Body of pancreas (body of pancreas male adult (37 years))
		if tissue == "pancreas_body_37m":
			download_ref("BSS00123", tissue)
		#Body of pancreas (body of pancreas male adult (54 years))
		if tissue == "pancreas_body_54m":
			download_ref("BSS00124", tissue)
		#Islet precursor cell (islet precursor cell)
		if tissue == "pancreas_islet_precursor_cell":
			download_ref("BSS00758", tissue)
		#Pancreas (pancreas female adult (30 years))
		if tissue == "pancreas_30f":
			download_ref("BSS01406", tissue)
		#Pancreas (pancreas male adult (34 years))
		if tissue == "pancreas_34m":
			download_ref("BSS01407", tissue)
			
		##Placenta & EEM
		#Amnion (amnion male embryo (16 weeks))
		if tissue == "placenta_and_eem_amnion_16_weeks_m":
			download_ref("BSS00074", tissue)
		#Amnion stem cell (amniotic stem cell)
		if tissue == "placenta_and_eem_amnion_stem_cell":
			download_ref("BSS00076", tissue)
		#Chorion (chorion)
		if tissue == "placenta_and_eem_chorion":
			download_ref("BSS00209", tissue)
		#Chorion (chorion female embryo (40 weeks))
		if tissue == "placenta_and_eem_chorion_40_weeks_f":
			download_ref("BSS00211", tissue)
		#Chorion (chorion male embryo (16 weeks))
		if tissue == "placenta_and_eem_chorion_16_weeks_m":
			download_ref("BSS00212", tissue)
		#Chorionic villus (chorionic villus embryo (16 weeks))
		if tissue == "placenta_and_eem_chorionic_villus_16_weeks":
			download_ref("BSS00214", tissue)
		#Chorionic villus (chorionic villus female embryo (40 weeks))
		if tissue == "placenta_and_eem_chorionic_villus_40_weeks_f":
			download_ref("BSS00215", tissue)
		#Chorionic villus (chorionic villus male embryo (16 weeks))
		if tissue == "placenta_and_eem_chorionic_villus_16_weeks_m":
			download_ref("BSS00216", tissue)
		#Chorionic villus (chorionic villus male embryo (38 weeks))
		if tissue == "placenta_and_eem_chorionic_villus_38_weeks_m":
			download_ref("BSS00217", tissue)
		#Trophoblast (HTR-8 SVneo)
		if tissue == "placenta_and_eem_trophoblast_htr-8_svneo":
			download_ref("BSS00714", tissue)
		#Placenta (placenta embryo (102 days))
		if tissue == "placenta_and_eem_placenta_102_days":
			download_ref("BSS01430", tissue)
		#Placenta (placenta embryo (16 weeks))
		if tissue == "placenta_and_eem_placenta_16_weeks":
			download_ref("BSS01431", tissue)
		#Placenta (placenta embryo (53 days))
		if tissue == "placenta_and_eem_placenta_53_days":
			download_ref("BSS01432", tissue)
		#Placenta (placenta embryo (56 days) and embryo (59 days))
		if tissue == "placenta_and_eem_placenta_56_days_59_days":
			download_ref("BSS01433", tissue)
		#Placenta (placenta female embryo (101 day) and male embryo (105 days))
		if tissue == "placenta_and_eem_placenta_101_days_f_105_days_m":
			download_ref("BSS01435", tissue)
		#Placenta (placenta female embryo (105 days))
		if tissue == "placenta_and_eem_placenta_105_days_f":
			download_ref("BSS01436", tissue)
		#Placenta (placenta female embryo (108 days))
		if tissue == "placenta_and_eem_placenta_108_days_f":
			download_ref("BSS01437", tissue)
		#Placenta (placenta female embryo (113 days))
		if tissue == "placenta_and_eem_placenta_113_days_f":
			download_ref("BSS01438", tissue)
		#Placenta (placenta female embryo (85 days))
		if tissue == "placenta_and_eem_placenta_85_days_f":
			download_ref("BSS01440", tissue)
		#Placenta (placenta male embryo (16 weeks))
		if tissue == "placenta_and_eem_placenta_16_weeks_m":
			download_ref("BSS01441", tissue)
		#Placenta (placenta male embryo (85 days))
		if tissue == "placenta_and_eem_placenta_85_days_m":
			download_ref("BSS01443", tissue)
		#Placenta (placenta male embryo (91 day))
		if tissue == "placenta_and_eem_placenta_91_days_m":
			download_ref("BSS01444", tissue)
		#Placenta (placental basal plate female embryo (40 weeks))
		if tissue == "placenta_and_eem_placental_basal_plate_40_weeks_f":
			download_ref("BSS01446", tissue)
		#Placenta (placental basal plate male embryo (38 weeks))
		if tissue == "placenta_and_eem_placental_basal_plate_38_weeks_m":
			download_ref("BSS01448", tissue)
		#Trophoblast (trophoblast cell embryo (17 weeks) and embryo (18 weeks))
		if tissue == "placenta_and_eem_trophoblast_cell_17_weeks_18_weeks":
			download_ref("BSS01852", tissue)
		#Trophoblast (trophoblast cell embryo (21 week))
		if tissue == "placenta_and_eem_trophoblast_cell_21_weeks":
			download_ref("BSS01853", tissue)
		#Trophoblast (trophoblast cell embryo (23 weeks))
		if tissue == "placenta_and_eem_trophoblast_cell_23_weeks":
			download_ref("BSS01855", tissue)
		#Trophoblast (trophoblast cell embryo (39 weeks) and embryo (40 weeks))
		if tissue == "placenta_and_eem_trophoblast_cell_39_weeks_40_weeks":
			download_ref("BSS01856", tissue)
		#Trophoblast (trophoblast female embryo (20 weeks))
		if tissue == "placenta_and_eem_trophoblast_20_weeks_f":
			download_ref("BSS01859", tissue)
		#Trophoblast (trophoblast female embryo (40 weeks))
		if tissue == "placenta_and_eem_trophoblast_40_weeks_f":
			download_ref("BSS01860", tissue)
		#Trophoblast (umbilical cord embryo (59 days) and male embryo (76 days))
		if tissue == "placenta_and_eem_umbilical_cord_59_days_76_days_m":
			download_ref("BSS01867", tissue)
			
		##PNS
		#Spinal cord (spinal cord female embryo (108 days))
		if tissue == "pns_spinal_cord_108_days_f":
			download_ref("BSS01613", tissue)
		#Spinal cord (spinal cord female embryo (113 days))
		if tissue == "pns_spinal_cord_113_days_f":
			download_ref("BSS01614", tissue)
		#Spinal cord (spinal cord female embryo (59 days) and male embryo (72 days))
		if tissue == "pns_spinal_cord_59_days_f_72_days_m":
			download_ref("BSS01617", tissue)
		#Spinal cord (spinal cord female embryo (87 days))
		if tissue == "pns_spinal_cord_87_days_f":
			download_ref("BSS01618", tissue)
		#Spinal cord (spinal cord female embryo (89 days))
		if tissue == "pns_spinal_cord_89_days_f":
			download_ref("BSS01619", tissue)
		#Spinal cord (spinal cord male embryo (105 days))
		if tissue == "pns_spinal_cord_105_days_m":
			download_ref("BSS01620", tissue)
		#Spinal cord (spinal cord male embryo (96 days))
		if tissue == "pns_spinal_cord_96_days_m":
			download_ref("BSS01621", tissue)
		#Tibial nerve (tibial nerve female adult (51 year))
		if tissue == "pns_tibial_nerve_51f":
			download_ref("BSS01840", tissue)
		#Tibial nerve (tibial nerve female adult (53 years))
		if tissue == "pns_tibial_nerve_53f":
			download_ref("BSS01841", tissue)
		#Tibial nerve (tibial nerve male adult (37 years))
		if tissue == "pns_tibial_nerve_37m":
			download_ref("BSS01842", tissue)
			
		##Reproductive
		#Prostate gland (prostate gland male adult (37 years))
		if tissue == "reproductive_prostate_gland_37m":
			download_ref("BSS01456", tissue)
		#Prostate gland (prostate male adult (54 years))
		if tissue == "reproductive_prostate_gland_54m":
			download_ref("BSS01459", tissue)
		#Uterus (uterus female adult (53 years))
		if tissue == "reproductive_uterus_53f":
			download_ref("BSS01884", tissue)
		#Uterus (vagina female adult (51 year))
		if tissue == "reproductive_vagina_51f":
			download_ref("BSS01886", tissue)
		#Uterus (vagina female adult (53 years))
		if tissue == "reproductive_vagina_53f":
			download_ref("BSS01887", tissue)
			
		##Sm. Muscle
		#Colon muscle (muscle layer of colon female adult (56 years))
		if tissue == "sm_muscle_colon_56f":
			download_ref("BSS01285", tissue)
		#Colon muscle (muscle layer of colon female adult (77 years))
		if tissue == "sm_muscle_colon_77f":
			download_ref("BSS01286", tissue)
		#Duodenum muscle (muscle layer of duodenum male adult (59 years))
		if tissue == "sm_muscle_duodenum_59m":
			download_ref("BSS01287", tissue)
		#Duodenum muscle (muscle layer of duodenum male adult (73 years))
		if tissue == "sm_muscle_duodenum_73m":
			download_ref("BSS01288", tissue)
		#Rectum muscle (rectal smooth muscle tissue female adult (50 years))
		if tissue == "sm_muscle_rectum_50f":
			download_ref("BSS01475", tissue)
		#Brain vasculature smooth cell (smooth muscle cell of the brain vasculature female)
		if tissue == "sm_muscle_brain_vasculature_smooth_cell":
			download_ref("BSS01606", tissue)
		#Stomach muscle (stomach smooth muscle female adult (84 years))
		if tissue == "sm_muscle_stomach_84f":
			download_ref("BSS01659", tissue)
		#Stomach muscle (stomach smooth muscle male adult (59 years))
		if tissue == "sm_muscle_stomach_59m":
			download_ref("BSS01660", tissue)
			
		##Spleen
		#Spleen (spleen embryo (112 days))
		if tissue == "spleen_112_days":
			download_ref("BSS01625", tissue)
		#Spleen (spleen female adult (30 years))
		if tissue == "spleen_30f":
			download_ref("BSS01628", tissue)
		#Spleen (spleen female adult (53 years))
		if tissue == "spleen_53f":
			download_ref("BSS01630", tissue)
		#Spleen (spleen male adult (34 years))
		if tissue == "spleen_34m":
			download_ref("BSS01631", tissue)
		#Spleen (spleen male adult (54 years))
		if tissue == "spleen_54m":
			download_ref("BSS01633", tissue)
		#Spleen (spleen male child (3 years))
		if tissue == "spleen_3m":
			download_ref("BSS01634", tissue)
			
		##Stromal
		#Skin fibroblast (AG04449)
		if tissue == "stromal_skin_fibroblast_ag04449":
			download_ref("BSS00061", tissue)
		#Lung fibroblast (AG04450)
		if tissue == "stromal_lung_fibroblast_ag04450":
			download_ref("BSS00062", tissue)
		#Skin fibroblast (AG08395)
		if tissue == "stromal_skin_fibroblast_ag08395":
			download_ref("BSS00063", tissue)
		#Lung fibroblast (AG08396)
		if tissue == "stromal_lung_fibroblast_ag08396":
			download_ref("BSS00064", tissue)
		#Skin fibroblast (AG09309)
		if tissue == "stromal_skin_fibroblast_ag08396":
			download_ref("BSS00066", tissue)
		#Gingival fibroblast (AG09319)
		if tissue == "stromal_gingival_fibroblast_ag09319":
			download_ref("BSS00067", tissue)
		#Skin fibroblast (AG10803)
		if tissue == "stromal_skin_fibroblast_ag10803":
			download_ref("BSS00068", tissue)
		#Skin fibroblast (AG20443)
		if tissue == "stromal_skin_fibroblast_ag20443":
			download_ref("BSS00069", tissue)
		#Skin fibroblast (BJ)
		if tissue == "stromal_skin_fibroblast_bj":
			download_ref("BSS00113", tissue)
		#Pericyte (brain pericyte)
		if tissue == "stromal_brain_pericyte":
			download_ref("BSS00144", tissue)
		#Cardiac fibroblast (cardiac fibroblast)
		if tissue == "stromal_cardiac_fibroblast":
			download_ref("BSS00166", tissue)
		#Cardiac fibroblast (cardiac fibroblast female)
		if tissue == "stromal_cardiac_fibroblast_f":
			download_ref("BSS00167", tissue)
		#Cardiac fibroblast (cardiac fibroblast female embryo (94 days) and female embryo (98 days))
		if tissue == "stromal_cardiac_fibroblast_94_days_f_98_days_f":
			download_ref("BSS00168", tissue)
		#Skin fibroblast (EH)
		if tissue == "stromal_skin_fibroblast_eh":
			download_ref("BSS00275", tissue)
		#Skin fibroblast (EL)
		if tissue == "stromal_skin_fibroblast_el":
			download_ref("BSS00276", tissue)
		#Skin fibroblast (ELR)
		if tissue == "stromal_skin_fibroblast_elr":
			download_ref("BSS00278", tissue)
		#Breast fibroblast (fibroblast of breast female adult (17 years))
		if tissue == "stromal_breast_fibroblast_17f":
			download_ref("BSS00332", tissue)
		#Breast fibroblast (fibroblast of breast female adult (26 years))
		if tissue == "stromal_breast_fibroblast_26f":
			download_ref("BSS00333", tissue)
		#Dermis fibroblast (fibroblast of dermis)
		if tissue == "stromal_dermis_fibroblast":
			download_ref("BSS00334", tissue)
		#Dermis fibroblast (fibroblast of dermis female adult)
		if tissue == "stromal_dermis_fibroblast_f":
			download_ref("BSS00335", tissue)
		#Dermis fibroblast (fibroblast of dermis NONE and female adult)
		if tissue == "stromal_dermis_fibroblast_none_f":
			download_ref("BSS00337", tissue)
		#Gingival fibroblast (fibroblast of gingiva)
		if tissue == "stromal_gingival_fibroblast":
			download_ref("BSS00338", tissue)
		#Lung fibroblast (fibroblast of lung)
		if tissue == "stromal_lung_fibroblast":
			download_ref("BSS00339", tissue)
		#Lung fibroblast (fibroblast of lung female child (11 year) and male adult (45 years))
		if tissue == "stromal_lung_fibroblast_11f_45m":
			download_ref("BSS00341", tissue)
		#Lung fibroblast (fibroblast of lung male adult (45 years))
		if tissue == "stromal_lung_fibroblast_45m":
			download_ref("BSS00342", tissue)
		#Mammary fibroblast (fibroblast of mammary gland female)
		if tissue == "stromal_mammary_fibroblast_f":
			download_ref("BSS00343", tissue)
		#Peridontal ligament fibroblast (fibroblast of peridontal ligament male)
		if tissue == "stromal_peridontal_ligament_fibroblast_m":
			download_ref("BSS00344", tissue)
		#Pulmonary artery fibroblast (fibroblast of pulmonary artery)
		if tissue == "stromal_pulmonary_artery_fibroblast":
			download_ref("BSS00345", tissue)
		#Skin fibroblast (fibroblast of skin of abdomen male embryo (97 days))
		if tissue == "stromal_skin_fibroblast_abdomen_97_days_m":
			download_ref("BSS00346", tissue)
		#Aorta fibroblast (fibroblast of the aortic adventitia female)
		if tissue == "stromal_aorta_fibroblast_f":
			download_ref("BSS00347", tissue)
		#Conjunctiva fibroblast (fibroblast of the conjunctiva)
		if tissue == "stromal_conjunctiva_fibroblast":
			download_ref("BSS00349", tissue)
		#Villous mesenchyme fibroblast (fibroblast of villous mesenchyme)
		if tissue == "stromal_villous_mesenchyme_fibroblast":
			download_ref("BSS00350", tissue)
		#Foreskin fibroblast (foreskin fibroblast male newborn)
		if tissue == "stromal_foreskin_fibroblast_newborn_m":
			download_ref("BSS00353", tissue)
		#Skin fibroblast (GM03348)
		if tissue == "stromal_skin_fibroblast_gm03348":
			download_ref("BSS00390", tissue)
		#Skin fibroblast (GM04503)
		if tissue == "stromal_skin_fibroblast_gm04503":
			download_ref("BSS00393", tissue)
		#Skin fibroblast (GM04504)
		if tissue == "stromal_skin_fibroblast_gm04504":
			download_ref("BSS00394", tissue)
		#Skin fibroblast (GM23248)
		if tissue == "stromal_skin_fibroblast_gm23248":
			download_ref("BSS00476", tissue)
		#Foreskin fibroblast (HFF-Myc originated from foreskin fibroblast)
		if tissue == "stromal_foreskin_fibroblast_hff-myc_foreskin_fibroblast_origin":
			download_ref("BSS00697", tissue)
		#Lung fibroblast (IMR-90)
		if tissue == "stromal_lung_fibroblast_imr-90":
			download_ref("BSS00720", tissue)
		#Skin fibroblast (skin fibroblast male embryo (97 days))
		if tissue == "stromal_lung_fibroblast_97_days_m":
			download_ref("BSS01583", tissue)
		#Bone marrow stromal cell (stromal cell of bone marrow male)
		if tissue == "stromal_bone_marrow_m":
			download_ref("BSS01661", tissue)
		#Lung fibroblast (WI38)
		if tissue == "stromal_lung_fibroblast_wi38":
			download_ref("BSS01891", tissue)
			
		##Thymus
		#Thymus (thymus female embryo)
		if tissue == "thymus_embryo_f":
			download_ref("BSS01818", tissue)
		#Thymus (thymus female embryo (105 days))
		if tissue == "thymus_105_days_f":
			download_ref("BSS01819", tissue)
		#Thymus (thymus female embryo (110 days))
		if tissue == "thymus_110_days_f":
			download_ref("BSS01820", tissue)
		#Thymus (thymus female embryo (113 days))
		if tissue == "thymus_113_days_f":
			download_ref("BSS01821", tissue)
		#Thymus (thymus female embryo (147 days))
		if tissue == "thymus_147_days_f":
			download_ref("BSS01823", tissue)
		#Thymus (thymus female embryo (98 days))
		if tissue == "thymus_98_days_f":
			download_ref("BSS01824", tissue)
		#Thymus (thymus male child (3 years))
		if tissue == "thymus_3m":
			download_ref("BSS01825", tissue)
		#Thymus (thymus male embryo (104 days))
		if tissue == "thymus_104_days_m":
			download_ref("BSS01826", tissue)
		#Thymus (thymus male embryo (108 days))
		if tissue == "thymus_108_days_m":
			download_ref("BSS01827", tissue)
		#Thymus (thymus male embryo (113 days))
		if tissue == "thymus_113_days_m":
			download_ref("BSS01828", tissue)
		#Thymus (thymus male embryo (127 days))
		if tissue == "thymus_127_days_m":
			download_ref("BSS01829", tissue)
			
		##Urinary
		#Urinary bladder (urinary bladder male adult (34 years))
		if tissue == "urinary_bladder_34m":
			download_ref("BSS01876", tissue)
		#Urinary bladder (urinary bladder male embryo (76 days))
		if tissue == "urinary_bladder_76_days_m":
			download_ref("BSS01878", tissue)
		#Urothelium cell (urothelium cell line)
		if tissue == "urinary_urothelium_cell_line":
			download_ref("BSS01879", tissue)
		
		
		#User-provided files
		elif tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_hg19.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_hg19.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_hg19.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_hg19.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_hg19.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_hg19.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_hg19.bed.gz"))
			user_files_index += 1
		
###Mouse###
#GRCm39/mm39
elif args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39":
	#This TSS dataset was manually lifted over from the original mm10 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm39.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm39.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm39.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm39.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm39.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm39.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm39.bed.gz"))
			user_files_index += 1

#GRCm38/mm10
elif args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38":
	if not os.path.exists(os.path.join(args.ref_dir, "refTSS_v3.1_mouse_coordinate.mm10.bed.gz")):
		download_ref_genome("http://reftss.clst.riken.jp/datafiles/current/mouse/refTSS_v3.1_mouse_coordinate.mm10.bed.gz", "refTSS_v3.1_mouse_coordinate.mm10.bed.gz")
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm10.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm10.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm10.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm10.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm10.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm10.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm10.bed.gz"))
			user_files_index += 1
		
	
#GRCm37/mm9
elif args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "grcm37":
	#This TSS dataset was manually lifted over from the original mm10 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm9.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm9.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm9.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm9.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm9.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm9.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm9.bed.gz"))
			user_files_index += 1
	
###Import all of the bed files to be worked on (input and histone mark ChIP-seq) using pybedtools
##Human
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	bed_ref = "hg38"
	bed_tss_path = os.path.join(args.ref_dir, "refTSS_v3.1_human_coordinate.hg38.bed.gz")
elif args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	bed_ref = "hg19"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_human_coordinate.lifted.sorted.hg19.bed.gz"

##Mouse
if args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39":
	bed_ref = "mm39"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_mouse_coordinate.lifted.sorted.mm39.bed.gz"
elif args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38":
	bed_ref = "mm10"
	bed_tss_path = os.path,join(args.ref_dir, "refTSS_v3.1_mouse_coordinate.mm10.bed.gz")
elif args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "mm9":
	bed_ref = "mm9"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_mouse_coordinate.lifted.sorted.mm9.bed.gz"