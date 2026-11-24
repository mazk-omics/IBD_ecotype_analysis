# 本文件固定服务器正式结果和 EcoTyper 原始模型路径；不自动搜索相似目录。
# 只保存不含凭据的路径和模型身份，正式运行前逐项核对。

project_root <- "/home/mazekai/IBD_EcoTyper"
formal_root <- file.path(project_root, "04_CIBERSORTx/MSCCR_SignatureD_TPM_EcoTyper_formal_20261101")
formal_export_root <- file.path(formal_root, "05_ecotyper_output")

dataset_name <- "MSCCR_SIGNATURED_IBD2025_TPM_REGION3_20261101"
fractions_name <- "Custom"
ecotyper_root <- file.path(project_root, "software/ecotyper")
model_root <- file.path(ecotyper_root, "EcoTyper", dataset_name, fractions_name)
rank_data_file <- file.path(model_root, "Analysis/rank_selection/rank_data.txt")
state_discovery_root <- file.path(model_root, "Cell_States/discovery")
ecotype_discovery_root <- file.path(model_root, "Ecotypes/discovery")
advanced_qc_root <- file.path(model_root, "Analysis/advanced_QC_filter")
hires_root <- file.path(ecotyper_root, "CIBERSORTx/hires", dataset_name, fractions_name)

status_file <- file.path(formal_root, "status.txt")
runtime_config_file <- file.path(formal_root, "config_discovery_bulk.retry2.resume.sanitized.yml")
initial_config_file <- file.path(formal_root, "config_discovery_bulk.sanitized.yml")
fractions_file <- file.path(formal_root, "02_fractions/output/CIBERSORTx_Results.txt")
cell_types_file <- file.path(formal_root, "cell_types.txt")
annotation_file <- file.path(project_root, "03_reference/MSCCR_SignatureD_IBD2025_TPM_20261101/MSCCR_CD_UC_2025_annotation.tsv")
mixture_tpm_file <- file.path(project_root, "03_reference/MSCCR_SignatureD_IBD2025_TPM_20261101/MSCCR_CD_UC_2025_Ensembl75_protein_coding_TPM.txt")

scrna_validation_root <- file.path(project_root, "06_validation/MSCCR_SignatureD_HiRes_expression_purification_20261109")
scrna_round2_root <- file.path(project_root, "06_validation/Cell2021_signature_matrix_validation_20261025/06_multidataset_validation_20261031")
eligible_pool_root <- file.path(scrna_round2_root, "01_sampling_design/eligible_pools")
qualification_file <- file.path(RUN_ROOT, "inputs/dataset_donor_sample_eligibility.tsv")
ensembl75_identity_file <- file.path(project_root, "03_reference/web_reference_three_pilots/msccr_gene_identity_20261006/MSCCR_Ensembl75_HGNC_identity.tsv")

expected_samples <- 2025L
expected_cell_types <- 20L
output_root <- file.path(RUN_ROOT, "output")
