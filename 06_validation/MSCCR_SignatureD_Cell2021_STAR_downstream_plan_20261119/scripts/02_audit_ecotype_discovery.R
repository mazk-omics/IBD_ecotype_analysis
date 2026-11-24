#!/usr/bin/env Rscript
# 本脚本只读复核正式模型的 Jaccard、ecotype 分区、abundance 和 assignment。
# 官方 1-Jaccard 为主结果；STAR Euclidean 和 k-medoids 只作为敏感性，不改正式模型。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取冻结模型和前序审计
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 02_audit_ecotype_discovery.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

required_pkgs <- c("cluster", "data.table", "digest", "jsonlite")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)

freeze_file <- file.path(output_root, "00_freeze/formal_model/formal_freeze_summary.json")
state_qc_file <- file.path(output_root, "01_state_qc_audit/state_qc_impact_summary.json")
if (!file.exists(freeze_file) || !file.exists(state_qc_file)) stop("Run Phase 0 and Phase 1 first", call. = FALSE)
freeze <- jsonlite::read_json(freeze_file, simplifyVector = TRUE)
state_qc <- jsonlite::read_json(state_qc_file, simplifyVector = TRUE)
if (!identical(freeze$status, "PASS_FORMAL_IDENTITY")) stop("Formal input freeze did not pass", call. = FALSE)
if (!identical(state_qc$status, "PASS_READ_ONLY_QC_AUDIT")) stop("State-QC audit did not pass", call. = FALSE)

audit_dir <- file.path(output_root, "02_ecotype_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

jaccard_file <- file.path(ecotype_discovery_root, "jaccard_matrix.txt")
binary_file <- file.path(ecotype_discovery_root, "binary_classification_all_states.txt")
initial_ecotype_file <- file.path(ecotype_discovery_root, "initial_ecotypes.txt")
ecotype_file <- file.path(ecotype_discovery_root, "ecotypes.txt")
saved_state_abundance_file <- file.path(ecotype_discovery_root, "combined_state_abundances.txt")
saved_ecotype_abundance_file <- file.path(ecotype_discovery_root, "ecotype_abundance.txt")
saved_assignment_p_file <- file.path(ecotype_discovery_root, "assignment_p_vals.txt")
saved_initial_assignment_file <- file.path(ecotype_discovery_root, "initial_ecotype_assignment.txt")
required_files <- c(jaccard_file, binary_file, initial_ecotype_file, ecotype_file, saved_state_abundance_file, saved_ecotype_abundance_file, saved_assignment_p_file, saved_initial_assignment_file)
if (any(!file.exists(required_files))) stop("Missing ecotype audit input: ", paste(required_files[!file.exists(required_files)], collapse = "; "), call. = FALSE)

# Step2. 从 binary classification 重算 85×85 Jaccard
binary <- read.delim(binary_file, row.names = 1, check.names = FALSE)
clusters <- as.matrix(binary)
if (!is.logical(clusters)) clusters <- clusters == "TRUE"
if (anyNA(clusters)) stop("Binary state classification contains NA", call. = FALSE)

state_counts <- rowSums(clusters)
intersection <- tcrossprod(clusters * 1L)
union <- outer(state_counts, state_counts, "+") - intersection
jaccard_recalculated <- intersection / union

p_hyper <- matrix(NA_real_, nrow(clusters), nrow(clusters), dimnames = list(rownames(clusters), rownames(clusters)))
for (i in seq_len(nrow(clusters))) {
  for (j in i:nrow(clusters)) {
    p_value <- 1 - phyper(intersection[i, j], state_counts[[i]], ncol(clusters) - state_counts[[i]], state_counts[[j]])
    p_hyper[i, j] <- p_hyper[j, i] <- p_value
  }
}
jaccard_recalculated[p_hyper > 0.01] <- 0
dimnames(jaccard_recalculated) <- list(rownames(clusters), rownames(clusters))

jaccard_saved <- as.matrix(read.delim(jaccard_file, row.names = 1, check.names = FALSE))
storage.mode(jaccard_saved) <- "double"
if (!identical(dimnames(jaccard_recalculated), dimnames(jaccard_saved))) stop("Saved Jaccard state order changed", call. = FALSE)
jaccard_max_error <- max(abs(jaccard_recalculated - jaccard_saved))
jaccard_exact <- is.finite(jaccard_max_error) && jaccard_max_error < 1e-12
if (!jaccard_exact) stop("Jaccard matrix does not reproduce; max error = ", signif(jaccard_max_error, 6), call. = FALSE)
cat("JACCARD_REPRODUCED", nrow(jaccard_saved), "states", ncol(clusters), "samples; max error", format(jaccard_max_error, scientific = TRUE), "\n")

# Step3. 重算官方 1-Jaccard 层次聚类和 79 个成员
official_distance <- as.dist(1 - jaccard_saved)
official_hc <- hclust(official_distance, method = "average")
official_k <- 2:(nrow(jaccard_saved) - 1L)
official_silhouette <- numeric(length(official_k))
for (i in seq_along(official_k)) {
  one_cluster <- cutree(official_hc, k = official_k[[i]])
  official_silhouette[[i]] <- summary(cluster::silhouette(one_cluster, official_distance))$avg.width
}
official_k_table <- data.frame(K = official_k, Silhouette = official_silhouette)
official_selected_k <- official_k[[which.max(official_silhouette)]]
official_initial_cluster <- cutree(official_hc, k = official_selected_k)
official_initial_label <- sprintf("IE%02d", official_initial_cluster)
names(official_initial_label) <- names(official_initial_cluster)

initial_saved <- read.delim(initial_ecotype_file, row.names = 1, check.names = FALSE)
initial_match <- identical(rownames(initial_saved), names(official_initial_label)) && identical(as.character(initial_saved$InitialEcotype), unname(official_initial_label))
if (!initial_match) stop("Initial official ecotype partition does not reproduce", call. = FALSE)

min_states <- 3L
official_sizes <- table(official_initial_label)
official_retained_initial <- names(official_sizes)[official_sizes >= min_states]
official_member_ids <- names(official_initial_label)[official_initial_label %in% official_retained_initial]
official_member_initial <- official_initial_label[official_member_ids]
official_ecotype_map <- paste0("E", seq_along(unique(official_member_initial)))
names(official_ecotype_map) <- unique(official_member_initial)
official_member_label <- unname(official_ecotype_map[official_member_initial])
names(official_member_label) <- official_member_ids

ecotypes <- read.delim(ecotype_file, check.names = FALSE)
saved_member_label <- setNames(as.character(ecotypes$Ecotype), ecotypes$ID)
member_identity_match <- setequal(names(saved_member_label), names(official_member_label))
member_partition_match <- member_identity_match && identical(saved_member_label[sort(names(saved_member_label))], official_member_label[sort(names(official_member_label))])
if (!member_partition_match) stop("Filtered official ecotype membership does not reproduce", call. = FALSE)

official_initial_silhouette <- max(official_silhouette)
official_member_distance <- as.dist(1 - jaccard_saved[names(saved_member_label), names(saved_member_label), drop = FALSE])
official_member_numeric <- as.integer(sub("E", "", saved_member_label))
official_filtered_silhouette <- summary(cluster::silhouette(official_member_numeric, official_member_distance))$avg.width
saved_initial_silhouette <- read.delim(file.path(ecotype_discovery_root, "silhouette_initial.txt"), check.names = FALSE)[[1L]][[1L]]
saved_filtered_silhouette <- read.delim(file.path(ecotype_discovery_root, "silhouette.txt"), check.names = FALSE)[[1L]][[1L]]
if (abs(official_initial_silhouette - saved_initial_silhouette) >= 1e-12 || abs(official_filtered_silhouette - saved_filtered_silhouette) >= 1e-12) stop("Saved official silhouette does not reproduce", call. = FALSE)
cat("OFFICIAL_PARTITION_REPRODUCED K", official_selected_k, "members", length(official_member_label), "ecotypes", length(unique(official_member_label)), "silhouette", signif(official_filtered_silhouette, 6), "\n")

# Step4. 重算 state/ecotype abundance 和正式 assignment
rank_data <- data.table::fread(rank_data_file, colClasses = "character", check.names = FALSE)
all_state_abundance <- NULL
assigned_member_samples <- character()
for (i in seq_len(nrow(rank_data))) {
  cell_type <- rank_data[[1L]][[i]]
  selected_rank <- trimws(rank_data[[2L]][[i]])
  state_dir <- file.path(state_discovery_root, cell_type, selected_rank)
  one_abundance <- read.delim(file.path(state_dir, "state_abundances.txt"), row.names = 1, check.names = FALSE)
  rownames(one_abundance) <- paste0(cell_type, "_", rownames(one_abundance))
  all_state_abundance <- rbind(all_state_abundance, one_abundance)

  one_assignment <- read.delim(file.path(state_dir, "state_assignment.txt"), check.names = FALSE)
  member_states <- ecotypes$State[ecotypes$CellType == cell_type]
  assigned_member_samples <- union(assigned_member_samples, one_assignment$ID[one_assignment$State %in% member_states])
}

all_state_abundance <- all_state_abundance[ecotypes$ID, , drop = FALSE]
saved_state_abundance <- as.matrix(read.delim(saved_state_abundance_file, row.names = 1, check.names = FALSE))
storage.mode(saved_state_abundance) <- "double"
state_abundance_error <- max(abs(as.matrix(all_state_abundance) - saved_state_abundance))
if (!identical(dimnames(as.matrix(all_state_abundance)), dimnames(saved_state_abundance)) || state_abundance_error >= 1e-12) stop("Combined state abundance does not reproduce", call. = FALSE)

ecotype_abundance <- do.call(rbind, lapply(unique(as.character(ecotypes$Ecotype)), function(one_ecotype) {
  colMeans(all_state_abundance[ecotypes$ID[ecotypes$Ecotype == one_ecotype], , drop = FALSE])
}))
rownames(ecotype_abundance) <- unique(as.character(ecotypes$Ecotype))
ecotype_abundance <- apply(ecotype_abundance, 2L, function(x) x / sum(x))
saved_ecotype_abundance <- as.matrix(read.delim(saved_ecotype_abundance_file, row.names = 1, check.names = FALSE))
storage.mode(saved_ecotype_abundance) <- "double"
ecotype_abundance_error <- max(abs(ecotype_abundance - saved_ecotype_abundance))
if (!identical(dimnames(ecotype_abundance), dimnames(saved_ecotype_abundance)) || ecotype_abundance_error >= 1e-12) stop("Ecotype abundance does not reproduce", call. = FALSE)

assignment_p <- matrix(1, nrow = nrow(ecotype_abundance), ncol = ncol(ecotype_abundance), dimnames = dimnames(ecotype_abundance))
for (i in seq_len(nrow(ecotype_abundance))) {
  member_rows <- ecotypes$ID[ecotypes$Ecotype == rownames(ecotype_abundance)[[i]]]
  for (j in seq_len(ncol(ecotype_abundance))) {
    assignment_p[i, j] <- tryCatch(t.test(all_state_abundance[member_rows, j], all_state_abundance[!rownames(all_state_abundance) %in% member_rows, j])$p.value, error = function(e) 1)
  }
}
saved_assignment_p <- as.matrix(read.delim(saved_assignment_p_file, row.names = 1, check.names = FALSE))
storage.mode(saved_assignment_p) <- "double"
assignment_p_error <- max(abs(assignment_p - saved_assignment_p))
if (!identical(dimnames(assignment_p), dimnames(saved_assignment_p)) || assignment_p_error >= 1e-10) stop("Assignment P values do not reproduce", call. = FALSE)

winner_index <- max.col(t(ecotype_abundance), ties.method = "first")
winner_ecotype <- rownames(ecotype_abundance)[winner_index]
winner_ties <- apply(ecotype_abundance, 2L, function(x) sum(x == max(x)))
winner_p <- assignment_p[cbind(winner_index, seq_len(ncol(assignment_p)))]
winner_q <- p.adjust(winner_p, method = "BH")
assigned_to_member_state <- colnames(ecotype_abundance) %in% assigned_member_samples
formal_assignment <- ifelse(winner_q < 0.25 & assigned_to_member_state, winner_ecotype, "Unassigned")

saved_initial_assignment <- read.delim(saved_initial_assignment_file, check.names = FALSE)
rownames(saved_initial_assignment) <- saved_initial_assignment$ID
saved_initial_assignment <- saved_initial_assignment[colnames(ecotype_abundance), , drop = FALSE]
assignment_match <- identical(as.character(saved_initial_assignment$MaxEcotype), winner_ecotype) &&
  max(abs(saved_initial_assignment$AssignmentP - winner_p)) < 1e-10 &&
  max(abs(saved_initial_assignment$AssignmentQ - winner_q)) < 1e-10 &&
  identical(as.logical(saved_initial_assignment$AssignedToEcotypeStates), assigned_to_member_state) &&
  identical(as.character(saved_initial_assignment$Ecotype), formal_assignment)
if (!assignment_match) stop("Formal ecotype assignment does not reproduce", call. = FALSE)
cat("ASSIGNMENT_REPRODUCED", sum(formal_assignment != "Unassigned"), "assigned", sum(formal_assignment == "Unassigned"), "unassigned\n")

# Step5. 计算 STAR Euclidean、Louvain 和 k-medoids 非阻断敏感性
star_distance <- dist(jaccard_saved, method = "euclidean")
star_hc <- hclust(star_distance, method = "average")
star_silhouette <- numeric(length(official_k))
for (i in seq_along(official_k)) {
  one_cluster <- cutree(star_hc, k = official_k[[i]])
  star_silhouette[[i]] <- summary(cluster::silhouette(one_cluster, star_distance))$avg.width
}
star_selected_k <- official_k[[which.max(star_silhouette)]]
star_cluster <- cutree(star_hc, k = star_selected_k)
star_sizes <- table(star_cluster)
star_member_ids <- names(star_cluster)[star_cluster %in% as.integer(names(star_sizes)[star_sizes >= min_states])]

pam_k <- 2:min(30L, nrow(jaccard_saved) - 1L)
pam_silhouette <- numeric(length(pam_k))
pam_clusters <- vector("list", length(pam_k))
for (i in seq_along(pam_k)) {
  pam_fit <- cluster::pam(official_distance, k = pam_k[[i]], diss = TRUE)
  pam_clusters[[i]] <- pam_fit$clustering
  pam_silhouette[[i]] <- pam_fit$silinfo$avg.width
}
pam_selected_index <- which.max(pam_silhouette)
pam_selected_k <- pam_k[[pam_selected_index]]
pam_cluster <- pam_clusters[[pam_selected_index]]
pam_sizes <- table(pam_cluster)
pam_member_ids <- names(pam_cluster)[pam_cluster %in% as.integer(names(pam_sizes)[pam_sizes >= min_states])]

networktoolbox_archive <- file.path(RUN_ROOT, "software_sources/NetworkToolbox_1.4.0.tar.gz")
networktoolbox_louvain <- file.path(RUN_ROOT, "software_sources/NetworkToolbox_1.4.0_louvain.R")
if (!file.exists(networktoolbox_archive) || !file.exists(networktoolbox_louvain)) stop("Missing frozen NetworkToolbox 1.4.0 source", call. = FALSE)
source(networktoolbox_louvain, local = FALSE)

louvain_gamma <- sort(unique(c(seq(0.01, 5, by = 0.01), seq(3.2012, 3.2030, by = 0.0001))))
louvain_cores <- 8L
louvain_clusters <- parallel::mclapply(louvain_gamma, function(one_gamma) louvain(jaccard_saved, gamma = one_gamma)$community, mc.cores = louvain_cores)
louvain_k <- vapply(louvain_clusters, function(x) length(unique(x)), integer(1))
louvain_evaluable <- louvain_k >= 2L & louvain_k <= 30L
if (!all(2:30 %in% louvain_k)) stop("Louvain gamma grid does not cover every K from 2 to 30", call. = FALSE)
louvain_silhouette <- rep(NA_real_, length(louvain_gamma))
for (i in which(louvain_evaluable)) louvain_silhouette[[i]] <- summary(cluster::silhouette(louvain_clusters[[i]], official_distance))$avg.width
louvain_selected_index <- which.max(louvain_silhouette)
louvain_selected_gamma <- louvain_gamma[[louvain_selected_index]]
louvain_selected_k <- louvain_k[[louvain_selected_index]]
louvain_cluster <- louvain_clusters[[louvain_selected_index]]
names(louvain_cluster) <- rownames(jaccard_saved)
louvain_sizes <- table(louvain_cluster)
louvain_member_ids <- names(louvain_cluster)[louvain_cluster %in% as.integer(names(louvain_sizes)[louvain_sizes >= min_states])]

common_star <- intersect(names(saved_member_label), star_member_ids)
star_contingency <- table(saved_member_label[common_star], star_cluster[common_star])
star_n <- sum(star_contingency)
star_index <- sum(star_contingency * (star_contingency - 1) / 2)
star_row_pairs <- sum(rowSums(star_contingency) * (rowSums(star_contingency) - 1) / 2)
star_col_pairs <- sum(colSums(star_contingency) * (colSums(star_contingency) - 1) / 2)
star_expected <- star_row_pairs * star_col_pairs / (star_n * (star_n - 1) / 2)
star_ari <- (star_index - star_expected) / ((star_row_pairs + star_col_pairs) / 2 - star_expected)

common_pam <- intersect(names(saved_member_label), pam_member_ids)
pam_contingency <- table(saved_member_label[common_pam], pam_cluster[common_pam])
pam_n <- sum(pam_contingency)
pam_index <- sum(pam_contingency * (pam_contingency - 1) / 2)
pam_row_pairs <- sum(rowSums(pam_contingency) * (rowSums(pam_contingency) - 1) / 2)
pam_col_pairs <- sum(colSums(pam_contingency) * (colSums(pam_contingency) - 1) / 2)
pam_expected <- pam_row_pairs * pam_col_pairs / (pam_n * (pam_n - 1) / 2)
pam_ari <- (pam_index - pam_expected) / ((pam_row_pairs + pam_col_pairs) / 2 - pam_expected)

common_louvain <- intersect(names(saved_member_label), louvain_member_ids)
louvain_contingency <- table(saved_member_label[common_louvain], louvain_cluster[common_louvain])
louvain_n <- sum(louvain_contingency)
louvain_index <- sum(louvain_contingency * (louvain_contingency - 1) / 2)
louvain_row_pairs <- sum(rowSums(louvain_contingency) * (rowSums(louvain_contingency) - 1) / 2)
louvain_col_pairs <- sum(colSums(louvain_contingency) * (colSums(louvain_contingency) - 1) / 2)
louvain_expected <- louvain_row_pairs * louvain_col_pairs / (louvain_n * (louvain_n - 1) / 2)
louvain_ari <- (louvain_index - louvain_expected) / ((louvain_row_pairs + louvain_col_pairs) / 2 - louvain_expected)

alternatives <- data.frame(
  Method = c("OFFICIAL_1_MINUS_JACCARD", "STAR_EUCLIDEAN_ROWS", "LOUVAIN_NETWORKTOOLBOX_1.4.0", "K_MEDOIDS_1_MINUS_JACCARD"),
  Status = c("PASS_EXACT_REPRODUCTION", "PASS_NONBLOCKING_SENSITIVITY", "PASS_NONBLOCKING_SENSITIVITY", "PASS_NONBLOCKING_SENSITIVITY"),
  SelectedK = c(official_selected_k, star_selected_k, louvain_selected_k, pam_selected_k),
  SelectedGamma = c(NA_real_, NA_real_, louvain_selected_gamma, NA_real_),
  RetainedStates = c(length(official_member_ids), length(star_member_ids), length(louvain_member_ids), length(pam_member_ids)),
  RetainedClusters = c(length(unique(saved_member_label)), length(unique(star_cluster[star_member_ids])), length(unique(louvain_cluster[louvain_member_ids])), length(unique(pam_cluster[pam_member_ids]))),
  MaxSilhouette = c(max(official_silhouette), max(star_silhouette), max(louvain_silhouette, na.rm = TRUE), max(pam_silhouette)),
  CommonOfficialStates = c(length(saved_member_label), length(common_star), length(common_louvain), length(common_pam)),
  ARIvsOfficial = c(1, star_ari, louvain_ari, pam_ari),
  stringsAsFactors = FALSE
)
print(alternatives)

# Step6. 保存审计结果，不写入正式模型目录
data.table::fwrite(official_k_table, file.path(audit_dir, "official_k_silhouette.tsv"), sep = "\t")
data.table::fwrite(data.frame(K = official_k, Silhouette = star_silhouette), file.path(audit_dir, "star_euclidean_k_silhouette.tsv"), sep = "\t")
data.table::fwrite(data.frame(Gamma = louvain_gamma, K = louvain_k, Silhouette = louvain_silhouette, Evaluable = louvain_evaluable), file.path(audit_dir, "louvain_gamma_silhouette.tsv"), sep = "\t")
data.table::fwrite(data.frame(K = pam_k, Silhouette = pam_silhouette), file.path(audit_dir, "kmedoids_k_silhouette.tsv"), sep = "\t")
data.table::fwrite(alternatives, file.path(audit_dir, "alternative_partition_summary.tsv"), sep = "\t")
data.table::fwrite(data.frame(ID = colnames(ecotype_abundance), MaxEcotype = winner_ecotype, AssignmentP = winner_p, AssignmentQ = winner_q, AssignedToEcotypeStates = assigned_to_member_state, Ecotype = formal_assignment), file.path(audit_dir, "assignment_recalculated.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "R_sessionInfo.txt"))

summary <- list(
  status = "PASS_ECOTYPE_AUDIT",
  candidate_states = nrow(jaccard_saved),
  official_selected_k = official_selected_k,
  official_retained_states = length(official_member_label),
  official_ecotypes = length(unique(official_member_label)),
  official_filtered_silhouette = official_filtered_silhouette,
  jaccard_max_absolute_error = jaccard_max_error,
  ecotype_abundance_max_absolute_error = ecotype_abundance_error,
  assignment_p_max_absolute_error = assignment_p_error,
  assigned_samples = sum(formal_assignment != "Unassigned"),
  unassigned_samples = sum(formal_assignment == "Unassigned"),
  abundance_winner_ties = sum(winner_ties > 1L),
  star_euclidean_selected_k = star_selected_k,
  star_euclidean_ari_vs_official = star_ari,
  louvain_package_version = "1.4.0",
  louvain_source_sha256 = digest::digest(networktoolbox_louvain, algo = "sha256", file = TRUE, serialize = FALSE),
  louvain_selected_gamma = louvain_selected_gamma,
  louvain_selected_k = louvain_selected_k,
  louvain_ari_vs_official = louvain_ari,
  kmedoids_selected_k = pam_selected_k,
  kmedoids_ari_vs_official = pam_ari,
  model_changed = FALSE,
  automatic_rerun = FALSE,
  release_decision = "G3_PASS_OFFICIAL_PRIMARY"
)
jsonlite::write_json(summary, file.path(audit_dir, "ecotype_discovery_audit.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
cat("ECOTYPE_AUDIT_PASS official model and alternative partitions complete; model unchanged\n")
