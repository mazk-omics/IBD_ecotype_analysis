#!/usr/bin/env Rscript
# 本脚本按 Cell 2021 方法对每个母细胞类型的正式 H* abundance 做 Euclidean UMAP。
# UMAP 只用于显示冻结 state，不参与重新聚类、state 命名或 validation 判定。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 固定环境、输入和随机种子
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 14_plot_cell2021_hstar_umap.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

umap_lib <- file.path(project_root, "software/Rlibs/cell2021_umap_R403_0.2.0.0")
if (!dir.exists(umap_lib)) stop("Cell 2021 umap library is missing", call. = FALSE)
.libPaths(c(umap_lib, .libPaths()))

required_pkgs <- c("data.table", "digest", "jsonlite", "umap", "ggplot2")
pkg_ok <- vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(pkg_ok)) stop("Missing R packages: ", paste(required_pkgs[!pkg_ok], collapse = ", "), call. = FALSE)
if (as.character(packageVersion("umap")) != "0.2.0.0") stop("umap version is not 0.2.0.0", call. = FALSE)
library(data.table)
library(ggplot2)

result_dir <- file.path(output_root, "10_state_annotation/hstar_umap_cell2021_exact")
figure_dir <- file.path(result_dir, "individual_cell_types")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
completion_file <- file.path(result_dir, "hstar_umap_complete.json")
if (file.exists(completion_file)) stop("H* UMAP already complete; do not overwrite", call. = FALSE)

umap_seed <- 20261122L
rank_data <- fread(rank_data_file, colClasses = "character", check.names = FALSE)
rank_lookup <- setNames(trimws(rank_data[[2L]]), rank_data[[1L]])
annotation <- fread(annotation_file, colClasses = "character", check.names = FALSE)
if (nrow(annotation) != expected_samples || anyDuplicated(annotation$ID)) stop("Formal annotation identity failed", call. = FALSE)

palette <- c(S01 = "#4477AA", S02 = "#EE6677", S03 = "#228833", S04 = "#CCBB44", S05 = "#66CCEE", S06 = "#AA3377")
coordinate_parts <- list()
manifest_parts <- list()

# Step2. 每个母细胞类型独立运行相同默认参数 UMAP
for (cell_index in seq_along(rank_lookup)) {
  cell_type <- names(rank_lookup)[[cell_index]]
  rank <- rank_lookup[[cell_type]]
  type_dir <- file.path(state_discovery_root, cell_type, rank)
  abundance_file <- file.path(type_dir, "state_abundances.txt")
  assignment_file <- file.path(type_dir, "state_assignment.txt")
  if (!file.exists(abundance_file) || !file.exists(assignment_file)) stop("Missing H* UMAP input: ", cell_type, call. = FALSE)

  abundance <- fread(abundance_file, check.names = FALSE)
  state_ids <- abundance[[1L]]
  abundance <- as.matrix(abundance[, -1L, with = FALSE])
  rownames(abundance) <- state_ids
  storage.mode(abundance) <- "double"
  if (nrow(abundance) < 2L || ncol(abundance) != expected_samples || any(!is.finite(abundance)) || any(abundance < 0)) stop("Invalid H* abundance: ", cell_type, call. = FALSE)
  if (anyDuplicated(colnames(abundance)) || !setequal(colnames(abundance), annotation$ID)) stop("H* sample identity failed: ", cell_type, call. = FALSE)

  assignment <- fread(assignment_file, header = FALSE, skip = 1L, colClasses = "character")
  if (ncol(assignment) != 4L) stop("Unexpected formal assignment layout: ", cell_type, call. = FALSE)
  setnames(assignment, c("RowIndex", "SampleID", "State", "InitialState"))
  if (anyDuplicated(assignment$SampleID) || any(!assignment$SampleID %in% colnames(abundance)) || any(!assignment$State %in% state_ids)) stop("Formal assignment identity failed: ", cell_type, call. = FALSE)

  hstar_all <- t(abundance)
  hstar <- hstar_all[assignment$SampleID, , drop = FALSE]
  if (!identical(rownames(hstar), assignment$SampleID) || nrow(hstar) != nrow(assignment)) stop("Assigned-only H* identity failed: ", cell_type, call. = FALSE)
  distance <- as.matrix(stats::dist(hstar, method = "euclidean"))
  if (!identical(rownames(distance), rownames(hstar)) || any(!is.finite(distance))) stop("Euclidean distance identity failed: ", cell_type, call. = FALSE)
  set.seed(umap_seed)
  embedding <- umap::umap(distance, input = "dist", random_state = umap_seed)$layout
  if (!identical(rownames(embedding), rownames(hstar)) || ncol(embedding) != 2L || any(!is.finite(embedding))) stop("Invalid UMAP coordinates: ", cell_type, call. = FALSE)

  coordinates <- data.table(SampleID = rownames(embedding), UMAP1 = embedding[, 1L], UMAP2 = embedding[, 2L])
  coordinates[, `:=`(CellType = cell_type, Rank = rank, HArgmaxState = state_ids[max.col(hstar, ties.method = "first")])]
  coordinates <- merge(coordinates, assignment[, .(SampleID, FormalState = State)], by = "SampleID", all.x = TRUE, sort = FALSE)
  coordinates <- merge(coordinates, annotation[, .(SampleID = ID, PatientID, Disease, Inflammation, Region3)], by = "SampleID", all.x = TRUE, sort = FALSE)
  if (anyNA(coordinates$FormalState) || any(coordinates$FormalState != coordinates$HArgmaxState)) stop("Assigned-only state identity failed: ", cell_type, call. = FALSE)
  coordinates[, `:=`(FormalAssigned = TRUE, DisplayState = factor(FormalState, levels = state_ids))]
  coordinate_parts[[cell_type]] <- coordinates

  manifest_parts[[cell_type]] <- data.table(
    CellType = cell_type, Rank = rank, States = nrow(abundance), AllSamples = ncol(abundance), UmapSamples = nrow(hstar), FormalAssigned = nrow(assignment),
    Seed = umap_seed, Metric = "precomputed_euclidean", UmapInput = "dist", RVersion = R.version.string, UmapVersion = as.character(packageVersion("umap")),
    AbundanceSHA256 = digest::digest(abundance_file, file = TRUE, algo = "sha256"), AssignmentSHA256 = digest::digest(assignment_file, file = TRUE, algo = "sha256")
  )

  plot_type <- ggplot(coordinates, aes(UMAP1, UMAP2, colour = DisplayState)) +
    geom_point(size = 0.55, alpha = 0.70) +
    scale_colour_manual(values = palette, drop = FALSE) +
    coord_equal() +
    labs(title = cell_type, subtitle = paste0(nrow(assignment), " formally assigned biopsies"), colour = "State") +
    theme_classic(base_size = 11) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "right")
  ggsave(file.path(figure_dir, paste0(cell_type, "_Hstar_UMAP.png")), plot_type, width = 6.2, height = 5.2, dpi = 320, bg = "white")
  ggsave(file.path(figure_dir, paste0(cell_type, "_Hstar_UMAP.pdf")), plot_type, width = 6.2, height = 5.2, device = cairo_pdf)
  cat("HSTAR_UMAP", cell_type, nrow(abundance), "states", nrow(assignment), "assigned\n")
}

# Step3. 保存坐标和 20 个母细胞类型总览图
coordinates_all <- rbindlist(coordinate_parts, use.names = TRUE)
manifest <- rbindlist(manifest_parts)
if (nrow(coordinates_all) != sum(manifest$UmapSamples) || nrow(manifest) != expected_cell_types || any(manifest$UmapSamples != manifest$FormalAssigned)) stop("H* UMAP output dimensions failed", call. = FALSE)
fwrite(coordinates_all, file.path(result_dir, "hstar_umap_coordinates.tsv.gz"), sep = "\t", quote = FALSE)
fwrite(manifest, file.path(result_dir, "hstar_umap_manifest.tsv"), sep = "\t", quote = FALSE)

overview <- ggplot(coordinates_all, aes(UMAP1, UMAP2, colour = DisplayState)) +
  geom_point(size = 0.16, alpha = 0.55) +
  scale_colour_manual(values = palette, drop = FALSE) +
  facet_wrap(~CellType, scales = "free", ncol = 4) +
  labs(title = expression("Cell 2021 " * H^"*" * " abundance UMAP"), subtitle = "Assigned biopsies only; precomputed Euclidean distance; umap 0.2.0.0 defaults", colour = "State") +
  theme_classic(base_size = 9) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), strip.text = element_text(size = 8), legend.position = "bottom")
ggsave(file.path(result_dir, "all_cell_types_Hstar_UMAP.png"), overview, width = 14, height = 17, dpi = 320, bg = "white")
ggsave(file.path(result_dir, "all_cell_types_Hstar_UMAP.pdf"), overview, width = 14, height = 17, device = cairo_pdf)

completion <- list(
  status = "PASS_CELL2021_HSTAR_UMAP_EXACT", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  cell_types = nrow(manifest), all_samples_per_cell_type = expected_samples, assigned_coordinates = nrow(coordinates_all),
  r_version = R.version.string, umap_version = as.character(packageVersion("umap")), seed = umap_seed, metric = "precomputed_euclidean", umap_input = "dist", configuration = "umap.defaults"
)
jsonlite::write_json(completion, completion_file, pretty = TRUE, auto_unbox = TRUE)

print(manifest[, .(CellType, States, FormalAssigned)])
cat("PASS_CELL2021_HSTAR_UMAP_EXACT\n")
