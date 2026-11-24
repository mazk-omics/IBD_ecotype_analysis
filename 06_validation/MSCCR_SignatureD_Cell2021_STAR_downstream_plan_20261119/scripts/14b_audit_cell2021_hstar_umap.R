#!/usr/bin/env Rscript
# 本脚本审计 H* UMAP 坐标、正式 assignment 和 H* 最大分量的一致性。
# 只读取已完成坐标与 manifest，不重跑 UMAP。
options(stringsAsFactors = FALSE, warn = 1)

# Step1. 读取已完成输出
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 14b_audit_cell2021_hstar_umap.R RUN_ROOT", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(RUN_ROOT, "config/run_paths.R"))

if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Missing audit packages", call. = FALSE)
library(data.table)

result_dir <- file.path(output_root, "10_state_annotation/hstar_umap_cell2021_exact")
coordinates <- fread(file.path(result_dir, "hstar_umap_coordinates.tsv.gz"))
manifest <- fread(file.path(result_dir, "hstar_umap_manifest.tsv"))
completion <- jsonlite::read_json(file.path(result_dir, "hstar_umap_complete.json"), simplifyVector = TRUE)
if (!identical(completion$status, "PASS_CELL2021_HSTAR_UMAP_EXACT")) stop("H* UMAP completion status failed", call. = FALSE)
if (nrow(coordinates) != sum(manifest$UmapSamples) || nrow(manifest) != expected_cell_types || any(manifest$UmapSamples != manifest$FormalAssigned)) stop("H* UMAP dimensions failed", call. = FALSE)
if (any(!is.finite(coordinates$UMAP1)) || any(!is.finite(coordinates$UMAP2))) stop("H* UMAP contains non-finite coordinates", call. = FALSE)

# Step2. 核对正式 assignment 与 H* 最大分量
assignment_audit <- coordinates[FormalAssigned == TRUE, .(
  Assigned = .N, Match = sum(FormalState == HArgmaxState), Mismatch = sum(FormalState != HArgmaxState), MatchRate = mean(FormalState == HArgmaxState)
), by = CellType]
assignment_audit <- merge(manifest[, .(CellType, States, AllSamples, UmapSamples, FormalAssignedManifest = FormalAssigned)], assignment_audit, by = "CellType", all = TRUE, sort = FALSE)
if (anyNA(assignment_audit) || any(assignment_audit$Assigned != assignment_audit$FormalAssignedManifest) || any(assignment_audit$Mismatch != 0L)) stop("Formal assignment does not match H* argmax", call. = FALSE)
fwrite(assignment_audit, file.path(result_dir, "hstar_umap_assignment_identity_audit.tsv"), sep = "\t", quote = FALSE)

audit_completion <- list(
  status = "PASS_CELL2021_HSTAR_UMAP_EXACT_AUDIT", generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  cell_types = nrow(assignment_audit), coordinates = nrow(coordinates), formal_assignments = sum(assignment_audit$Assigned), mismatches = sum(assignment_audit$Mismatch)
)
jsonlite::write_json(audit_completion, file.path(result_dir, "hstar_umap_audit_complete.json"), pretty = TRUE, auto_unbox = TRUE)

print(assignment_audit)
cat("PASS_CELL2021_HSTAR_UMAP_EXACT_AUDIT\n")
