# 本兼容层仅扩展新版 patchwork 的 panel 名称识别，不改变测量或审计标准。
.nature_alignment_panel_rows <- function(grob) {
  panel_name <- grob$layout$name
  is_panel <- grepl("^panel(?:-[0-9]+)?$", panel_name, perl = TRUE) |
    grepl("^panel;.*-[0-9]+$", panel_name, perl = TRUE)
  rows <- grob$layout[is_panel, , drop = FALSE]
  if (nrow(rows) < 2L) stop("At least two patchwork panel cells are required for alignment QA", call. = FALSE)
  rows[order(rows$t, rows$l), , drop = FALSE]
}

# 新版 patchwork 的 null 单位不能从 grob$widths 直接换算，改为读取绘制后的 panel viewport。
write_patchwork_panel_layout <- function(
  plot,
  manifest_path,
  width_in,
  height_in,
  panel_ids = NULL,
  row_groups = NULL,
  column_groups = NULL,
  exemptions = list()
) {
  for (package in c("patchwork", "grid", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) stop(sprintf("Package '%s' is required for patchwork alignment QA", package), call. = FALSE)
  }

  grob <- patchwork::patchworkGrob(plot)
  probe_path <- tempfile(fileext = ".pdf")
  grDevices::pdf(probe_path, width = width_in, height = height_in, useDingbats = FALSE)
  device_open <- TRUE
  on.exit({
    if (device_open) grDevices::dev.off()
    unlink(probe_path)
  }, add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(grob)
  grid::grid.force()

  viewport_index <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  viewport_names <- grep("^(panel;|panel-[0-9]+[.])", viewport_index$name, value = TRUE)
  if (length(viewport_names) < 2L) stop("At least two patchwork panel viewports are required for alignment QA", call. = FALSE)

  panel_boxes <- vector("list", length(viewport_names))
  for (panel_index in seq_along(viewport_names)) {
    grid::seekViewport(viewport_names[[panel_index]])
    lower <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    upper <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    grid::upViewport(0)
    panel_boxes[[panel_index]] <- c(left = lower$x * 72, bottom = lower$y * 72, right = upper$x * 72, top = upper$y * 72)
  }
  panel_boxes <- do.call(rbind, panel_boxes)
  panel_boxes <- panel_boxes[order(-panel_boxes[, "bottom"], panel_boxes[, "left"]), , drop = FALSE]

  if (is.null(panel_ids)) panel_ids <- sprintf("panel-%02d", seq_len(nrow(panel_boxes)))
  panel_ids <- as.character(panel_ids)
  if (length(panel_ids) != nrow(panel_boxes) || anyDuplicated(panel_ids) || any(!nzchar(panel_ids))) stop("panel_ids must be unique, non-empty, and match the measured patchwork panels", call. = FALSE)

  row_levels <- sort(unique(round(panel_boxes[, "bottom"], 6)), decreasing = TRUE)
  column_levels <- sort(unique(round(panel_boxes[, "left"], 6)))
  panels <- lapply(seq_len(nrow(panel_boxes)), function(panel_index) {
    row_id <- match(round(panel_boxes[panel_index, "bottom"], 6), row_levels) - 1L
    column_id <- match(round(panel_boxes[panel_index, "left"], 6), column_levels) - 1L
    list(
      id = panel_ids[[panel_index]], bbox_pt = unname(panel_boxes[panel_index, ]), grid_id = "patchwork-grid-1",
      row_start = row_id, row_stop = row_id + 1L, col_start = column_id, col_stop = column_id + 1L
    )
  })

  manifest <- list(
    schema_version = 1L, backend = "r-patchwork",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = panels, exemptions = exemptions
  )
  if (!is.null(row_groups)) manifest$row_groups <- unname(.nature_alignment_groups(row_groups))
  if (!is.null(column_groups)) manifest$column_groups <- unname(.nature_alignment_groups(column_groups))

  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(pattern = ".panel-layout-", tmpdir = dirname(manifest_path), fileext = ".json")
  jsonlite::write_json(manifest, temporary_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  if (!file.rename(temporary_path, manifest_path)) {
    unlink(temporary_path)
    stop("Could not activate the patchwork alignment manifest atomically", call. = FALSE)
  }

  grDevices::dev.off()
  device_open <- FALSE
  unlink(probe_path)
  invisible(manifest)
}
