# 当前服务器 grid 版本把 patchwork 的 null panel units 转成 0；本文件只替换 layout 测量步骤。
# 先 source 官方 panel_alignment.R，再 source 本文件；审计规则和严格阈值保持不变。

.resolve_patchwork_units_pt <- function(units, total_pt, axis_name) {
  unit_types <- vapply(seq_along(units), function(index) grid::unitType(units[index]), character(1L))
  converted <- if (axis_name == "width") grid::convertWidth(units, "pt", valueOnly = TRUE) else grid::convertHeight(units, "pt", valueOnly = TRUE)
  null_units <- unit_types == "null"
  if (!any(null_units)) return(converted)

  null_weights <- vapply(which(null_units), function(index) as.numeric(units[index]), numeric(1L))
  available_pt <- total_pt - sum(converted[!null_units])
  if (!is.finite(available_pt) || available_pt <= 0 || sum(null_weights) <= 0) stop("Could not resolve patchwork null units", call. = FALSE)
  converted[null_units] <- available_pt * null_weights / sum(null_weights)
  converted
}

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
  if (!is.numeric(width_in) || length(width_in) != 1L || !is.finite(width_in) || width_in <= 0 || !is.numeric(height_in) || length(height_in) != 1L || !is.finite(height_in) || height_in <= 0) stop("width_in and height_in must be positive finite numbers", call. = FALSE)

  grob <- patchwork::patchworkGrob(plot)
  panel_rows <- .nature_alignment_panel_rows(grob)
  if (is.null(panel_ids)) panel_ids <- vapply(seq_len(nrow(panel_rows)), function(index) if (index <= length(letters)) letters[index] else paste0("panel-", index), character(1L))
  panel_ids <- as.character(panel_ids)
  if (length(panel_ids) != nrow(panel_rows) || anyDuplicated(panel_ids) || any(!nzchar(panel_ids))) stop("panel_ids must be unique, non-empty, and match the measured patchwork panels", call. = FALSE)

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

  widths_pt <- .resolve_patchwork_units_pt(grob$widths, width_in * 72, "width")
  heights_pt <- .resolve_patchwork_units_pt(grob$heights, height_in * 72, "height")
  if (any(!is.finite(widths_pt)) || any(!is.finite(heights_pt))) stop("Patchwork panel dimensions could not be converted to physical points", call. = FALSE)
  x_edges <- c(0, cumsum(widths_pt))
  top_edges <- c(0, cumsum(heights_pt))
  total_height_pt <- sum(heights_pt)

  panels <- lapply(seq_len(nrow(panel_rows)), function(index) {
    row <- panel_rows[index, ]
    left <- x_edges[row$l]
    right <- x_edges[row$r + 1L]
    top <- total_height_pt - top_edges[row$t]
    bottom <- total_height_pt - top_edges[row$b + 1L]
    list(
      id = panel_ids[index], bbox_pt = unname(c(left, bottom, right, top)), grid_id = "patchwork-grid-1",
      row_start = as.integer(row$t - 1L), row_stop = as.integer(row$b), col_start = as.integer(row$l - 1L), col_stop = as.integer(row$r)
    )
  })

  manifest <- list(schema_version = 1L, backend = "r-patchwork", figure = list(width_pt = width_in * 72, height_pt = height_in * 72), panels = panels, exemptions = exemptions)
  normalized_rows <- .nature_alignment_groups(row_groups)
  normalized_columns <- .nature_alignment_groups(column_groups)
  if (!is.null(normalized_rows)) manifest$row_groups <- normalized_rows
  if (!is.null(normalized_columns)) manifest$column_groups <- normalized_columns

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
