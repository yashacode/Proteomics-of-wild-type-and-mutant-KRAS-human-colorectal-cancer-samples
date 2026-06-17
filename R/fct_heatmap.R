#' Group-averaged, row z-scored expression heatmap
#'
#' Builds a heatmap from an expression matrix plus row/column annotations
#' (the components of a SummarizedExperiment, kept as plain objects so the
#' app has no Bioconductor runtime dependency). Selected features (genes)
#' are matched against a row-annotation column, samples are grouped by one
#' or more column-annotation columns, expression is averaged within each
#' group, and each row (gene) is z-scored across groups.
#'
#' @param expr Numeric expression matrix (features x samples), row-aligned
#'   to `row_data` and column-aligned to `col_data`.
#' @param row_data data.frame of feature annotations (one row per `expr` row).
#' @param col_data data.frame of sample annotations (one row per `expr` col).
#' @param genes Character vector of feature/gene names to display.
#' @param group_by Character vector of `col_data` column name(s) to group
#'   samples by. Multiple columns are combined into a single grouping.
#' @param gene_col `row_data` column to match `genes` against (default
#'   "Genes"). Values may contain multiple symbols separated by ";".
#' @param cluster_rows Order rows by hierarchical clustering (default TRUE).
#' @param interactive If TRUE, return a plotly object.
#' @param low_color,mid_color,high_color Diverging color scale endpoints.
#' @param plotly_source plotly source id when interactive.
#'
#' @return A ggplot (or plotly) object.
#' @noRd
plot_heatmap <- function(expr,
                         row_data,
                         col_data,
                         genes,
                         group_by,
                         gene_col = "Genes",
                         cluster_rows = TRUE,
                         interactive = FALSE,
                         low_color = "#0072B2",
                         mid_color = "white",
                         high_color = "#D55E00",
                         plotly_source = "heatmap",
                         click_input_id = NULL) {

  # ---- package checks ----
  required <- "ggplot2"
  if (interactive) required <- c(required, "plotly")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
  }

  # ---- input validation ----
  mat <- as.matrix(expr)
  rd <- as.data.frame(row_data)
  cd <- as.data.frame(col_data)
  if (nrow(mat) != nrow(rd)) stop("`expr` rows must align with `row_data`")
  if (ncol(mat) != nrow(cd)) stop("`expr` columns must align with `col_data`")

  if (!gene_col %in% colnames(rd)) {
    stop("gene_col '", gene_col, "' not found in row_data")
  }
  bad_groups <- group_by[!group_by %in% colnames(cd)]
  if (length(bad_groups) > 0) {
    stop("group_by column(s) not in col_data: ", paste(bad_groups, collapse = ", "))
  }

  genes <- unique(trimws(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) {
    stop("No genes provided")
  }

  # ---- match genes against rowData (handles ";"-separated symbols) ----
  row_symbols <- strsplit(as.character(rd[[gene_col]]), ";", fixed = TRUE)
  genes_lc <- tolower(genes)

  # For each requested gene, average matching rows into one sample-level vector
  per_gene <- lapply(genes, function(g) {
    hits <- which(vapply(row_symbols, function(s) tolower(g) %in% tolower(trimws(s)), logical(1)))
    if (length(hits) == 0) return(NULL)
    colMeans(mat[hits, , drop = FALSE], na.rm = TRUE)
  })
  matched <- !vapply(per_gene, is.null, logical(1))
  if (!any(matched)) {
    stop("None of the requested genes were found in rowData$", gene_col)
  }
  gene_mat <- do.call(rbind, per_gene[matched])
  rownames(gene_mat) <- genes[matched]

  # ---- group samples and average ----
  group_factor <- if (length(group_by) == 1) {
    as.character(cd[[group_by]])
  } else {
    apply(cd[, group_by, drop = FALSE], 1, paste, collapse = " | ")
  }
  groups <- sort(unique(group_factor))

  group_mat <- vapply(groups, function(g) {
    rowMeans(gene_mat[, group_factor == g, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(gene_mat)))
  if (is.null(dim(group_mat))) {
    group_mat <- matrix(group_mat, nrow = nrow(gene_mat),
                        dimnames = list(rownames(gene_mat), groups))
  }

  # ---- z-score each row across groups ----
  z <- t(apply(group_mat, 1, function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(x - mean(x, na.rm = TRUE))
    (x - mean(x, na.rm = TRUE)) / s
  }))
  colnames(z) <- colnames(group_mat)

  # ---- row ordering ----
  row_order <- rownames(z)
  if (cluster_rows && nrow(z) > 2) {
    row_order <- rownames(z)[stats::hclust(stats::dist(z))$order]
  }

  # ---- long form for ggplot ----
  plot_df <- data.frame(
    gene = factor(rep(rownames(z), times = ncol(z)), levels = row_order),
    group = factor(rep(colnames(z), each = nrow(z)), levels = colnames(z)),
    z = as.vector(z),
    stringsAsFactors = FALSE
  )
  plot_df$tooltip <- paste0(
    "Gene: ", plot_df$gene,
    "<br>Group: ", plot_df$group,
    "<br>z-score: ", round(plot_df$z, 2)
  )

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = group, y = gene, fill = z, text = tooltip)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = low_color, mid = mid_color, high = high_color,
      midpoint = 0, name = "Row z-score"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    )

  if (interactive) {
    p <- plotly::ggplotly(p, tooltip = "text", source = plotly_source)

    # Heatmap clicks carry no `key`; push the clicked row position (y) to a
    # dedicated Shiny input with a nonce + priority "event" so re-clicking
    # the same cell still fires. The caller maps y -> gene via gene_order.
    if (!is.null(click_input_id)) {
      p <- htmlwidgets::onRender(p, sprintf(
        "function(el){
           if (el.removeAllListeners) el.removeAllListeners('plotly_click');
           el.on('plotly_click', function(d){
             var pt = d.points[0];
             if (pt && pt.y !== undefined) {
               Shiny.setInputValue('%s', {value: pt.y, nonce: Math.random()}, {priority: 'event'});
             }
           });
         }",
        click_input_id
      ))
    }
  }

  # Heatmap-trace clicks don't carry a `key`, so expose the y-axis order
  # (bottom -> top) for callers to resolve a clicked row index to a gene.
  attr(p, "gene_order") <- row_order

  p
}
