plot_volcano <- function(
    df,

    fc_col = "logFC",
    p_col = "padj",
    feature_col = "Genes",

    fc_cutoff = 1,
    p_cutoff = 0.05,

    label_features = NULL,
    max_labels = 15,

    up_color = "#D55E00",
    down_color = "#0072B2",
    ns_color = "grey80",

    up_label = "Up",
    down_label = "Down",
    ns_label = "NS",

    point_size = 2,
    alpha = 0.8,

    tooltip_cols = NULL,

    interactive = FALSE,
    plotly_source = "volcano",
    click_input_id = NULL,

    title = NULL,
    y_axis_label = NULL
) {

  # =========================
  # package checks
  # =========================

  required_packages <- c(
    "ggplot2",
    "dplyr",
    "ggrepel"
  )

  if (interactive) {
    required_packages <- c(required_packages, "plotly")
  }

  missing_pkgs <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_pkgs) > 0) {

    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", ")
    )
  }

  # =========================
  # input validation
  # =========================

  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame")
  }

  required_cols <- c(
    fc_col,
    p_col,
    feature_col
  )

  missing_cols <- required_cols[
    !required_cols %in% colnames(df)
  ]

  if (length(missing_cols) > 0) {

    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(df[[fc_col]])) {
    stop(fc_col, " must be numeric")
  }

  if (!is.numeric(df[[p_col]])) {
    stop(p_col, " must be numeric")
  }

  if (any(df[[p_col]] <= 0, na.rm = TRUE)) {

    stop(
      "P-values must be > 0 for -log10 transformation"
    )
  }

  if (!is.numeric(fc_cutoff) ||
      length(fc_cutoff) != 1) {

    stop(
      "fc_cutoff must be a single numeric value"
    )
  }

  if (!is.numeric(p_cutoff) ||
      length(p_cutoff) != 1) {

    stop(
      "p_cutoff must be a single numeric value"
    )
  }

  if (!is.null(label_features) &&
      !is.character(label_features)) {

    stop(
      "label_features must be NULL or character vector"
    )
  }

  if (!is.logical(interactive) ||
      length(interactive) != 1) {

    stop(
      "interactive must be TRUE or FALSE"
    )
  }

  # =========================
  # tooltip validation
  # =========================

  if (!is.null(tooltip_cols)) {

    bad_tooltip_cols <- tooltip_cols[
      !tooltip_cols %in% colnames(df)
    ]

    if (length(bad_tooltip_cols) > 0) {

      stop(
        "tooltip_cols missing from df: ",
        paste(
          bad_tooltip_cols,
          collapse = ", "
        )
      )
    }
  }

  # =========================
  # tooltip construction
  # =========================

  tooltip_text <- paste0(
    feature_col, ": ",
    df[[feature_col]],

    "<br>",
    fc_col, ": ",
    round(df[[fc_col]], 3),

    "<br>",
    p_col, ": ",
    signif(df[[p_col]], 3)
  )

  if (!is.null(tooltip_cols)) {

    for (col in tooltip_cols) {

      tooltip_text <- paste0(
        tooltip_text,

        "<br>",
        col, ": ",
        df[[col]]
      )
    }
  }

  # =========================
  # preprocessing
  # =========================

  plot_df <- df |>
    dplyr::mutate(

      .fc = .data[[fc_col]],
      .p = .data[[p_col]],
      .feature = .data[[feature_col]],

      .neglog10p = -log10(.p),

      .tooltip = tooltip_text,

      regulation = dplyr::case_when(

        .p < p_cutoff &
          .fc >= fc_cutoff ~ up_label,

        .p < p_cutoff &
          .fc <= -fc_cutoff ~ down_label,

        TRUE ~ ns_label
      )
    )

  # Split into a non-interactive background layer (NS) and an
  # interactive layer (significant). Only the significant points
  # carry tooltip/key aesthetics, so plotly skips hover-text
  # generation for the bulk of points -> much faster conversion.
  ns_df  <- plot_df |> dplyr::filter(regulation == ns_label)
  sig_df <- plot_df |> dplyr::filter(regulation != ns_label)

  # =========================
  # label selection
  # =========================

  if (is.null(label_features)) {

    label_df <- plot_df |>
      dplyr::filter(
        regulation != ns_label
      ) |>
      dplyr::arrange(.p) |>
      head(max_labels)

  } else {

    label_df <- plot_df |>
      dplyr::filter(
        .feature %in% label_features
      )
  }

  # =========================
  # ggplot construction
  # =========================

  p <- ggplot2::ggplot(

    mapping = ggplot2::aes(
      x = .fc,
      y = .neglog10p
    )

  ) +

    # Background NS points: no tooltip/key -> not interactive, fast
    ggplot2::geom_point(

      data = ns_df,

      color = ns_color,
      size = point_size,
      alpha = alpha
    ) +

    # Significant points: interactive (hover + click via key)
    ggplot2::geom_point(

      data = sig_df,

      ggplot2::aes(
        color = regulation,
        text = .tooltip,
        key = .feature
      ),

      size = point_size,
      alpha = alpha
    ) +

    ggplot2::scale_color_manual(

      values = stats::setNames(
        c(up_color, down_color),
        c(up_label, down_label)
      )
    ) +

    ggplot2::geom_vline(

      xintercept = c(
        -fc_cutoff,
        fc_cutoff
      ),

      linetype = "dashed",
      color = "grey50"
    ) +

    ggplot2::geom_hline(

      yintercept = -log10(p_cutoff),

      linetype = "dashed",
      color = "grey50"
    ) +

    ggrepel::geom_text_repel(

      data = label_df,

      ggplot2::aes(
        label = .feature
      ),

      size = 3,
      max.overlaps = Inf
    ) +

    ggplot2::labs(

      title = title,

      x = "log2 Fold Change",

      # Use a plain string (an expression() y label is dropped by
      # plotly during ggplotly conversion).
      y = if (is.null(y_axis_label)) "-log10(p)" else y_axis_label,

      color = NULL
    ) +

    ggplot2::theme_minimal()

  # =========================
  # optional plotly conversion
  # =========================

  if (interactive) {

    p <- plotly::ggplotly(
      p,
      tooltip = "text",
      source = plotly_source
    )

    # Render points on the GPU via WebGL -> much faster for large
    # scatter plots (and lighter in the browser).
    p <- plotly::toWebGL(p)

    # Push each click to a dedicated Shiny input with a nonce + priority
    # "event", so re-clicking the same point still fires (plotly's own
    # plotly_click input dedupes identical consecutive values).
    if (!is.null(click_input_id)) {
      p <- htmlwidgets::onRender(p, sprintf(
        "function(el){
           if (el.removeAllListeners) el.removeAllListeners('plotly_click');
           el.on('plotly_click', function(d){
             var pt = d.points[0];
             if (!pt) return;
             // plotly stores the `key` aesthetic in the trace's data.key
             // array (indexed by pointNumber), not directly on the point.
             var k = null;
             if (pt.key !== undefined && pt.key !== null) {
               k = pt.key;
             } else if (pt.data && pt.data.key) {
               k = pt.data.key[pt.pointNumber];
             }
             if (k !== null && k !== undefined) {
               Shiny.setInputValue('%s', {value: k, nonce: Math.random()}, {priority: 'event'});
             }
           });
         }",
        click_input_id
      ))
    }
  }

  return(p)
}
